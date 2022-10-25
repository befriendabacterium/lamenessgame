library(scales)
library(dplyr)

# READ IN RAW DATA FROM MS FORMS --------------------------------------------------------

studydata_raw<-readxl::read_excel('inputs/lamenessstudydata_010721.xlsx')
#detach(tibble)
newcolnames<-readRDS('inputs/newcolnames.RDS')

# FORMATTING: LIKERT DATA --------------------------------------------------------------

# subset out and tidy likert columns
firstcol<-which(colnames(studydata_raw)=="Points - How strongly do you agree with the following statements?")
#cut out section of df where the likert data starts
likertdata<-studydata_raw[,firstcol:ncol(studydata_raw)]
#remove unwanted columns
likertdata<-likertdata[,-grep("Points|Feedback|Please share any general feedback/main thoughts after playing the game below|Finally, how did you find out about this study?",colnames(likertdata))]
#remove weird bits from some of the statements
colnames(likertdata)<-gsub("\r\n\r\n","",colnames(likertdata))
#convert columns to factor
likertdata<-as.data.frame(apply(likertdata,2,factor))

#count up the responses
likertdata_formatted<-apply(likertdata,2,function(x){plyr::count(x)})
likertdata_formatted

#rename columns with loop
for (l in 1:length(likertdata_formatted)){
  colnames(likertdata_formatted[[l]])[1]<-"Statement"
  colnames(likertdata_formatted[[l]])[2]<-names(likertdata_formatted[l])
}

#function for merge
my_merge <- function(df1, df2){                                # Create own merging function
  merge(df1, df2, by = "Statement", all=T)
}

#apply function to list
likertdata_formatted<-Reduce(my_merge, likertdata_formatted)                                    # Apply Reduce to own function

#transpose dataframe
likertdata_formatted<-t(likertdata_formatted)

#set colnames
colnames(likertdata_formatted)<-likertdata_formatted[1,]

#remove first row (colnames)
likertdata_formatted<-as.data.frame(likertdata_formatted[-1,])

#remove NA col (did not answer question cos non-farmer or missed it)
likertdata_formatted<-likertdata_formatted[,!is.na(colnames(likertdata_formatted))]

#vector of likert categories
likert_categories_ordered<-c("Strongly disagree", "Disagree", "Neutral", "Agree", "Strongly agree")

#re-order columns
likertdata_formatted<-likertdata_formatted[,match(likert_categories_ordered,colnames(likertdata_formatted))]

likertdata_formatted

write.csv(likertdata_formatted, 'outputs/processed_data/likertdata_formatted.csv')

# FORMATTING: PARTICIPANT PERFORMANCE DATA --------------------------------------------------------------

#REDUNDANT COLUMN REMOVAL (EASIER TO WORK WITH)
#make a vector of columns to remove
colstoremove<-grep(c('Email|Name|Total.points|Quiz.feedback|Points...|Feedback...'), colnames(studydata_raw))
#remove them
studydata_formatted<-studydata_raw[,-colstoremove]

#RENAMING OF COLUMNS (EASIER TO WORK WITH)
#rename columns
colnames(studydata_formatted)<-newcolnames

#SETTING CLASSES OF COLUMNS
#check classes of columns
apply(studydata_formatted,1,class)
#make numeric columns numeric
studydata_formatted[,c(5,6,7,8,20)]<-sapply(studydata_formatted[,c(5,6,7,8,20)],as.numeric)
#recheck classes of columns
apply(studydata_formatted,1,class)

#remove participant ID36 as game didn't work properly for them
studydata_formatted<-studydata_formatted[-which(studydata_formatted$ID==36),]

# PARSE FARMING EXPERIENCE  ---------------------------------------------------------------------

farmingexp_df<-studydata_formatted %>% select(farmingexp_YN, sheeproles_type, sheeplamenessexp_prev,sheepexp_years)

#change column classes
farmingexp_df<-farmingexp_df %>%
  mutate(across(c(farmingexp_YN,
                  sheeproles_type, 
                  sheeplamenessexp_prev),
                factor))


#combine sheep roles type cos poor sampling at tails
levels(farmingexp_df$sheeproles_type)

levels(farmingexp_df$sheeproles_type)[agrep("Farmer",levels(farmingexp_df$sheeproles_type))]<-'Farmer/Stockperson/Vet'
levels(farmingexp_df$sheeproles_type)[agrep("Stockman/woman/person",levels(farmingexp_df$sheeproles_type))]<-'Farmer/Stockperson/Vet'
levels(farmingexp_df$sheeproles_type)[agrep("Veterinarian",levels(farmingexp_df$sheeproles_type))]<-'Farmer/Stockperson/Vet'
levels(farmingexp_df$sheeproles_type)[levels(farmingexp_df$sheeproles_type)%in%("Other;")]<-'Other'
levels(farmingexp_df$sheeproles_type)

#reorder levels of lameness prevalence
farmingexp_df$sheeplamenessexp_prev<-factor(farmingexp_df$sheeplamenessexp_prev, levels=c("Under 2%", "Between 2 and 5%", "Between 5 and 10%","Over 10%"))
#combine lameness prevs cos poor sampling at tails
levels(farmingexp_df$sheeplamenessexp_prev)[levels(farmingexp_df$sheeplamenessexp_prev)%in%c("Under 2%", "Between 2 and 5%")]<-"Under 5%"
levels(farmingexp_df$sheeplamenessexp_prev)[levels(farmingexp_df$sheeplamenessexp_prev)%in%c("Between 5 and 10%","Over 10%")]<-"Over 5%"

# PARSE SYMPTOMS -------------------------------------------------------------

symptoms<-c(
          'Uneven posture',
          'Shortened stride on one leg when walking',
          'Pair of legs which were moving at different speeds',
          'Nodding of head',
          'Not weight bearing on affected leg when standing',
          'Not weight bearing on affected leg when walking',
          'Reluctance to move',
          'Slower walking pace',
          'Other')

symptomslookedfor_df<-data.frame(matrix(nrow=nrow(studydata_formatted), ncol=length(symptoms)))
colnames(symptomslookedfor_df)<-symptoms
for (s in 1:length(symptoms)){
  print(s)
  temp<-grepl(symptoms[s],studydata_formatted$symptoms_lookedfor)  
  symptomslookedfor_df[,s]<-temp
}

symptomslookedfor_df

colnames(symptomslookedfor_df)[which(colnames(symptomslookedfor_df)%in%c("Pair of legs which were moving at different speeds","Not weight bearing on affected leg when standing","Not weight bearing on affected leg when walking"))]<-c("Pair of legs moving at different speeds","Not weight bearing on affected leg (standing)", "Not weight bearing on affected leg (walking)")
symptoms[which(symptoms%in%c("Pair of legs which were moving at different speeds","Not weight bearing on affected leg when standing","Not weight bearing on affected leg when walking"))]<-c("Pair of legs moving at different speeds","Not weight bearing on affected leg (standing)", "Not weight bearing on affected leg (walking)")

# PARSE FARMING EXPERIENCE/ROLES -------------------------------------------------------------

roles<-c(
  'Farmer',
  'Stockman/woman/person',
  'Veterinarian',
  'Other')

roles_df<-data.frame(matrix(nrow=nrow(studydata_formatted), ncol=length(roles)))
colnames(roles_df)<-roles
for (r in 1:length(roles)){
  print(r)
  temp<-grepl(roles[r],studydata_formatted$sheeproles_type)  
  roles_df[,r]<-temp
}

palette<-c("#56B4E9","#CC79A7") #non-farmer blue, farmer pink
farmingexp_cols<-palette[as.factor(studydata_formatted$farmingexp_YN)]

# ADDING VARIABLES TO THE DATAFRAME--------------------------------------------------------------

#arcsin-transformed recall
studydata_formatted$recall_asin<-asin(sqrt(0.01*studydata_formatted$recall))
#arcsin-transformed accuracy
#studydata_formatted$accuracy_asin<-asin(sqrt(0.01*studydata_formatted$accuracy))
#timespent
timespent<-(600-studydata_formatted$timeremaining_s)/60
studydata_formatted<-cbind(studydata_formatted,timespent)
#symptoms looked for
studydata_formatted<-cbind(studydata_formatted,symptomslookedfor_df)
colnames(studydata_formatted)[48:56]<-c('UP','SS','LS','NH','WBS','WBW','RM','SW','O')
#colours
studydata_formatted$farmingexp_cols<-farmingexp_cols

# STUDY PARTICIPANTS ------------------------------------------------------

plyr::count(as.factor(studydata_formatted$sheeproles_type))

layout_matrix_1 <- matrix(1:1, ncol = 2) 
layout(layout_matrix_1)

accuracyrecall<-c(studydata_formatted$accuracy,studydata_formatted$recall)
accuracyrecall_group<-c(rep('Accuracy',63),rep('Recall',63))

write.csv(studydata_formatted,'outputs/processed_data/studydata_formatted.csv')

# FIGURE: LIKERT GRAPH ----------------------------------------------------

likertdata

sapply(likertdata, levels)

#make levels of each column the same
for(i in 1:ncol(likertdata)) { 
  likertdata[,i] <- factor(likertdata[,i], levels=likert_categories_ordered)
}

sapply(likertdata, levels)
likert_plot<-likert::likert(likertdata)
library(plyr)

grDevices::tiff(paste('outputs/figures/likertplot.tiff', sep=''), res=300, units='in', width=10, height=8)
plot(likert_plot)
likert::likert.options()

dev.off()

# FIGURE: ACCURACY VS RECALL ----------------------------------------------

grDevices::tiff(paste('outputs/figures/accuracyvsrecall.tiff', sep=''), res=300, units='in', width=10, height=8)
par(mar=c(4,8,4,4))
beeswarm::beeswarm(accuracyrecall~accuracyrecall_group, 
                   pch = 19, cex=2.5, pwcol=rep(farmingexp_cols,2),
                   xlab='', ylab='',
                   xlim=c(0.5,2.5), 
                   las=2, xaxt='n',yaxt='n',
                   method='compactswarm')

axis(1, c(1,2), c('Accuracy', 'Recall'), line=1, tick=F, cex.axis=2)
axis(2, seq(0,100,10), seq(0,100,10), las=2, cex.axis=1.5)
title(ylab='Score (%)', cex.lab=2, line=4)

mtext('A', side=3, cex=2, line=0.5, adj=-0.08)
dev.off()

# FIGURE: FARMING EXPERIENCE ----------------------------------------------

grDevices::tiff(paste('outputs/figures/farmingexperience.tiff', sep=''), res=300, units='in', width=10, height=6)

axis.indices<-1
plot.indices<-2:3
panel.rows<-2
panel.cols<-2

layout.matrix<-matrix(plot.indices,
                      nrow=2,
                      ncol=1, byrow=F)
#rep rows
layout.matrix<-layout.matrix[, rep(1, each=2)]
#rep cols
layout.matrix<-layout.matrix[rep(1:2, each=2),]
#add first plot columns
layout.matrix<-cbind(matrix(1,nrow=4,ncol=2),layout.matrix)
#add axis column
#layout.matrix<-cbind(rep(axis.indices,nrow(layout.matrix)),layout.matrix)
#check
layout.matrix
layout(layout.matrix)

par(mar=c(7,10,2,0.5))

#plot.new()
#text(0.5,0.5,'Recall/percentage of lame sheep identified (arcsine-transformed)', cex=2, srt=90)

#FIGURE A
boxplot(studydata_formatted$recall~farmingexp_df$farmingexp_YN,
        pch = 19, cex=2, pwcol=studydata_formatted$farmingexp_cols,ylim=c(0,100),las=2,cex.axis=1.75,
        xlab = '', ylab='',
        xaxt='n', notch=T, cex.lab=2, tick=F)
text(0.65,97.5,'A', cex=5)

#add points
beeswarm::beeswarm(studydata_formatted$recall~farmingexp_df$farmingexp_YN,
                   pch = 19, cex=2, pwcol=studydata_formatted$farmingexp_cols,
                   xlab = 'Have you ever worked in farming or a related field (e.g. farm vet)?', ylab='',
                   notch=F, add=T)

axis(1, 1:2, c("No", "Yes"), cex.axis=1.75, padj=1)
mtext('Experience in farming/related field?', side=1, cex=1.25, line=5)
mtext('Recall/percentage of lame sheep identified', side=2, cex=1.25, line=5)


#FIGURE B
#reset margins
par(mar=c(5,5,0.5,2))
boxplot(studydata_formatted$recall~farmingexp_df$sheeplamenessexp_prev,
        pch = 19, cex=2, pwcol=studydata_formatted$farmingexp_cols, ylim=c(0,100),las=2, cex.axis=1.75,
        xlab = '', ylab='',
        xaxt='n',
        notch=T)
text(0.55,95,'B', cex=3)
beeswarm::beeswarm(studydata_formatted$recall~farmingexp_df$sheeplamenessexp_prev,
                   pch = 19, cex=2, pwcol=studydata_formatted$farmingexp_cols,
                   xaxt='n',
                   xlab = 'What do you think was the average level of lameness in the flock(s) with which you worked/work, over one year?', ylab='',
                   notch=F, add=T)

axis(1, 1:2, c("Under 5%", "Over 5%"), cex.axis=1.25)
mtext('Perceived annual prevalence of lameness experienced', side=1, cex=1, line=3)

#FIGURE C
plot(studydata_formatted$recall~farmingexp_df$sheepexp_years,
     pch = 19, cex=2, col=studydata_formatted$farmingexp_cols,ylim=c(0,100), las=2, cex.axis=1.75,
     xaxt='n',
     xlab = '', ylab='')
abline(lm(studydata_formatted$recall~farmingexp_df$sheepexp_years), lwd=5)

axis(1, seq(0,30,5),seq(0,30,5), cex.axis=1.25)
mtext('Years working with sheep', side=1, cex=1, line=3)
text(1,95,'C', cex=3)

dev.off()

# FIGURE: SYMPTOMS LOOKED FOR ----------------------------------------------

grDevices::tiff(paste('outputs/figures/symptoms.tiff', sep=''), res=300, units='in', width=10, height=12)

axis.indices<-1
plot.indices<-2:10

layout.matrix<-matrix(plot.indices,
                        nrow=3,
                        ncol=3, byrow=T)
#rep rows
layout.matrix<-layout.matrix[, rep(1:3, each=6)]
#rep cols
layout.matrix<-layout.matrix[rep(1:3, each=6),]
#add axis column
layout.matrix<-cbind(rep(axis.indices,nrow(layout.matrix)),layout.matrix)
#check
layout.matrix

layout(layout.matrix)

par(mar=c(5,0.5,0.5,1))

plot.new()
text(0.5,0.5,'Recall/percentage of lame sheep identified', cex=3, srt=90)

#recall by symptoms
for (s in 1:length(symptoms)){
  #recall by moving type
  boxplot(studydata_formatted$recall~as.factor(symptomslookedfor_df[,s]),
                                pch = 19, cex=2, pwcol=studydata_formatted$farmingexp_cols,
                                xlab = symptoms[s], ylab='',
                                xaxt='n', yaxt='n', ylim=c(0,100), cex.lab=1.25,
                                notch=T)
  
  beeswarm::beeswarm(studydata_formatted$recall~as.factor(symptomslookedfor_df[,s]),
                                pch = 19, cex=1, pwcol=studydata_formatted$farmingexp_cols,
                                add=T)
  axis(1, 1:2, c("No", "Yes"))
  axis(2,seq(0,100,20),NA,tick = T)
  means<-tapply(studydata_formatted$recall,as.factor(symptomslookedfor_df[,s]), mean)
  mids<-1:nlevels(as.factor(symptomslookedfor_df[,s]))
  segments(mids-0.25,means,mids+0.25,means, lwd=2, lty=3, col="black")
}

dev.off()

# FIGURE: USER ENGAGEMENT ---------------------------------------------------

grDevices::tiff(paste('outputs/figures/UE.tiff', sep=''), res=300, units='in', width=12, height=6)

HCI_df<-studydata_formatted %>% select(timesplayed, tutorial, cpu_setup, controlsprobs_YN, strategy_type, moving_type, timespent)

#change cpu setup levels
HCI_df$cpu_setup[agrep('Mouse',HCI_df$cpu_setup)]<-'Mouse'
HCI_df$cpu_setup[agrep('Track-pad',HCI_df$cpu_setup)]<-'Trackpad'
HCI_df$cpu_setup[agrep('Laptop or Desktop computer;',HCI_df$cpu_setup)]<-'Unknown'
HCI_df$cpu_setup<-factor(HCI_df$cpu_setup, levels=c('Unknown','Trackpad','Mouse'))

#change column classes
HCI_df<-HCI_df %>%
  mutate(across(c(tutorial,
                  cpu_setup,
                  controlsprobs_YN, 
                  strategy_type,
                  moving_type),
                factor))

levels(HCI_df$strategy_type)<-c('Up-close','Zoomer','Other')
levels(HCI_df$moving_type)<-c('Other','Randomly','Semi-randomly','Linear')
HCI_df$moving_type<-factor(HCI_df$moving_type, levels=c('Linear','Semi-randomly','Randomly','Other'))


axis.indices<-1
plot.indices<-1:6
panel.rows<-2
panel.cols<-2

layout.matrix<-matrix(plot.indices,
                      nrow=2,
                      ncol=3, byrow=T)
#rep rows
# layout.matrix<-layout.matrix[, rep(1:2, each=2)]
# #rep cols
# layout.matrix<-layout.matrix[rep(1:2, each=2),]
# #add first plot columns
# layout.matrix<-cbind(matrix(1,nrow=4,ncol=4),layout.matrix)
#add axis column
#layout.matrix<-cbind(rep(axis.indices,nrow(layout.matrix)),layout.matrix)
#check
layout.matrix
layout(layout.matrix)

par(mar=c(5,5,0.5,0.5))

#FIGURE A
boxplot(studydata_formatted$recall~HCI_df$controlsprobs_YN,
        pch = 19, cex=2, pwcol=studydata_formatted$farmingexp_cols, 
        xaxt='n',xlab='',
        ylim=c(0,100), ylab='',las=2,
        notch=T)
beeswarm::beeswarm(studydata_formatted$recall~HCI_df$controlsprobs_YN,
                   pch = 19, cex=1, pwcol=studydata_formatted$farmingexp_cols,
                   xlab = '', ylab='',
                   notch=F, add=T)
axis(1, 1:2, c("No", "Yes"), cex.axis=1.25, padj=0)
mtext("Problems with controls?", side=1, cex=1, line=3.5)
text(0.6,95,'A', cex=3)

#FIGURE B
boxplot(studydata_formatted$recall~HCI_df$strategy_type,
        pch = 19, cex=2, pwcol=studydata_formatted$farmingexp_cols,
        xaxt='n',xlab='',
        ylim=c(0,100), ylab='',las=2,
        notch=T)
beeswarm::beeswarm(studydata_formatted$recall~HCI_df$strategy_type,
                   pch = 19, cex=1, pwcol=studydata_formatted$farmingexp_cols,
                   xlab = '', ylab='',
                   notch=F, add=T)
axis(1, 1:nlevels(HCI_df$strategy_type), levels(HCI_df$strategy_type), cex.axis=1.25, padj=0)
mtext("Observing type", side=1, cex=1, line=3.5)
text(0.6,95,'B', cex=3)

#FIGURE C
boxplot(studydata_formatted$recall~HCI_df$moving_type,
        pch = 19, cex=2, pwcol=studydata_formatted$farmingexp_cols,
        xaxt='n',xlab='',
        ylim=c(0,100), ylab='',las=2,
        notch=T)
beeswarm::beeswarm(studydata_formatted$recall~HCI_df$moving_type,
                   pch = 19, cex=1, pwcol=studydata_formatted$farmingexp_cols,
                   xlab = 'What do you think was the average level of lameness in the flock(s) with which you worked/work, over one year?', ylab='',
                   notch=F, add=T)
axis(1, 1:nlevels(HCI_df$moving_type), levels(HCI_df$moving_type), cex.axis=1, padj=0)
mtext("Moving type", side=1, cex=1, line=3.5)
text(0.6,95,'C', cex=3)

#FIGURE D
boxplot(studydata_formatted$recall~HCI_df$tutorial,
        pch = 19, cex=2, pwcol=studydata_formatted$farmingexp_cols,
        xaxt='n',xlab='',
        ylim=c(0,100), ylab='',las=2,
        notch=T)
beeswarm::beeswarm(studydata_formatted$recall~HCI_df$tutorial,
                   pch = 19, cex=1, pwcol=studydata_formatted$farmingexp_cols,
                   xlab = '', ylab='',
                   notch=F, add=T)
axis(1, 1:nlevels(as.factor(HCI_df$tutorial)), c('No', 'Yes \n(observed walking)','Yes \n(did not observe)'), cex.axis=1, padj=0.5)
mtext("Completed pre-game tutorial?", side=1, cex=1, line=3.5)
text(0.6,95,'D', cex=3)

#FIGURE E
boxplot(studydata_formatted$recall~HCI_df$cpu_setup,
        pch = 19, cex=2, pwcol=studydata_formatted$farmingexp_cols,
        xaxt='n',xlab='',
        ylim=c(0,100), ylab='',las=2,
        notch=T)
beeswarm::beeswarm(studydata_formatted$recall~HCI_df$cpu_setup,
                   pch = 19, cex=1, pwcol=studydata_formatted$farmingexp_cols,
                   xlab = '', ylab='',
                   notch=F, add=T)
axis(1, 1:nlevels(HCI_df$cpu_setup), levels(HCI_df$cpu_setup), cex.axis=1.25, padj=0)
mtext("Computer set-up", side=1, cex=1, line=3.5)
text(0.6,95,'E', cex=3)

#FIGURE F
plot(studydata_formatted$recall~HCI_df$timespent,
     col=studydata_formatted$farmingexp_cols, pch=19, xlim=c(0,10), las=2,
     xaxt='n',xlab='',
     ylim=c(0,100), ylab='',las=2,
     cex.lab=1.25)

axis(1, seq(0,10,1),seq(0,10,1), cex.axis=1.25)
mtext('Time spent playing (minutes)', side=1, cex=1, line=3.5)

m1<-lm(studydata_formatted$recall~HCI_df$timespent)
abline(m1)

text(0.6,95,'F', cex=3)

dev.off()

# POWER ANALYSIS ----------------------------------------------------------

#Cohen suggests that f values of 0.1, 0.25, and 0.4 represent small, medium, and large
#power analysis for a simple linear regression with stringent alpha and beta values
pwr_a1<-pwr::pwr.f2.test(u=1, v=61, sig.level = 0.05, power=0.95) 
pwr_a1
saveRDS(pwr_a1, 'outputs/models/power_analysis.RDS')
effectsize::f2_to_eta2(pwr_a1$f2)

# LINEAR MODELS  -------------------------------------

layout.matrix<-matrix(1:4,nrow=2,ncol=2)
layout(layout.matrix)

farmingexperience_model<-lm(studydata_formatted$recall_asin~farmingexp_df$farmingexp_YN)
saveRDS(farmingexperience_model, 'outputs/models/farmingexperience_model.RDS')
summary(farmingexperience_model)
plot(farmingexperience_model)

symptom_model.1<-lm(studydata_formatted$recall_asin~symptomslookedfor_df$`Pair of legs moving at different speeds`)
saveRDS(symptom_model.1, 'outputs/models/symptom_limp_model.RDS')
summary(symptom_model.1)
plot(symptom_model.1)

symptom_model.2<-lm(studydata_formatted$recall_asin~symptomslookedfor_df$`Not weight bearing on affected leg (standing)`)
saveRDS(symptom_model.2, 'outputs/models/symptom_raisedleg_model.RDS')
summary(symptom_model.2)
plot(symptom_model.2)

userengagement_model<-lm(studydata_formatted$recall_asin~HCI_df$timespent)
saveRDS(userengagement_model, 'outputs/models/userengagement_model.RDS')
summary(userengagement_model)
plot(userengagement_model)

#model comparison
model_comparison<-anova(farmingexperience_model,symptom_model.1,userengagement_model)
saveRDS(model_comparison, 'outputs/models/model_comparison.RDS')

