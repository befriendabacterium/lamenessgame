library(scales)
library(dplyr)
library(ggplot2)

# LOAD CUSTOM FUNCTIONS ----------------------------------------------------------

actualrange<-function(x){range(x)[2]-range(x)[1]}
#calculate lower hinge (smallest data value that is larger than the first quartile: https://webhelp.esri.com/arcgisdesktop/9.3/body.cfm?tocVisable=1&ID=478&TopicName=Box%20plot%20graphs)
lower_hinge<-function(x){min(x[x>quantile(x)[2]])}
#calculate lower hinge (largest data value that is smaller than the first quartile: https://webhelp.esri.com/arcgisdesktop/9.3/body.cfm?tocVisable=1&ID=478&TopicName=Box%20plot%20graphs)
upper_hinge<-function(x){min(x[x<quantile(x)[4]])}

# READ IN RAW DATA FROM MS FORMS --------------------------------------------------------

studydata_raw<-readxl::read_excel('inputs/lamenessstudydata_010721.xlsx')
#detach(tibble)
newcolnames<-readRDS('inputs/newcolnames.RDS')

# FORMATTING: PARTICIPANT PERFORMANCE DATA --------------------------------------------------------------

#REDUNDANT COLUMN REMOVAL (EASIER TO WORK WITH)
#make a vector of columns to remove
colstoremove<-grep(c('Email|Name|Total.points|Quiz.feedback|Points...|Feedback...|Finally...'), colnames(studydata_raw))
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

# ARCSINE-TRANSFORM SCORES --------------------------------------------------------------

#arcsin-transformed accuracy
studydata_formatted$accuracy_asin<-asin(sqrt(0.01*studydata_formatted$accuracy))
#arcsin-transformed recall
studydata_formatted$recall_asin<-asin(sqrt(0.01*studydata_formatted$recall))

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

#reorder levels of lameness prevalence
studydata_formatted$sheeplamenessexp_prev<-factor(studydata_formatted$sheeplamenessexp_prev, levels=c("Under 2%", "Between 2 and 5%", "Between 5 and 10%","Over 10%"))
#combine lameness prevs cos poor sampling at tails
#levels(studydata_formatted$sheeplamenessexp_prev)[levels(studydata_formatted$sheeplamenessexp_prev)%in%c("Under 2%", "Between 2 and 5%")]<-"Under 5%"
#levels(studydata_formatted$sheeplamenessexp_prev)[levels(studydata_formatted$sheeplamenessexp_prev)%in%c("Between 5 and 10%","Over 10%")]<-"Over 5%"

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

# PARSE LAMENESS SIGNS (SYMPTOMS) -------------------------------------------------------------

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

#add to main dataframe
studydata_formatted<-cbind(studydata_formatted,symptomslookedfor_df)
colnames(studydata_formatted)[47:55]<-c('UP','SS','LS','NH','WBS','WBW','RM','SW','O')

#remove the original symptoms_lookedfor column to remove redundancy in the dataframe
studydata_formatted<-studydata_formatted[,-which(colnames(studydata_formatted)=='symptoms_lookedfor')]

# PARSE USER ENGAGEMENT ---------------------------------------------------

#calculate and add timespent
studydata_formatted$timespent<-(600-studydata_formatted$timeremaining_s)/60
#remove the original timeremaining_s column to remove redundancy in the dataframe
studydata_formatted<-studydata_formatted[,-which(colnames(studydata_formatted)=='timeremaining_s')]

#change cpu setup levels
studydata_formatted$cpu_setup[agrep('Mouse',studydata_formatted$cpu_setup)]<-'Mouse'
studydata_formatted$cpu_setup[agrep('Track-pad',studydata_formatted$cpu_setup)]<-'Trackpad'
studydata_formatted$cpu_setup[agrep('Laptop or Desktop computer;',studydata_formatted$cpu_setup)]<-'Unknown'
studydata_formatted$cpu_setup<-factor(studydata_formatted$cpu_setup, levels=c('Unknown','Trackpad','Mouse'))

#change columns to factor
studydata_formatted$tutorial<-as.factor(studydata_formatted$tutorial)
studydata_formatted$cpu_setup<-studydata_formatted$cpu_setup
studydata_formatted$controlsprobs_YN<-as.factor(studydata_formatted$controlsprobs_YN)
studydata_formatted$strategy_type<-as.factor(studydata_formatted$strategy_type)
studydata_formatted$moving_type<-as.factor(studydata_formatted$moving_type)

#edit levels
levels(studydata_formatted$strategy_type)<-c('Up-close','Zoomer','Other')
levels(studydata_formatted$moving_type)<-c('Other','Randomly','Semi-randomly','Linear')
studydata_formatted$moving_type<-factor(studydata_formatted$moving_type, levels=c('Linear','Semi-randomly','Randomly','Other'))

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

#write csv for the formatted likert data
write.csv(likertdata_formatted, 'outputs/processed_data/likertdata_formatted.csv')

#ugly way of identifying likert columns in study data
likert_cols<-c()
for (l in 1:ncol(studydata_formatted)){
  
  current<-agrep(gsub(' ','.',colnames(likertdata)[l]),colnames(studydata_formatted))
  #provided there's a match, bind it
  if(!any(is.na(current))){
    likert_cols<-c(likert_cols,current)
  }

  
}

# TIDYING UP --------------------------------------------------------------

#remove likert columns from the main studydata dataframe to keep things tidy (we already have them, in a better format, in the other dataframe)
studydata_formatted<-studydata_formatted[,-likert_cols]

#let's reorder the columns a little to make the formatted dataframe more human-readable
colnames(studydata_formatted)

#move farming experience columns to end
studydata_formatted<-
  relocate(studydata_formatted, 
           c('accuracy','accuracy_asin','recall','recall_asin'),
           .after = colnames(studydata_formatted)[ncol(studydata_formatted)])

#move farming experience columns to end
studydata_formatted<-
  relocate(studydata_formatted, 
         c('farmingexp_YN', 'sheepexp_years', 'sheeproles_type', 'sheeproles_details','sheeplamenessexp_prev'),
         .after = colnames(studydata_formatted)[ncol(studydata_formatted)])

#move signs looked for columns to end
studydata_formatted<-
  relocate(studydata_formatted, 
           c('UP','SS','LS','NH','WBS','WBW','RM','SW','O', 'symptoms_otherdetails'),
           .after = colnames(studydata_formatted)[ncol(studydata_formatted)])

#move user engagement columns to end
studydata_formatted<-
  relocate(studydata_formatted, 
           c('timesplayed','tutorial','cpu_setup','controlsprobs_YN','controlsprobs_descrip','strategy_type','strategy_otherdetails','moving_type','moving_otherdetails', 'timespent'),
           .after = colnames(studydata_formatted)[ncol(studydata_formatted)])

#move feedback column to end
studydata_formatted<-
  relocate(studydata_formatted, 
           c('Please.share.any.general.feedback.main.thoughts.after.playing.the.game.below.'),
           .after = colnames(studydata_formatted)[ncol(studydata_formatted)])

colnames(studydata_formatted)

#write the csv
write.csv(studydata_formatted,'outputs/processed_data/studydata_formatted.csv')

# POWER ANALYSIS ----------------------------------------------------------

#Cohen suggests that f values of 0.1, 0.25, and 0.4 represent small, medium, and large
#power analysis for a simple linear regression with stringent alpha and beta values
pwr_a1<-pwr::pwr.f2.test(u=1, v=61, sig.level = 0.05, power=0.95) 
pwr_a1
saveRDS(pwr_a1, 'outputs/models/power_analysis.RDS')
effectsize::f2_to_eta2(pwr_a1$f2)

# SKEWNESS TEST: ACCURACY VS RECALL ---------------------------------------

#test for skewness in accuracy and recall scores
moments::agostino.test(studydata_formatted$accuracy)
moments::agostino.test(studydata_formatted$recall)

# FIGURE: ACCURACY VS RECALL ----------------------------------------------

accuracyvsrecall<-c(studydata_formatted$accuracy,studydata_formatted$recall)
accuracyvsrecall<-cbind(c(rep('accuracy',length(accuracyvsrecall)/2), rep('recall',length(accuracyvsrecall)/2)),
                        accuracyvsrecall)
accuracyvsrecall<-as.data.frame(accuracyvsrecall)
colnames(accuracyvsrecall)<-c('score_type','score')
accuracyvsrecall$score<-as.numeric(accuracyvsrecall$score)

grDevices::tiff(paste('outputs/figures/accuracyvsrecall.tiff', sep=''), res=300, units='in', width=8, height=8)
par(mar=c(4,6,4,4))

#INITIATE PLOT
plot(1~1,
     xlim=c(0,1),ylim=c(0,100),
    xlab='', ylab='',
    xaxt='n',yaxt='n',
    type='n')

#swarm locations
swarm_locs<-c(0.25,0.75)

#ADD MEAN LINES
#calculate means per group
means<-tapply(accuracyvsrecall$score,as.factor(accuracyvsrecall$score_type), mean)
#calculate means per group
cis<-tapply(accuracyvsrecall$score,as.factor(accuracyvsrecall$score_type), Rmisc::CI)
#calculate lower and upper hinges using custom functions
lowerhinges<-tapply(accuracyvsrecall$score,as.factor(accuracyvsrecall$score_type), lower_hinge)
upperhinges<-tapply(accuracyvsrecall$score,as.factor(accuracyvsrecall$score_type), upper_hinge)
notches_matrix<-cbind(unlist(lapply(cis,min))<lowerhinges,unlist(lapply(cis,max))<lowerhinges)
#would it have notches in a box plot? i.e. is confidence interval (range) is more than the inter-quartile range
notches<-rowSums(notches_matrix)
#initiate mean colour column by adding 1 to true/false notches (so false is lty=1)
mean_cols<-notches+1
#add means with appropriate colour to plot
segments(swarm_locs-0.125,means,swarm_locs+0.125,means, lwd=3, col=mean_cols)

#ADD BEESWARM POINTS
beeswarm::beeswarm(accuracyvsrecall$score~accuracyvsrecall$score_type, 
                   at=swarm_locs, pwcol=rep('darkgrey',length(accuracyvsrecall$score)),
                   pch = 19, cex=1.5,
                   xlim=swarm_locs, add=T,
                   las=2, xaxt='n',yaxt='n',
                   method='compactswarm')

axis(1, swarm_locs, c('Accuracy', 'Recall'), line=0, tick=T, cex.axis=1.5)
axis(2, seq(0,100,10), seq(0,100,10), las=2, cex.axis=1.5)
title(ylab='Score (%)', cex.lab=1.5, line=4)

dev.off()

# FIGURE: FARMING EXPERIENCE ----------------------------------------------

grDevices::tiff(paste('outputs/figures/farmingexperience.tiff', sep=''), res=300, units='in', width=9, height=6)

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

#FIGURE A
#initiate plot
plot(1~1,
     xlim=c(0,1),ylim=c(0,100),
     xlab='', ylab='',
     xaxt='n',yaxt='n',
     type='n')
#because this is one of the tested relationships, frame the plot with bold border
  box(lwd=3)

text(0.05,97.5,'A', cex=5)
#swarm locations
swarm_locs<-c(0.25,0.75)


#ADD MEAN LINES
#calculate means per group
means<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$farmingexp_YN), mean)
#calculate means per group
cis<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$farmingexp_YN), Rmisc::CI)
#calculate lower and upper hinges using custom functions
lowerhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$farmingexp_YN), lower_hinge)
upperhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$farmingexp_YN), upper_hinge)
notches_matrix<-cbind(unlist(lapply(cis,min))<lowerhinges,unlist(lapply(cis,max))<lowerhinges)
#would it have notches in a box plot? i.e. is confidence interval (range) is more than the inter-quartile range
notches<-rowSums(notches_matrix)
#initiate mean colour column by adding 1 to true/false notches (so false is lty=1)
mean_cols<-notches+1
#add means with appropriate colour to plot
segments(swarm_locs-0.125,means,swarm_locs+0.125,means, lwd=3, col=mean_cols)


#ADD BEESWARM POINTS
beeswarm::beeswarm(studydata_formatted$recall~studydata_formatted$farmingexp_YN,
                   at=swarm_locs, pwcol=rep('darkgrey',length(studydata_formatted$farmingexp_YN)),
                   pch = 19, cex=2,
                   xlab = 'Have you ever worked in farming or a related field (e.g. farm vet)?', ylab='',
                   notch=F, add=T)

axis(1, swarm_locs, c("No", "Yes"), cex.axis=1.5, padj=1)
axis(2, seq(0,100,10), seq(0,100,10), las=2, cex.axis=1.5)
mtext('Experience in farming/related field?', side=1, cex=1.25, line=5)
mtext('Recall/percentage of lame sheep identified', side=2, cex=1.25, line=5)

#FIGURE B
#reset margins
par(mar=c(5,5,0.5,2))
plot(studydata_formatted$recall~studydata_formatted$sheepexp_years,
     col=rep('darkgrey',length(studydata_formatted$sheepexp_years)),
     pch = 19, cex=2,ylim=c(0,100), las=2, cex.axis=1.5,
     xaxt='n',
     xlab = '', ylab='')

text(0.5,95,'B', cex=3)
axis(1, seq(0,30,5),seq(0,30,5), cex.axis=1.5)
mtext('Years working with sheep', side=1, cex=1, line=3)


#FIGURE C

#reset margins
par(mar=c(5,5,0.5,2))
#initiate plot
plot(1~1,
     xlim=c(0,2),ylim=c(0,100),
     xlab='', ylab='',
     xaxt='n',yaxt='n',
     type='n')
text(0.025,95,'C', cex=3)
#swarm locations
swarm_locs<-c(0.25,0.75, 1.25, 1.75)

#ADD MEAN LINES
#calculate means per group
means<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$sheeplamenessexp_prev), mean)
#calculate means per group
cis<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$sheeplamenessexp_prev), Rmisc::CI)
#calculate lower and upper hinges using custom functions
lowerhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$sheeplamenessexp_prev), lower_hinge)
upperhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$sheeplamenessexp_prev), upper_hinge)
notches_matrix<-cbind(unlist(lapply(cis,min))<lowerhinges,unlist(lapply(cis,max))<lowerhinges)
#would it have notches in a box plot? i.e. is confidence interval (range) is more than the inter-quartile range
notches<-rowSums(notches_matrix)
#initiate mean colour column by adding 1 to true/false notches (so false is lty=1)
mean_cols<-notches+1
#add means with appropriate colour to plot
segments(swarm_locs-0.125,means,swarm_locs+0.125,means, lwd=3, col=mean_cols)

#ADD BEESWARM POINTS
beeswarm::beeswarm(studydata_formatted$recall~studydata_formatted$sheeplamenessexp_prev,
                   at=swarm_locs, pwcol=rep('darkgrey',length(studydata_formatted$farmingexp_YN)),
                   pch = 19, cex=2,
                   xaxt='n',
                   xlab = 'What do you think was the average level of lameness in the flock(s) with which you worked/work, over one year?', ylab='',
                   notch=F, add=T)

axis(1, swarm_locs, c("Under 2%", "Between 2 and 5%", "Between 5 and 10%","Over 10%"), cex.axis=0.8)
axis(2, seq(0,100,20), seq(0,100,20), las=2, cex.axis=1.5)
mtext('Perceived annual prevalence of lameness experienced', side=1, cex=1, line=3)

dev.off()

# MODEL(S): FARMING EXPERIENCE --------------------------------------------

#make a test counter to keep track of number of tests and adjust p values accordingly
tests<-0

#make model on arcsine-transformed data
farmingexperience_model<-lm(studydata_formatted$recall_asin~studydata_formatted$farmingexp_YN)
#add 1 to the test counter
tests<-tests+1

#summarise and adjust p values
summary(farmingexperience_model)
p.adjust(summary(farmingexperience_model)$coefficients[,"Pr(>|t|)"][2], method = 'bonferroni', n = tests)

#save model
saveRDS(farmingexperience_model, 'outputs/models/farmingexperience_model.RDS')

#plot model
layout.matrix<-matrix(1:4,nrow=2,ncol=2)
layout(layout.matrix)
plot(farmingexperience_model)

#summarise and adjust p values
summary(farmingexperience_model)
p.adjust(summary(farmingexperience_model)$coefficients[,"Pr(>|t|)"][2], method = 'bonferroni', n = 1)

# FIGURE: SYMPTOMS LOOKED FOR ----------------------------------------------

grDevices::tiff(paste('outputs/figures/symptoms.tiff', sep=''), res=300, units='in', width=9, height=12)

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
text(0.5,0.5,'Recall/percentage of lame sheep identified', cex=2.5, srt=90)

#recall by symptoms
for (s in 1:length(symptoms)){
  
  #initiate plot
  plot(1~1,
       xlim=c(0,1),ylim=c(0,100),
       xlab='', ylab='',
       xaxt='n',yaxt='n',
       type='n')
  #if the symptom looked at is one of the tested relationships, frame the plot with bold border
  if (s%in%c(1,3,5)){
    box(lwd=3)
  }
  #swarm locations
  swarm_locs<-c(0.25,0.75)
  
  #ADD MEAN LINES
  #calculate means per group
  means<-tapply(studydata_formatted$recall,as.factor(symptomslookedfor_df[,s]), mean)
  #calculate means per group
  cis<-tapply(studydata_formatted$recall,as.factor(symptomslookedfor_df[,s]), Rmisc::CI)
  #calculate lower and upper hinges using custom functions
  lowerhinges<-tapply(studydata_formatted$recall,as.factor(symptomslookedfor_df[,s]), lower_hinge)
  upperhinges<-tapply(studydata_formatted$recall,as.factor(symptomslookedfor_df[,s]), upper_hinge)
  notches_matrix<-cbind(unlist(lapply(cis,min))<lowerhinges,unlist(lapply(cis,max))<lowerhinges)
  #would it have notches in a box plot? i.e. is confidence interval (range) is more than the inter-quartile range
  notches<-rowSums(notches_matrix)
  #initiate mean colour column by adding 1 to true/false notches (so false is lty=1)
  mean_cols<-notches+1
  #add means with appropriate colour to plot
  segments(swarm_locs-0.125,means,swarm_locs+0.125,means, lwd=3, col=mean_cols)
  
  #ADD BEESWARM POINTS
  beeswarm::beeswarm(studydata_formatted$recall~symptomslookedfor_df[,s],
                     at=swarm_locs, pwcol=rep('darkgrey',length(studydata_formatted$farmingexp_YN)),
                     pch = 19, cex=2,
                     xaxt='n',
                     xlab = '', ylab='',
                     notch=F, add=T)

  axis(1, swarm_locs, c("No", "Yes"))
  axis(2,seq(0,100,20),NA,tick = T)
  mtext(symptoms[s], side=1, cex=0.85, line=3)
  
}

dev.off()

# MODEL(S): SYMPTOMS LOOKED FOR --------------------------------------------

#make model (uneven posture) on arcsine-transformed data
symptom_model.1<-lm(recall_asin~UP, data=studydata_formatted)
#add 1 to the test counter
tests<-tests+1
#summarise and adjust p values
summary(symptom_model.1)
symptom_model.1_adjp<-p.adjust(summary(symptom_model.1)$coefficients[,"Pr(>|t|)"][2], method = 'bonferroni', n = tests)
#save model
saveRDS(symptom_model.1, 'outputs/models/symptom_unevenposture_model.RDS')
saveRDS(symptom_model.1_adjp, 'outputs/models/symptom_unevenposture_adjp.RDS')
#plot model
layout.matrix<-matrix(1:4,nrow=2,ncol=2)
layout(layout.matrix)
plot(symptom_model.1)


#make model (limp) on arcsine-transformed data
symptom_model.2<-lm(recall_asin~LS, data=studydata_formatted)
#add 1 to the test counter
tests<-tests+1
#summarise and adjust p values
summary(symptom_model.2)
symptom_model.2_adjp<-p.adjust(summary(symptom_model.2)$coefficients[,"Pr(>|t|)"][2], method = 'bonferroni', n = tests)
#save model and adjusted p value
saveRDS(symptom_model.2, 'outputs/models/symptom_limp_model.RDS')
saveRDS(symptom_model.2_adjp, 'outputs/models/symptom_limp_adjp.RDS')
#plot model
layout.matrix<-matrix(1:4,nrow=2,ncol=2)
layout(layout.matrix)
plot(symptom_model.1)


#make model (raised leg)
symptom_model.3<-lm(recall_asin~WBS, data=studydata_formatted)
#add 1 to the test counter
tests<-tests+1
#summarise and adjust p values
summary(symptom_model.3)
symptom_model.3_adjp<-p.adjust(summary(symptom_model.3)$coefficients[,"Pr(>|t|)"][2], method = 'bonferroni', n = tests)
#save model
saveRDS(symptom_model.3, 'outputs/models/symptom_raisedleg_model.RDS')
saveRDS(symptom_model.3_adjp, 'outputs/models/symptom_raisedleg_adjp.RDS')
#plot model
layout.matrix<-matrix(1:4,nrow=2,ncol=2)
layout(layout.matrix)
plot(symptom_model.2)

# FIGURE: USER ENGAGEMENT ---------------------------------------------------

grDevices::tiff(paste('outputs/figures/UE.tiff', sep=''), res=300, units='in', width=11.5, height=5)

axis.indices<-1
plot.indices<-2:7
panel.rows<-2
panel.cols<-2

layout.matrix<-matrix(plot.indices,
                      nrow=2,
                      ncol=3, byrow=T)
#rep cols
layout.matrix<-layout.matrix[, rep(1:ncol(layout.matrix), each=6)]
# #rep rows
# layout.matrix<-layout.matrix[rep(1:2, each=2),]
# #add first plot columns
# layout.matrix<-cbind(matrix(1,nrow=4,ncol=4),layout.matrix)
#add axis column
layout.matrix<-cbind(rep(axis.indices,nrow(layout.matrix)),layout.matrix)
#add big panel for final plot
layout.matrix<-cbind(layout.matrix,matrix(8,nrow=2,ncol=13))

#check
layout.matrix
layout(layout.matrix)

#add margin text
par(mar=c(0,0,0,0))
plot.new()
text(0.5,0.5,'Recall/percentage of lame sheep identified', cex=1.5, srt=90)

par(mar=c(4.5,2,0.5,0))

#FIGURE A: INITIATE PLOT
plot(1~1,
     xlim=c(0,2.5),ylim=c(0,100),
     xlab='', ylab='',
     xaxt='n',yaxt='n',
     type='n')
axis(2,seq(0,100,20),seq(0,100,20), las=2)
text(0.05,95,'A', cex=2)
#swarm locations
swarm_locs<-c(0.25,0.75,1.25,1.75,2.25)

#ADD MEAN LINES
#calculate means per group
means<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$timesplayed), mean)
#calculate means per group
cis<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$timesplayed), Rmisc::CI)
#calculate lower and upper hinges using custom functions
lowerhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$timesplayed), lower_hinge)
upperhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$timesplayed), upper_hinge)
notches_matrix<-cbind(unlist(lapply(cis,min))<lowerhinges,unlist(lapply(cis,max))<lowerhinges)
#would it have notches in a box plot? i.e. is confidence interval (range) is more than the inter-quartile range
notches<-rowSums(notches_matrix)
#initiate mean colour column by adding 1 to true/false notches (so false is lty=1)
mean_cols<-notches+1
#if NA (no variance, make red too)
mean_cols[is.na(mean_cols)]<-2

#add means with appropriate colour to plot
segments(swarm_locs-0.125,means,swarm_locs+0.125,means, lwd=3, col=mean_cols)

#ADD BEESWARM POINTS
beeswarm::beeswarm(studydata_formatted$recall~studydata_formatted$timesplayed,
                   at=swarm_locs, pwcol=rep('darkgrey',length(studydata_formatted$timesplayed)),
                   pch = 19, cex=1,
                   xlab = '', ylab='',
                   notch=F, add=T)

axis(1, swarm_locs, c(1,2,3,4,5)) #n.b. +1 to labels to make more intuitive
mtext("Times played before submitting score", side=1, cex=0.7, line=2.5)

par(mar=c(4.5,1,0.5,1))

#FIGURE B: INITIATE PLOT
plot(1~1,
     xlim=c(0,1),ylim=c(0,100),
     xlab='', ylab='',
     xaxt='n',yaxt='n',
     type='n')
text(0.05,95,'B', cex=2)
#swarm locations
swarm_locs<-c(0.25,0.75)

#ADD MEAN LINES
#calculate means per group
means<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$controlsprobs_YN), mean)
#calculate means per group
cis<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$controlsprobs_YN), Rmisc::CI)
#calculate lower and upper hinges using custom functions
lowerhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$controlsprobs_YN), lower_hinge)
upperhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$controlsprobs_YN), upper_hinge)
notches_matrix<-cbind(unlist(lapply(cis,min))<lowerhinges,unlist(lapply(cis,max))<lowerhinges)
#would it have notches in a box plot? i.e. is confidence interval (range) is more than the inter-quartile range
notches<-rowSums(notches_matrix)
#initiate mean colour column by adding 1 to true/false notches (so false is lty=1)
mean_cols<-notches+1
#add means with appropriate colour to plot
segments(swarm_locs-0.125,means,swarm_locs+0.125,means, lwd=3, col=mean_cols)

#ADD BEESWARM POINTS
beeswarm::beeswarm(studydata_formatted$recall~studydata_formatted$controlsprobs_YN,
                   at=swarm_locs, pwcol=rep('darkgrey',length(studydata_formatted$controlsprobs_YN)),
                   pch = 19, cex=1,
                   xlab = '', ylab='',
                   notch=F, add=T)

axis(1, swarm_locs, c("No", "Yes"))
mtext("Problems with controls?", side=1, cex=0.7, line=2.5)

par(mar=c(4.5,1,0.5,1))

#FIGURE C: INITIATE PLOT
plot(1~1,
     xlim=c(0,1.5),ylim=c(0,100),
     xlab='', ylab='',
     xaxt='n',yaxt='n',
     type='n')
text(0.05,95,'C', cex=2)
#swarm locations
swarm_locs<-c(0.25,0.75,1.25)

#ADD MEAN LINES
#calculate means per group
means<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$strategy_type), mean)
#calculate means per group
cis<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$strategy_type), Rmisc::CI)
#calculate lower and upper hinges using custom functions
lowerhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$strategy_type), lower_hinge)
upperhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$strategy_type), upper_hinge)
notches_matrix<-cbind(unlist(lapply(cis,min))<lowerhinges,unlist(lapply(cis,max))<lowerhinges)
#would it have notches in a box plot? i.e. is confidence interval (range) is more than the inter-quartile range
notches<-rowSums(notches_matrix)
#initiate mean colour column by adding 1 to true/false notches (so false is lty=1)
mean_cols<-notches+1
#add means with appropriate colour to plot
segments(swarm_locs-0.125,means,swarm_locs+0.125,means, lwd=3, col=mean_cols)

#ADD BEESWARM POINTS TO PLOT
beeswarm::beeswarm(studydata_formatted$recall~studydata_formatted$strategy_type,
                   at=swarm_locs,pwcol=rep('darkgrey',length(studydata_formatted$strategy_type)),
                   pch = 19, cex=1,
                   xlab = '', ylab='',
                   notch=F, add=T)
axis(1, swarm_locs, levels(studydata_formatted$strategy_type), cex.axis=1, padj=0)
axis(2,seq(0,100,20),NA,tick = T)

mtext("Observing type", side=1, cex=0.7, line=2.5)

par(mar=c(4.5,2,0.5,0))

#FIGURE D: INITIATE PLOT
plot(1~1,
     xlim=c(0,2),ylim=c(0,100),
     xlab='', ylab='',
     xaxt='n',yaxt='n',
     type='n')
axis(2,seq(0,100,20),seq(0,100,20), las=2)

text(0.05,95,'D', cex=2)
#swarm locations
swarm_locs<-c(0.25,0.75,1.25,1.75)

#ADD MEAN LINES
#calculate means per group
means<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$moving_type), mean)
#calculate means per group
cis<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$moving_type), Rmisc::CI)
#calculate lower and upper hinges using custom functions
lowerhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$moving_type), lower_hinge)
upperhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$moving_type), upper_hinge)
notches_matrix<-cbind(unlist(lapply(cis,min))<lowerhinges,unlist(lapply(cis,max))<lowerhinges)
#would it have notches in a box plot? i.e. is confidence interval (range) is more than the inter-quartile range
notches<-rowSums(notches_matrix)
#initiate mean colour column by adding 1 to true/false notches (so false is lty=1)
mean_cols<-notches+1
#add means with appropriate colour to plot
segments(swarm_locs-0.125,means,swarm_locs+0.125,means, lwd=3, col=mean_cols)

beeswarm::beeswarm(studydata_formatted$recall~studydata_formatted$moving_type,
                   at=swarm_locs, pwcol=rep('darkgrey',length(studydata_formatted$moving_type)),
                   pch = 19, cex=1,
                   notch=F, add=T)
axis(1, swarm_locs, levels(studydata_formatted$moving_type), cex.axis=0.65, padj=0.2)
axis(2,seq(0,100,20),NA,tick = T)
mtext("Moving type", side=1, cex=0.7, line=2.5)

par(mar=c(4.5,1,0.5,1))

#FIGURE E: INITIATE PLOT
plot(1~1,
     xlim=c(0,1.5),ylim=c(0,100),
     xlab='', ylab='',
     xaxt='n',yaxt='n',
     type='n')
text(0.05,95,'E', cex=2)
#swarm locations
swarm_locs<-c(0.25,0.75,1.25)

#ADD MEAN LINES
#calculate means per group
means<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$tutorial), mean)
#calculate means per group
cis<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$tutorial), Rmisc::CI)
#calculate lower and upper hinges using custom functions
lowerhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$tutorial), lower_hinge)
upperhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$tutorial), upper_hinge)
notches_matrix<-cbind(unlist(lapply(cis,min))<lowerhinges,unlist(lapply(cis,max))<lowerhinges)
#would it have notches in a box plot? i.e. is confidence interval (range) is more than the inter-quartile range
notches<-rowSums(notches_matrix)
#initiate mean colour column by adding 1 to true/false notches (so false is lty=1)
mean_cols<-notches+1
#add means with appropriate colour to plot
segments(swarm_locs-0.125,means,swarm_locs+0.125,means, lwd=3, col=mean_cols)

#ADD BEESWARM POINTS
beeswarm::beeswarm(studydata_formatted$recall~studydata_formatted$tutorial,
                   at=swarm_locs, pwcol=rep('darkgrey',length(studydata_formatted$tutorial)),
                   pch = 19, cex=1,
                   xlab = '', ylab='',
                   notch=F, add=T)
axis(1, swarm_locs, c('No', 'Yes \n(observed walking)','Yes \n(did not observe)'), cex.axis=0.6, padj=0.2)
mtext("Completed pre-game tutorial?", side=1, cex=0.7, line=2.5)

par(mar=c(4.5,1,0.5,1))

#FIGURE F: INITIATE PLOT
plot(1~1,
     xlim=c(0,1.5),ylim=c(0,100),
     xlab='', ylab='',
     xaxt='n',yaxt='n',
     type='n')
text(0.05,95,'F', cex=2)
#swarm locations
swarm_locs<-c(0.25,0.75,1.25)

#ADD MEAN LINES
#calculate means per group
means<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$cpu_setup), mean)
#calculate means per group
cis<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$cpu_setup), Rmisc::CI)
#calculate lower and upper hinges using custom functions
lowerhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$cpu_setup), lower_hinge)
upperhinges<-tapply(studydata_formatted$recall,as.factor(studydata_formatted$cpu_setup), upper_hinge)
notches_matrix<-cbind(unlist(lapply(cis,min))<lowerhinges,unlist(lapply(cis,max))<lowerhinges)
#would it have notches in a box plot? i.e. is confidence interval (range) is more than the inter-quartile range
notches<-rowSums(notches_matrix)
#initiate mean colour column by adding 1 to true/false notches (so false is lty=1)
mean_cols<-notches+1
#add means with appropriate colour to plot
segments(swarm_locs-0.125,means,swarm_locs+0.125,means, lwd=3, col=mean_cols)

#ADD BEESWARM POINTS
beeswarm::beeswarm(studydata_formatted$recall~studydata_formatted$cpu_setup,
                   at=swarm_locs, pwcol=rep('darkgrey',length(studydata_formatted$cpu_setup)),
                   pch = 19, cex=1,
                   xlab = '', ylab='',
                   notch=F, add=T)
axis(1, swarm_locs, levels(studydata_formatted$cpu_setup), cex.axis=1, padj=0)
axis(2,seq(0,100,20),NA,tick = T)
mtext("Computer set-up", side=1, cex=0.7, line=2.5)

#FIGURE G: INITIATE PLOT

plot(studydata_formatted$recall~studydata_formatted$timespent,
     col='darkgrey', pch=19, xlim=c(1,10), las=2,
     xaxt='n',xlab='',
     yaxt='n', ylab='',
     ylim=c(0,100),las=2,
     cex.lab=1.25)
text(1.25,95,'G', cex=5)
#frame plot as is a tested relationship
box(lwd=3)

#MAKE AND PLOT PREDICTIONS FROM A TEST MODEL
model<-lm(recall_asin~timespent, data=studydata_formatted)
#make new x values to predict from
newx <- seq(0,11, length.out=63)
#predict from the test model
preds <- predict(model, newdata = data.frame(timespent=newx), interval = 'confidence')
#back-transform predictions (from arcsine-square root)
preds<-(sin(preds)^2)*100
#fill in area between regression line and confidence interval
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = alpha('grey',0.5), border = NA)
#add fitted regression line
lines(preds[,1]~newx, lwd=3)

#ADD BEESWARM POINTS
points(studydata_formatted$recall~studydata_formatted$timespent,
     col='darkgrey', pch=19, xlim=c(1,10), las=2,
     xaxt='n',xlab='',
     ylim=c(0,100), ylab='',las=2,
     cex.lab=1.25)

axis(1, seq(0,10,1),seq(0,10,1), cex.axis=1.25)
axis(2,seq(0,100,20),NA,tick = T)
mtext('Time spent playing (minutes)', side=1, cex=1, line=2.5)

dev.off()

# MODEL(S): USER ENGAGEMENT --------------------------------------------

#make model
userengagement_model<-lm(studydata_formatted$recall_asin~studydata_formatted$timespent)
#add 1 to the test counter
tests<-tests+1
#summarise model and adjust p-values
summary(userengagement_model)
userengagement_adjp<-p.adjust(summary(userengagement_model)$coefficients[,"Pr(>|t|)"][2], method = 'bonferroni', n = 4)

#save model and adjusted p value
saveRDS(userengagement_model, 'outputs/models/userengagement_model.RDS')
saveRDS(userengagement_adjp, 'outputs/models/userengagement_adjp.RDS')

#plot model
layout.matrix<-matrix(1:4,nrow=2,ncol=2)
layout(layout.matrix)
plot(userengagement_model)

# ANALYSIS OF OPEN-FORM FEEDBACK ------------------------------------------

#done via qualitative method (thematic analysis) so there's no R code for this. Please refer to the manuscript, raw datasheet, and supplementary material if you want to look at the raw responses/results of thematic analysis in more detail

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

grDevices::tiff(paste('outputs/figures/likertplot.tiff', sep=''), res=300, units='in', width=14, height=8)

plot(likert_plot,
     text.size=2.5,
     plot.percents=T, plot.percent.low=F,plot.percent.neutral=F, plot.percent.high=F)

dev.off()

# TESTING FOR DIFFERENCES IN SYMPTOMS LOOKED FOR BY FARMING EXPERIENCE --------

symptomslookedfor_byfarmingexp_sum<-t(aggregate(symptomslookedfor_df, list(studydata_formatted$sheeplamenessexp_prev), function(x){sum(x)}))
#put yes first for plotting purposes
symptomslookedfor_byfarmingexp_sum<-symptomslookedfor_byfarmingexp_sum[,c(2,1)]
colnames(symptomslookedfor_byfarmingexp_sum)<-c('Farming experience','No farming experience')
symptomslookedfor_byfarmingexp_sum<-symptomslookedfor_byfarmingexp_sum[-1,]
rownames_keep<-rownames(symptomslookedfor_byfarmingexp_sum)
symptomslookedfor_byfarmingexp_sum<-apply(symptomslookedfor_byfarmingexp_sum,2,as.numeric)
rownames(symptomslookedfor_byfarmingexp_sum)<-rownames_keep
symptomslookedfor_byfarmingexp_sum<-as.table(as.matrix(symptomslookedfor_byfarmingexp_sum))
#save processed data
saveRDS(symptomslookedfor_byfarmingexp_sum,'outputs/processed_data/symptomslookedfor_byfarmingexp_sum.csv')

#check number of symptoms ticked per participant
colSums(symptomslookedfor_byfarmingexp_sum)
x<-colSums(symptomslookedfor_byfarmingexp_sum)
symptomsVSfarmingexp_chisq<-chisq.test(symptomslookedfor_byfarmingexp_sum)
saveRDS(symptomsVSfarmingexp_chisq,'outputs/models/symptomsVSfarmingexp_chisq.RDS')

dev.off()
grDevices::tiff(paste('outputs/figures/balloonplot.tiff', sep=''), res=300, units='in', width=12, height=6)
gplots::balloonplot(t(symptomslookedfor_byfarmingexp_sum),
            main ="", sorted=T,
            xlab ="", ylab="",
            label = FALSE, show.margins = FALSE)
dev.off()
