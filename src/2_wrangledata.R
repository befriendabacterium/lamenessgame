# LOAD NECESSARY PACKAGES -------------------------------------------------

library(scales)
library(dplyr)
library(ggplot2)

# READ IN RAW DATA FROM MS FORMS --------------------------------------------------------

#read in raw study data
studydata_raw<-readxl::read_excel('inputs/lamenessstudydata_010721.xlsx')
# #detach(tibble)
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

#remove not needed dataframe
rm(farmingexp_df)

# PARSE FARMING EXPERIENCE/ROLES (NOT USED) -------------------------------------------------------------

# roles<-c(
#   'Farmer',
#   'Stockman/woman/person',
#   'Veterinarian',
#   'Other')
# 
# roles_df<-data.frame(matrix(nrow=nrow(studydata_formatted), ncol=length(roles)))
# colnames(roles_df)<-roles
# for (r in 1:length(roles)){
#   print(r)
#   temp<-grepl(roles[r],studydata_formatted$sheeproles_type)  
#   roles_df[,r]<-temp
# }

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

#write symptoms looked for df (easier to use for plotting)
write.csv(symptomslookedfor_df,'outputs/processed_data/symptomslookedfor_df.csv', row.names = F)

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


# MAKE CONTINGENCY TABLE OF LAMENESS SIGNS LOOKED FOR BY FARMING EXPERIENCE -----------------------------------------------------------------------

symptomslookedfor_byfarmingexp_sum<-t(aggregate(symptomslookedfor_df, list(studydata_formatted$farmingexp_YN), function(x){sum(x)}))
#put yes first for plotting purposes
symptomslookedfor_byfarmingexp_sum<-symptomslookedfor_byfarmingexp_sum[,c(2,1)]
#rename columns more clearly
colnames(symptomslookedfor_byfarmingexp_sum)<-c('Farming experience','No farming experience')
#remove first row with column names (redundant)
symptomslookedfor_byfarmingexp_sum<-symptomslookedfor_byfarmingexp_sum[-1,]
#make a vector of row names
rownames_keep<-rownames(symptomslookedfor_byfarmingexp_sum)
#coerce columns to numeric
symptomslookedfor_byfarmingexp_sum<-apply(symptomslookedfor_byfarmingexp_sum,2,as.numeric)
#put row names back (they're lost in coercion)
rownames(symptomslookedfor_byfarmingexp_sum)<-rownames_keep
#coerce the dataframe to table (needed for baloon plot)
symptomslookedfor_byfarmingexp_sum<-as.table(as.matrix(symptomslookedfor_byfarmingexp_sum))
#check it
symptomslookedfor_byfarmingexp_sum
#save processed data
write.csv(symptomslookedfor_byfarmingexp_sum,'outputs/processed_data/symptomslookedfor_byfarmingexp_sum.csv')
saveRDS(symptomslookedfor_byfarmingexp_sum,'outputs/processed_data/symptomslookedfor_byfarmingexp_sum.RDS')

# FORMATTING: LIKERT DATA --------------------------------------------------------------

# subset out and tidy likert columns
firstcol<-which(colnames(studydata_raw)=="Points - How strongly do you agree with the following statements?")
#cut out section of df where the likert data starts
likertdata_formatted<-studydata_raw[,firstcol:ncol(studydata_raw)]
#remove unwanted columns
likertdata_formatted<-likertdata_formatted[,-grep("Points|Feedback|Please share any general feedback/main thoughts after playing the game below|Finally, how did you find out about this study?",colnames(likertdata_formatted))]
#remove weird bits from some of the statements
colnames(likertdata_formatted)<-gsub("\r\n\r\n","",colnames(likertdata_formatted))
#convert columns to factor
likertdata_formatted<-as.data.frame(apply(likertdata_formatted,2,factor))

write.csv(likertdata_formatted,'outputs/processed_data/likertdata_formatted.csv', row.names = F)
saveRDS(likertdata_formatted,'outputs/processed_data/likertdata_formatted.RDS')

#vector of likert categories
likert_categories_ordered<-c("Strongly disagree", "Disagree", "Neutral", "Agree", "Strongly agree")


# REDUNDANT FORMATTING CODE
# #count up the responses
# likertdata_formatted<-apply(likertdata,2,function(x){plyr::count(x)})
# likertdata_formatted
# 
# #rename columns with loop
# for (l in 1:length(likertdata_formatted)){
#   colnames(likertdata_formatted[[l]])[1]<-"Statement"
#   colnames(likertdata_formatted[[l]])[2]<-names(likertdata_formatted[l])
# }
# 
# #function for merge
# my_merge <- function(df1, df2){                                # Create own merging function
#   merge(df1, df2, by = "Statement", all=T)
# }
# 
# #apply function to list
# likertdata_formatted<-Reduce(my_merge, likertdata_formatted)                                    # Apply Reduce to own function
# 
# #transpose dataframe
# likertdata_formatted<-t(likertdata_formatted)
# 
# #set colnames
# colnames(likertdata_formatted)<-likertdata_formatted[1,]
# 
# #remove first row (colnames)
# likertdata_formatted<-as.data.frame(likertdata_formatted[-1,])
# 
# #remove NA col (did not answer question cos non-farmer or missed it)
# likertdata_formatted<-likertdata_formatted[,!is.na(colnames(likertdata_formatted))]
# 
# #vector of likert categories
# likert_categories_ordered<-c("Strongly disagree", "Disagree", "Neutral", "Agree", "Strongly agree")
# 
# #re-order columns
# likertdata_formatted<-likertdata_formatted[,match(likert_categories_ordered,colnames(likertdata_formatted))]
# 
# likertdata_formatted
# 
# #write csv for the formatted likert data
# write.csv(likertdata_formatted, 'outputs/processed_data/likertdata_formatted.csv')

#ugly way of identifying likert columns in study data
likert_cols<-c()
for (l in 1:ncol(studydata_formatted)){
  
  current<-agrep(gsub(' ','.',colnames(likertdata_formatted)[l]),colnames(studydata_formatted))
  #provided there's a match, bind it
  if(!any(is.na(current))){
    likert_cols<-c(likert_cols,current)
  }
  
  
}

# TIDYING UP --------------------------------------------------------------

#remove reundant objects created through the pipeline 
rm(studydata_raw, newcolnames, colstoremove, current, firstcol, l, s, temp, rownames_keep)

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

#write the RDS (keeps level ordering for plotting)
saveRDS(studydata_formatted,'outputs/processed_data/studydata_formatted.RDS')

#clear whole workspace
rm(list = ls())
