# DESCRIPTION -------------------------------------------------------------

# Run this script to wrangle the data into a format more suited to analysis and plotting.

# READ IN RAW DATA FROM MS FORMS --------------------------------------------------------

#read in raw study data
studydata_raw<-readxl::read_excel('inputs/lamenessstudydata_010721.xlsx')
#read in R data file with the desired new column names
newcolnames<-readRDS('inputs/newcolnames.RDS')

# FORMATTING: PARTICIPANT PERFORMANCE DATA --------------------------------------------------------------
#Here we do some initial, basic formatting on the dataframe

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
#Here we arcsine-transform the outcome data, which is expressed in percentages, to make it more amenable to parametric analysis

#create a new column of arcsin-transformed accuracy scores
studydata_formatted$accuracy_asin<-asin(sqrt(0.01*studydata_formatted$accuracy))
#arcsin-transformed recall
studydata_formatted$recall_asin<-asin(sqrt(0.01*studydata_formatted$recall))

# PARSE FARMING/SHEEP FARMING EXPERIENCE  ---------------------------------------------------------------------
# Here we do some tidying on the farming/sheep farming experience data to make it easier to work with/plot sensibly

#make a vector of possible roles in which people may have worked with sheep
roles<-c(
  'Farmer',
  'Stockman/woman/person',
  'Veterinarian',
  'Other')

#create an empty dataframe with 63 rows (participants) and 4 cols (roles) 
roles_df<-data.frame(matrix(nrow=nrow(studydata_formatted), ncol=length(roles)))
#set col names to roles
colnames(roles_df)<-roles
#loop over each role and create a binary vector (column) that expresses whether or not the participant had worked in that role
for (r in 1:length(roles)){
  #use grep to make a temporary vector to say whether or not (TRUE/FALSE) the role appeared in each participants' answer string (roles are ; seperated)
  temp<-grepl(roles[r],studydata_formatted$sheeproles_type)
  #write this temporary vector to the column for the current role
  roles_df[,r]<-temp
}

#add a new column, 'Any', which expresses whether or not the participant had worked with sheep in any role
roles_df<-cbind(roles_df,Any=rowSums(roles_df)>=1)

#make row names participant IDs (participant IDs start at 8, not 1)
rownames(roles_df)<-studydata_formatted$ID
#change colnames to the questiontype_subtype format 
colnames(roles_df)<-c('sheeproles_farmer','sheeproles_stockperson','sheeproles_vet', 'sheeproles_other', 'sheeproles_any')

#write csv and RDS with row names
write.csv(roles_df, 'outputs/processed_data/roles_df.csv', row.names = T)
saveRDS(roles_df, 'outputs/processed_data/roles_df.RDS')

#add to main dataframe
studydata_formatted<-cbind(studydata_formatted,roles_df)

#remove the original sheeproles_type coluumn from the dataframe to remove redundancy in the dataframe
studydata_formatted<-studydata_formatted[,-which(colnames(studydata_formatted)=='sheeproles_type')]

#reorder levels of lameness prevalence so that they plot/model in order of increasing prevalence
studydata_formatted$sheeplamenessexp_prev<-factor(studydata_formatted$sheeplamenessexp_prev, levels=c("Under 2%", "Between 2 and 5%", "Between 5 and 10%","Over 10%"))

# PARSE LAMENESS SIGNS (SYMPTOMS) -------------------------------------------------------------

#create a vector of symptoms
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

#create an empty dataframe with 63 rows (participants) and 9 cols (symptoms) 
symptomslookedfor_df<-data.frame(matrix(nrow=nrow(studydata_formatted), ncol=length(symptoms)))
#set col names to symptoms
colnames(symptomslookedfor_df)<-symptoms
#loop over each symptom and create a binary vector (column) that expresses whether or not the participant looked for that symptom
for (s in 1:length(symptoms)){
  #use grep to make a temporary vector to say whether or not (TRUE/FALSE) the symptom appeared in each participants' answer string (symptoms are ; seperated)
  temp<-grepl(symptoms[s],studydata_formatted$symptoms_lookedfor)  
  #write this temporary vector to the column for the current symptom
  symptomslookedfor_df[,s]<-temp
}

#shorten the names of some symptoms so that fit on plot x-axis labels easier
colnames(symptomslookedfor_df)[which(colnames(symptomslookedfor_df)%in%c("Pair of legs which were moving at different speeds","Not weight bearing on affected leg when standing","Not weight bearing on affected leg when walking"))]<-c("Pair of legs moving at different speeds","Not weight bearing on affected leg (standing)", "Not weight bearing on affected leg (walking)")
#do same for symptoms vector
symptoms[which(symptoms%in%c("Pair of legs which were moving at different speeds","Not weight bearing on affected leg when standing","Not weight bearing on affected leg when walking"))]<-c("Pair of legs moving at different speeds","Not weight bearing on affected leg (standing)", "Not weight bearing on affected leg (walking)")

#make row names participant IDs (participant IDs start at 8, not 1)
rownames(symptomslookedfor_df)<-studydata_formatted$ID

#write symptoms looked for df (easier to use for plotting)
saveRDS(symptomslookedfor_df,'outputs/processed_data/symptomslookedfor_df.RDS')
write.csv(symptomslookedfor_df,'outputs/processed_data/symptomslookedfor_df.csv', row.names = T)

#change colnames of dataframe to shorthand (we'll use the symptoms vector to label)
colnames(symptomslookedfor_df)<-c('UP','SS','LS','NH','WBS','WBW','RM','SW','O')
#add to main dataframe
studydata_formatted<-cbind(studydata_formatted,symptomslookedfor_df)

#remove the original symptoms_lookedfor column to remove redundancy in the dataframe
studydata_formatted<-studydata_formatted[,-which(colnames(studydata_formatted)=='symptoms_lookedfor')]

# PARSE USER ENGAGEMENT ---------------------------------------------------

#calculate and add timespent
studydata_formatted$timespent<-(600-studydata_formatted$timeremaining_s)/60
#remove the original timeremaining_s column to remove redundancy in the dataframe
studydata_formatted<-studydata_formatted[,-which(colnames(studydata_formatted)=='timeremaining_s')]

#change cpu setup levels
#those strings containing 'mouse' become 'mouse'
studydata_formatted$cpu_setup[agrep('Mouse',studydata_formatted$cpu_setup)]<-'Mouse'
#those strings containing 'trackpad' become 'trackpad'
studydata_formatted$cpu_setup[agrep('Track-pad',studydata_formatted$cpu_setup)]<-'Trackpad'
#those remaining (i.e. only containing laptop or desktop computer) get 'unknown'
studydata_formatted$cpu_setup[agrep('Laptop or Desktop computer;',studydata_formatted$cpu_setup)]<-'Unknown'
#reorder levels for plotting
studydata_formatted$cpu_setup<-factor(studydata_formatted$cpu_setup, levels=c('Unknown','Trackpad','Mouse'))

#change all columns to factor
studydata_formatted$tutorial<-as.factor(studydata_formatted$tutorial)
studydata_formatted$cpu_setup<-as.factor(studydata_formatted$cpu_setup)
studydata_formatted$controlsprobs_YN<-as.factor(studydata_formatted$controlsprobs_YN)
studydata_formatted$strategy_type<-as.factor(studydata_formatted$strategy_type)
studydata_formatted$moving_type<-as.factor(studydata_formatted$moving_type)

#shorten level names for strategy type
levels(studydata_formatted$strategy_type)<-c('Up-close','Zoomer','Other')
#shorten level names for moving type
levels(studydata_formatted$moving_type)<-c('Other','Randomly','Semi-randomly','Linear')
#reorder levels for moving type for plotting
studydata_formatted$moving_type<-factor(studydata_formatted$moving_type, levels=c('Linear','Semi-randomly','Randomly','Other'))

# MAKE CONTINGENCY TABLE OF LAMENESS SIGNS LOOKED FOR BY FARMING EXPERIENCE -----------------------------------------------------------------------

#count up the number of symptoms (row sums of symptoms looked for) by farming experience (create contigency table)
symptomslookedfor_byfarmingexp_sum<-t(aggregate(symptomslookedfor_df, list(studydata_formatted$farmingexp_YN), function(x){sum(x)}))
#put yes (farming experience) first for plotting purposes
symptomslookedfor_byfarmingexp_sum<-symptomslookedfor_byfarmingexp_sum[,c(2,1)]
#rename columns more clearly
colnames(symptomslookedfor_byfarmingexp_sum)<-c('Farming experience','No farming experience')
#remove first row with column names (redundant)
symptomslookedfor_byfarmingexp_sum<-symptomslookedfor_byfarmingexp_sum[-1,]
#rename rows more clearly
rownames(symptomslookedfor_byfarmingexp_sum)<-symptoms
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
firstcol<-which(colnames(studydata_formatted)=="The.game.is.a.realistic.representation.of.recognising.sheep.lameness.in.the.field..")
#cut out section of df where the likert data starts + 20 questions + open form feedback question
likertdata_formatted<-studydata_formatted[,seq(firstcol, length.out=21)]
#identify open form feedback question
openformq<-which(colnames(likertdata_formatted)=='Please.share.any.general.feedback.main.thoughts.after.playing.the.game.below.')
#remove it
likertdata_formatted<-likertdata_formatted[,-openformq]
#replace periods in likert data with spaces
colnames(likertdata_formatted)<-gsub("\\."," ",colnames(likertdata_formatted))
#convert columns to factor
likertdata_formatted<-as.data.frame(apply(likertdata_formatted,2,factor))
#make row names to participant ID
rownames(likertdata_formatted)<-studydata_formatted$ID

#write csv including row names
write.csv(likertdata_formatted,'outputs/processed_data/likertdata_formatted.csv', row.names = T)
#write RDS (better for plotting)
saveRDS(likertdata_formatted,'outputs/processed_data/likertdata_formatted.RDS')

#vector of likert categories
likert_categories_ordered<-c("Strongly disagree", "Disagree", "Neutral", "Agree", "Strongly agree")

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
  dplyr::relocate(studydata_formatted, 
           c('accuracy','accuracy_asin','recall','recall_asin'),
           .after = colnames(studydata_formatted)[ncol(studydata_formatted)])

#move farming experience columns to end
studydata_formatted<-
  dplyr::relocate(studydata_formatted, 
           c('farmingexp_YN', 'sheeproles_farmer','sheeproles_stockperson','sheeproles_vet','sheeproles_other', 'sheeproles_any','sheepexp_years', 'sheeproles_details','sheeplamenessexp_prev'),
           .after = colnames(studydata_formatted)[ncol(studydata_formatted)])

#move signs looked for columns to end
studydata_formatted<-
  dplyr::relocate(studydata_formatted, 
           c('UP','SS','LS','NH','WBS','WBW','RM','SW','O', 'symptoms_otherdetails'),
           .after = colnames(studydata_formatted)[ncol(studydata_formatted)])

#move user engagement columns to end
studydata_formatted<-
  dplyr::relocate(studydata_formatted, 
           c('timesplayed','tutorial','cpu_setup','controlsprobs_YN','controlsprobs_descrip','strategy_type','strategy_otherdetails','moving_type','moving_otherdetails', 'timespent'),
           .after = colnames(studydata_formatted)[ncol(studydata_formatted)])

#move feedback column to end
studydata_formatted<-
  dplyr::relocate(studydata_formatted, 
           c('Please.share.any.general.feedback.main.thoughts.after.playing.the.game.below.'),
           .after = colnames(studydata_formatted)[ncol(studydata_formatted)])

colnames(studydata_formatted)

#write the csv
write.csv(studydata_formatted,'outputs/processed_data/studydata_formatted.csv', row.names = F)

#write the RDS (keeps level ordering for plotting)
saveRDS(studydata_formatted,'outputs/processed_data/studydata_formatted.RDS')

#clear whole workspace
rm(list = ls())
