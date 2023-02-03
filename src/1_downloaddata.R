# DOWNLOAD DATA FROM OSF --------------------------------------------------

#change this to '2_preanalysis' to run the code from step 3, to '3_end' to download the end result of running the code
whichpoint<-'1_start'

# FROM OSF-ARCHIEVED DATA (TO DOWNLOAD FROM OSF)
my_project <- osfr::osf_ls_files(osfr::osf_retrieve_node("a6qu4"))
data_folder <- my_project[which(my_project$name==whichpoint),]

#download all the folders (inputs and outputs)
osfr::osf_download(osfr::osf_ls_files(data_folder),getwd(), recurse = T, conflicts='overwrite')
#remove the osf project as no longer needed as we have a local copy
rm(my_project, data_folder)

