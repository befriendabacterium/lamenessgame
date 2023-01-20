# LOAD NECESSARY PACKAGES -------------------------------------------------

#install.packages('osfr')
library(osfr)

# FROM OSF-ARCHIEVED DATA (TO DOWNLOAD FROM OSF)
my_project <- osfr::osf_ls_files(osfr::osf_retrieve_node("a6qu4"))

#download all the folders (inputs and outputs)
osfr::osf_download(my_project,getwd(), recurse = T, conflicts='overwrite')
#remove the osf project as no longer needed as we have a local copy
rm(my_project)
