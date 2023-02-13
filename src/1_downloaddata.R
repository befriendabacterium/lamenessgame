# DESCRIPTION -------------------------------------------------------------

# Run this script to download the raw input data from OSF (https://osf.io/a6qu4/), plus some outputs that was produced manually (rather than via code).

# DOWNLOAD DATA FROM OSF --------------------------------------------------

# FROM OSF-ARCHIEVED DATA (TO DOWNLOAD FROM OSF)
my_project <- osfr::osf_ls_files(osfr::osf_retrieve_node("a6qu4"))

#download all the folders (inputs and outputs)
osfr::osf_download(my_project,getwd(), recurse = T, conflicts='overwrite')
#remove the osf project as no longer needed as we have a local copy
rm(my_project)

