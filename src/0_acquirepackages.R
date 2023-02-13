# START -------------------------------------------------------------------

#if Rtools package not installed, install. Needed for setting working directory
if('Rtools'%in%rownames(installed.packages())==F){
  install.packages('Rtools')
}

#if rstudioapi package not installed, install. Needed for setting working directory
if('rstudioapi'%in%rownames(installed.packages())==F){
  install.packages('rstudioapi')
}

#make a vector of scripts that we need to run
scripts<-list.files('src/', full.names = T)[-grep('acquirepackages.R',list.files('src/', full.names = T))]
#add the markdown script to this (it's outside of the 'src directory)
scripts<-c(scripts,'lamenessgame_MS.Rmd')

#parse the lines in the scripts to create searchable strings
file_lines<-unlist(sapply(scripts,readLines))

#find packages by looking for ::'s in scripts
pckgs <- lapply(file_lines, function(l) { 
  if(grepl("::", l)){
    gsub(".*?([[:alnum:]\\.]+)::.*","\\1", l) 
  } else {
    return(NULL)
  }
  
})  

#Create a vector of all required packages for this repo's code
required_packages<-unique(unlist(pckgs))
#add R Markdown to the list
required_packages<-c(required_packages,'rmarkdown')

# Create a vector of all packages are already installed
installed_packages<-rownames(installed.packages())
# Create a vector of packages that need to be installed on this computer to run the repo's code
needtoinstall_packages<-required_packages[which(!required_packages%in%installed_packages)]

#For all packages not installed, install them
install.packages(needtoinstall_packages)

#Although 99% of the time sourcing an external function via :: instead of library() is fine, if you load an RDS of an output produced by a package that is not currently loaded, it'll screw up, so best just to load the packages anyway, as below
invisible(lapply(required_packages, library, character.only = TRUE))

#Some links of code from which I derived this script here:
#https://statsandr.com/blog/an-efficient-way-to-install-and-load-r-packages/
#https://stackoverflow.com/questions/59725527/list-all-the-packages-required-in-a-script-assuming-packagefunction-in-r