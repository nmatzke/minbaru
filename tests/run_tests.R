
# UNIT TESTS IN R -- THE BARE MINIMUM
# 
# http://www.johnmyleswhite.com/notebook/2010/08/17/unit-testing-in-r-the-bare-minimum/
# 

library("testthat")

source("/drives/GDrive/__github/minbaru/R/read_bdist_v2.R")

# for cft()
source("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R") 

# for get_numstates_per_char()
source("/drives/GDrive/__github/BEASTmasteR/R/read_nexus_data2_v1.R")

# for list2str_fast() (used in read_nexus_data2_v1.R)
source("/drives/GDrive/__github/BEASTmasteR/R/basics_v1.R")
 
test_dir("tests", reporter = "Summary")

