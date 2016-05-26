# Tests for the "minbaru" package
library(testthat)

# Code to source
# For minbaru functions 
source("/drives/GDrive/__github/minbaru/R/read_bdist_v2.R")
source("/drives/GDrive/__github/minbaru/R/fast_bdist_v1.R")
source("/drives/GDrive/__github/minbaru/R/bdist_relevance_v1.R")

# for cft()
source("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_generics_v1.R") 

# for get_numstates_per_char()
source("/drives/GDrive/__github/BEASTmasteR/R/read_nexus_data2_v1.R")

# for list2str_fast() (used in read_nexus_data2_v1.R)
source("/drives/GDrive/__github/BEASTmasteR/R/basics_v1.R")


#######################################################
# Set the working directory
#######################################################
wd = "/drives/GDrive/__github/minbaru/inst/extdata/example_scripts/"
setwd(wd)


cat("TESTING THAT STANDARD (LITERAL) BDIST CALCULATIONS WORK...ex5_Dembo_Matzke_etal_Character_matrix_v3_simp.nex")
# Read in an example file, correcting it to standard
# character representation (as e.g. from NEXUS)
# Have to have a cutoff of e.g. 0.95 to match bdist.pl defaults
nexfn = "/drives/GDrive/__github/minbaru/inst/extdata/examples/ex5_Dembo_Matzke_etal_Character_matrix_v3_simp.nex"

chars_list = read_nexus_data2(file=nexfn, check_ambig_chars=TRUE, convert_ambiguous_to=NULL, printall="short", convert_ambiguous_to_IUPAC=FALSE) 

chardf = fix_NEXUS_charslist_for_bdist(chars_list=chars_list, remove_spaces=TRUE, alphabetize=TRUE)

# Extract and subset by completeness (should be done manually)
dim(chardf)
relevance_res = bdist_relevance(chardf=chardf, min_char_relevance=0.75, min_taxic_relevance=0.4, alphabetize=TRUE, subset_data=TRUE)
chardf = relevance_res$subboth_chardf
dim(chardf)

runslow = FALSE
resfn = "ex5_Dembo_Matzke_etal_Character_matrix_v3_simp_result.Rdata"
if (runslow)
	{
	res = fast_bdist(chardf, alphabetize=TRUE, grepl_w_brackets=TRUE, num_bootstraps=100, printflag=TRUE)
	# Loads to res
	save(res, file=resfn)
	} else {
	# Loads to res
	load(file=resfn)
	} # END if (runslow)
res$bootstrap_freqs_pLT05

bdists_minbaru = round(res$bdist_dij, digits=3)
bdists_minbaru


#######################################################
# Make some plots
#######################################################
Wood_2016_Fig1_order = c("P_robustus", "P_boisei", "H_rudolfensis", "H_habilis", "H_neanderthalensis", "H_ergaster", "H_sapiens", "H_heidelbergensis", "H_erectus", "A_africanus", "Pan_troglodytes",  "Gorilla_gorilla", "A_afarensis")
Wood_2016_Fig1_order_nums = match(x=Wood_2016_Fig1_order, table=rownames(res$bdist_cor_Ps))
new_order=Wood_2016_Fig1_order_nums

# Different orderings of taxa
bds_bootstrap_plot(res, taxon_to_order_on=NULL, new_order="heatmap")
bds_bootstrap_plot(res, taxon_to_order_on="H_sapiens", new_order=NULL)
bds_bootstrap_plot(res, taxon_to_order_on=NULL, new_order=Wood_2016_Fig1_order_nums)
heatmap(t(res$cormat))

source("/drives/GDrive/__github/minbaru/R/fast_bdist_v1.R")
bds_bootstrap_plot(res, taxon_to_order_on=NULL, new_order="heatmap", toplabels=TRUE)
bds_bootstrap_plot(res, taxon_to_order_on="H_sapiens", new_order=NULL, toplabels=TRUE)
bds_bootstrap_plot(res, taxon_to_order_on=NULL, new_order=Wood_2016_Fig1_order_nums, toplabels=TRUE)


#######################################################
# 3-dimensional Multidimensional Scaling (and plotting)
#######################################################
# metaMDS from the vegan package
library(vegan)
dists = as.dist(res$bdist_dij)
mds_result = metaMDS(comm=dists, k=3)
plot(mds_result, type="t")



#######################################################
# Principle components
#######################################################
# http://planspace.org/2013/02/03/pca-3d-visualization-and-clustering-in-r/
species = names(res$bdist_dij)
genera = c("Australopithecus", "Australopithecus", "Gorilla", "Homo", "Homo", "Homo", "Homo", "Homo", "Homo", "Homo", "Paranthropus", "Paranthropus", "Pan")
colors = c("blue", "blue", "green3", "red", "red", "red", "red", "red", "red", "red", "black", "black", "green")

pc <- princomp(x=dists, cor=TRUE, scores=TRUE)
summary(pc)
plot(pc,type="lines")
x = biplot(x=pc, var.axes=FALSE)
ev = eigen(dists)
plot(x=ev$vectors[,1], y=ev$vectors[,2])
text(x=ev$vectors[,1], y=ev$vectors[,2], labels=species, cex=0.5, col=colors)
library(rgl)
plot3d(pc$scores[,1:3], col=colors, size=15)

