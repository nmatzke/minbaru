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
wd = "/drives/GDrive/__github/minbaru/inst/extdata/examples/"
setwd(wd)

# Check that math works
cat("\n\n")
cat("=========================BEGINNING TESTS=========================\n\n")
cat("TESTING THAT MATH WORKS...")
expect_that(1 ^ 1, equals(1))
expect_that(2 ^ 2, equals(4))
expect_that(2 + 2 == 4, is_true())
expect_that(2 == 1, is_false())

# Classes, printing
expect_that(1, is_a('numeric'))
expect_that(print('Hello World!'), prints_text('Hello World!'))
expect_that(log('a'), throws_error())

# Speed (deprecated)
expect_that(factorial(16), takes_less_than(1))
cat("...DONE\n")


cat("TESTING THAT READING BDIST FILES WORKS (correcting data as in NEXUS)...")
# Read in an example file, correcting it to standard
# character representation (as e.g. from NEXUS)
fn = "/drives/GDrive/__github/minbaru/inst/extdata/examples/ex1_from_bdist_online.txt"
#min_char_relevance=0; min_taxic_relevance=0; correct_numeric_above_10=TRUE; correct_w_get_numstates_per_char=TRUE; return_list=TRUE; alphabetize=TRUE; printflag=TRUE
readres = read_bdist_to_chardf(fn=fn, min_char_relevance=0, min_taxic_relevance=0, correct_numeric_above_10=TRUE, correct_w_get_numstates_per_char=TRUE, return_list=TRUE, alphabetize=TRUE, printflag=TRUE)
readres

readres_txt = dput(readres$chardf)
correct_txt = structure(list(daphne = c("0", "0", "2", "0", "0", "1", "23"), 
    freddie = c("0", "01", "1", "0", "0", "1", "01"), scooby = c("1", 
    "0", "2", "0", "1", "1", "1"), shaggy = c("0", "1", "0", 
    "2", "2", "0", "01"), velma = c("0", "01", "2", "1", "12", 
    "0", "?")), .Names = c("daphne", "freddie", "scooby", "shaggy", 
"velma"), row.names = c(NA, 7L), class = "data.frame")

# Run the test
expect_identical(object=readres_txt, expected=correct_txt)
cat("...DONE\n")


cat("TESTING THAT READING BDIST FILES WORKS (NOT correcting data as in NEXUS)...")
# Read in an example file, leaving it as original 
# from the example (e.g., "10" different from "01"),
# and 12, 13, 14 is a valid character pattern
fn = "/drives/GDrive/__github/minbaru/inst/extdata/examples/ex1_from_bdist_online.txt"
readres = read_bdist_to_chardf(fn=fn, min_char_relevance=0, min_taxic_relevance=0, correct_numeric_above_10=FALSE, correct_w_get_numstates_per_char=FALSE, return_list=TRUE, alphabetize=TRUE, printflag=TRUE)
readres_txt = dput(readres$chardf)
correct_txt = structure(list(daphne = c("0", "0", "2", "12", "0", "1", "(23)"
), freddie = c("0", "(10)", "1", "12", "0", "1", "(01)"), scooby = c("1", 
"0", "2", "12", "1", "1", "1"), shaggy = c("0", "1", "0", "14", 
"2", "0", "(01)"), velma = c("0", "(01)", "2", "13", "(12)", 
"0", "?")), .Names = c("daphne", "freddie", "scooby", "shaggy", 
"velma"), row.names = c(NA, 7L), class = "data.frame")

# Run the test
expect_identical(object=readres_txt, expected=correct_txt)
cat("...DONE\n")




# Read in an example file, leaving it as original 
# from the example (e.g., "10" different from "01"),
# except that 12, 13, 14 are converted to B, C, D
# Then do "literal" (standard) bdist

cat("TESTING THAT STANDARD (LITERAL) BDIST CALCULATIONS WORK...ex1_from_bdist_online.txt")
# Read in an example file, correcting it to standard
# character representation (as e.g. from NEXUS)
# Have to have a cutoff of e.g. 0.95 to match bdist.pl defaults
fn = "/drives/GDrive/__github/minbaru/inst/extdata/examples/ex1_from_bdist_online.txt"
readres = read_bdist_to_chardf(fn=fn, min_char_relevance=0.95, min_taxic_relevance=0, correct_numeric_above_10=TRUE, correct_w_get_numstates_per_char=TRUE, return_list=TRUE, alphabetize=TRUE, printflag=TRUE, subset_data=FALSE)

# Extract and subset by completeness (should be done manually)
chardf = readres$chardf[,readres$taxic_relevance$relTF]
chardf = chardf[readres$stats_table$include_TF,]
#chardf = readres$raw_chardf

# chardf
# keep_taxa_TF = readres$taxic_relevance$TaxicRel_raw >= 0.95
# chardf = readres$chardf[,keep_taxa_TF]
# chardf
res = fast_bdist(chardf, alphabetize=TRUE)
res = fast_bdist(chardf, alphabetize=TRUE, num_bootstraps=100, printflag=TRUE)
res
bdists_minbaru = round(res$bdist_dij, digits=3)

bdist_pl_output_txt='	daphne	freddie	scooby	shaggy	velma
daphne	0.000	0.167	0.333	0.833	0.500
freddie	0.167	0.000	0.500	0.667	0.667
scooby	0.333	0.500	0.000	1.000	0.500
shaggy	0.833	0.667	1.000	0.000	0.333
velma	0.500	0.667	0.500	0.333	0.000'
textcon = textConnection(bdist_pl_output_txt)
bdists_bdist_pl = read.table(textcon)
close(textcon)
bdists_bdist_pl

bdists_minbaru - bdists_bdist_pl

# Run the test
expect_identical(object=bdists_minbaru, expected=bdists_bdist_pl)
cat("...DONE\n")








# Read in an example file, from the CRSQ article:
# Figure 1A of:
# Robinson, D. Ashley; Cavanaugh, David P. (1998). 
#  A Quantitative Approach to Baraminology With Examples from the 
#  Catarrhine Primates. Creation Research Society Quarterly, 34(4), 
#  196-208.

cat("TESTING THAT STANDARD (LITERAL) BDIST CALCULATIONS WORK...ex3b_CRSQ_1998_p198_Fig1A_data.txt")
# Read in an example file, correcting it to standard
# character representation (as e.g. from NEXUS)
# Have to have a cutoff of e.g. 0.95 to match bdist.pl defaults
fn = "/drives/GDrive/__github/minbaru/inst/extdata/examples/ex3b_CRSQ_1998_p198_Fig1A_data.txt"
readres = read_bdist_to_chardf(fn=fn, min_char_relevance=0, min_taxic_relevance=0, correct_numeric_above_10=TRUE, correct_w_get_numstates_per_char=TRUE, return_list=TRUE, alphabetize=TRUE, printflag=TRUE, subset_data=FALSE)

# Extract and subset by completeness (should be done manually)
chardf = readres$chardf[,readres$taxic_relevance$relTF]
chardf = chardf[readres$stats_table$include_TF,]
#chardf = readres$raw_chardf

# chardf
# keep_taxa_TF = readres$taxic_relevance$TaxicRel_raw >= 0.95
# chardf = readres$chardf[,keep_taxa_TF]
# chardf
res = fast_bdist(chardf, alphabetize=TRUE)
res = fast_bdist(chardf, alphabetize=TRUE, num_bootstraps=100, printflag=TRUE)
res
bdists_minbaru = round(res$bdist_dij, digits=3)

bdistmds_output_txt='	sp1	sp2	sp3	sp4	sp5
sp1	0.000	0.400	0.600	0.750	0.800
sp2	0.400	0.000	0.800	1.000	0.600
sp3	0.600	0.800	0.000	0.750	1.000
sp4	0.750	1.000	0.750	0.000	0.750
sp5	0.800	0.600	1.000	0.750	0.000'
textcon = textConnection(bdistmds_output_txt)
bdists_bdistmds = read.table(textcon)
close(textcon)
bdists_bdistmds

bdists_minbaru - bdists_bdistmds

# Run the test
expect_identical(object=bdists_minbaru, expected=bdists_bdistmds)
cat("...DONE\n")













cat("TESTING THAT STANDARD (LITERAL) BDIST CALCULATIONS WORK...ex4_Dembo_Matzke_etal_for_bdist_subset.txt")
# Read in an example file, correcting it to standard
# character representation (as e.g. from NEXUS)
# Have to have a cutoff of e.g. 0.95 to match bdist.pl defaults
fn = "/drives/GDrive/__github/minbaru/inst/extdata/examples/ex4_Dembo_Matzke_etal_for_bdist_subset.txt"
readres = read_bdist_to_chardf(fn=fn, min_char_relevance=0.75, min_taxic_relevance=0.0, correct_numeric_above_10=TRUE, correct_w_get_numstates_per_char=TRUE, return_list=TRUE, alphabetize=TRUE, printflag=TRUE, subset_data=FALSE)

# Extract and subset by completeness (should be done manually)
dim(readres$chardf)
bdist_chardf_orig = readres$chardf
dim(bdist_chardf_orig)

chardf = readres$chardf[,readres$taxic_relevance$relTF]
chardf = chardf[readres$stats_table$include_TF,]
dim(chardf)
bdist_chardf = chardf
res = fast_bdist(chardf, alphabetize=TRUE, grepl_w_brackets=TRUE)

runslow = FALSE
resfn = "ex4_Dembo_Matzke_etal_for_bdist_subset_result.Rdata"
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

bdist_pl_output_txt='	A_afarensis	A_africanus	Gorilla_gorilla	H_erectus	H_ergaster	H_habilis	H_heidelbergensis	H_neanderthalensis	H_rudolfensis	H_sapiens	P_boisei	P_robustus	Pan_troglodytes
A_afarensis	0.000	0.287	0.328	0.348	0.346	0.346	0.346	0.383	0.369	0.434	0.552	0.553	0.304
A_africanus	0.287	0.000	0.371	0.242	0.284	0.239	0.304	0.292	0.261	0.333	0.431	0.412	0.296
Gorilla_gorilla	0.328	0.371	0.000	0.360	0.428	0.392	0.312	0.340	0.506	0.468	0.456	0.619	0.155
H_erectus	0.348	0.242	0.360	0.000	0.111	0.176	0.063	0.162	0.120	0.092	0.391	0.373	0.303
H_ergaster	0.346	0.284	0.428	0.111	0.000	0.147	0.134	0.198	0.130	0.098	0.431	0.370	0.352
H_habilis	0.346	0.239	0.392	0.176	0.147	0.000	0.205	0.274	0.130	0.207	0.431	0.361	0.352
H_heidelbergensis	0.346	0.304	0.312	0.063	0.134	0.205	0.000	0.105	0.167	0.054	0.291	0.412	0.243
H_neanderthalensis	0.383	0.292	0.340	0.162	0.198	0.274	0.105	0.000	0.236	0.151	0.346	0.379	0.305
H_rudolfensis	0.369	0.261	0.506	0.120	0.130	0.130	0.167	0.236	0.000	0.196	0.447	0.386	0.419
H_sapiens	0.434	0.333	0.468	0.092	0.098	0.207	0.054	0.151	0.196	0.000	0.390	0.417	0.410
P_boisei	0.552	0.431	0.456	0.391	0.431	0.431	0.291	0.346	0.447	0.390	0.000	0.151	0.515
P_robustus	0.553	0.412	0.619	0.373	0.370	0.361	0.412	0.379	0.386	0.417	0.151	0.000	0.644
Pan_troglodytes	0.304	0.296	0.155	0.303	0.352	0.352	0.243	0.305	0.419	0.410	0.515	0.644	0.000'
textcon = textConnection(bdist_pl_output_txt)
bdists_bdist_pl = read.table(textcon)
close(textcon)
bdists_bdist_pl

bdists_minbaru - bdists_bdist_pl

# Run the test
expect_identical(object=bdists_minbaru, expected=bdists_bdist_pl)
cat("...DONE\n")


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
chardf_orig = relevance_res$subtaxa_chardf
dim(chardf_orig)
chardf = relevance_res$subboth_chardf
dim(chardf)

# Check for disagreements -- just re-codings I think
TF = bdist_chardf == chardf
cbind(bdist_chardf[TF==FALSE], chardf[TF==FALSE])

res = fast_bdist(chardf, alphabetize=TRUE, grepl_w_brackets=TRUE)

bdists_minbaru = round(res$bdist_dij, digits=3)
bdists_minbaru

bdist_pl_output_txt='	A_afarensis	A_africanus	Gorilla_gorilla	H_erectus	H_ergaster	H_habilis	H_heidelbergensis	H_neanderthalensis	H_rudolfensis	H_sapiens	P_boisei	P_robustus	Pan_troglodytes
A_afarensis	0.000	0.287	0.328	0.348	0.346	0.346	0.346	0.383	0.369	0.434	0.552	0.553	0.304
A_africanus	0.287	0.000	0.371	0.242	0.284	0.239	0.304	0.292	0.261	0.333	0.431	0.412	0.296
Gorilla_gorilla	0.328	0.371	0.000	0.360	0.428	0.392	0.312	0.340	0.506	0.468	0.456	0.619	0.155
H_erectus	0.348	0.242	0.360	0.000	0.111	0.176	0.063	0.162	0.120	0.092	0.391	0.373	0.303
H_ergaster	0.346	0.284	0.428	0.111	0.000	0.147	0.134	0.198	0.130	0.098	0.431	0.370	0.352
H_habilis	0.346	0.239	0.392	0.176	0.147	0.000	0.205	0.274	0.130	0.207	0.431	0.361	0.352
H_heidelbergensis	0.346	0.304	0.312	0.063	0.134	0.205	0.000	0.105	0.167	0.054	0.291	0.412	0.243
H_neanderthalensis	0.383	0.292	0.340	0.162	0.198	0.274	0.105	0.000	0.236	0.151	0.346	0.379	0.305
H_rudolfensis	0.369	0.261	0.506	0.120	0.130	0.130	0.167	0.236	0.000	0.196	0.447	0.386	0.419
H_sapiens	0.434	0.333	0.468	0.092	0.098	0.207	0.054	0.151	0.196	0.000	0.390	0.417	0.410
P_boisei	0.552	0.431	0.456	0.391	0.431	0.431	0.291	0.346	0.447	0.390	0.000	0.151	0.515
P_robustus	0.553	0.412	0.619	0.373	0.370	0.361	0.412	0.379	0.386	0.417	0.151	0.000	0.644
Pan_troglodytes	0.304	0.296	0.155	0.303	0.352	0.352	0.243	0.305	0.419	0.410	0.515	0.644	0.000'
textcon = textConnection(bdist_pl_output_txt)
bdists_bdist_pl = read.table(textcon)
close(textcon)
bdists_bdist_pl

bdists_minbaru - bdists_bdist_pl

# Are all of the differences below 0.025?
differences_TF = sum(abs(bdists_minbaru - bdists_bdist_pl) < 0.025)

# Run the test
expect_true(object=differences_TF)
cat("...DONE\n")


# (These minor differences are probably due to 
#  auto-recoding of the dataset by read_nexus_data2. The only 
#  one that is concerning is 
#  H_sapiens, row 83 subset, row 124 full
#  "01"  "012" --> Fixed in bdist input dataset
# 
head(bdist_chardf_orig)
head(chardf_orig)
TFx = bdist_chardf_orig == chardf_orig
sum(TFx == FALSE)

TF1 = bdist_chardf_orig == "01"
TF2 = chardf_orig == "012"
TF = (TF1 + TF2) == 2
sum(TF)

rownums_table = matrix(data=1:nrow(chardf_orig), nrow=nrow(chardf_orig), ncol=ncol(chardf_orig), byrow=FALSE)
colnums_table = matrix(data=1:ncol(chardf_orig), nrow=nrow(chardf_orig), ncol=ncol(chardf_orig), byrow=TRUE)
rownums_table[TF]
# 124 (83 in subset)
colnums_table[TF]
#10
bdist_chardf_orig[124,]
chardf_orig[124,]

#








cat("=========================ENDING TESTS=========================\n\n")
