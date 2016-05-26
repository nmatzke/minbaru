
# Fast version of bdist
# This assumes literal matching of either strings, or e.g. 0 matching 01, with grepl()
# It also assumes you have already removed characters and taxa below the desired
# thresholds
fast_bdist <- function(chardf, alphabetize=TRUE, grepl_w_brackets=TRUE, remove_brackets_from_chardf=TRUE, num_bootstraps=0, is_this_a_bootstrap_subrun=FALSE, printflag=TRUE)
	{
	defaults='
	# Read in an example file, leaving it as original 
	# from the example (e.g., "10" different from "01").
	# HOWEVER, 12, 13, 14 is NOT a valid character pattern, it has to be converted to e.g. B, C, D
	# (allowed states are 0-9,A-...)
	fn = "/drives/GDrive/__github/minbaru/inst/extdata/examples/ex1_from_bdist_online.txt"
	fn = "/drives/GDrive/z_IDnew4/Wood_baraminology/2006_BDIST/testfile_alphabetical.txt"
	readres = read_bdist_to_chardf(fn=fn, correct_numeric_above_10=TRUE, correct_w_get_numstates_per_char=FALSE, return_list=TRUE)
	chardf = readres$chardf[1:6,]
	alphabetize=TRUE
	grepl_w_brackets=TRUE
	remove_brackets_from_chardf=TRUE
	num_bootstraps=0
	is_this_a_bootstrap_subrun=FALSE
	'
	# If we are adding [] to multistate characters, we need to 
	# remove {} () for sure
	if (grepl_w_brackets == TRUE)
		{
		remove_brackets_from_chardf=TRUE
		} # END if (grepl_w_brackets == TRUE)

	if (remove_brackets_from_chardf == TRUE)
		{
		chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern="\\{", replacement="")
		chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern="\\}", replacement="")
		chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern="\\(", replacement="")
		chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern="\\)", replacement="")
		} # END if (remove_brackets_from_chardf == TRUE)
	
	# If you are keeping the raw {} brackets in chardf,
	# replace them with \\{ and \\}
	# (which translate as \\\\{ and \\\\} )
	if (remove_brackets_from_chardf == FALSE)
		{
		chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern="\\{", replacement="\\\\{")
		chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern="\\}", replacement="\\\\}")
		} # END if (grepl_w_brackets == TRUE)

	
	# Alphabetize if desired
	if (alphabetize == TRUE)
		{
		new_order = order(colnames(chardf))
		chardf = chardf[,new_order]
		} # END if (alphabetize == TRUE)
	species_names = names(chardf)
	
	# Rows are characters
	# Taxa are columns
	# (a classic chardf, i.e. character data.frame)
	numtaxa = ncol(chardf)
	numchars = nrow(chardf)
	
	# NumQs in each taxon
	TF_Qs1 = chardf == "?"
	TF_Qs2 = chardf == "-"
	TF_Qs = (TF_Qs1 + TF_Qs2) > 0
	numQs_by_taxon = colSums(TF_Qs)
	numchars_by_taxon = numchars - numQs_by_taxon
	
	# Vectorize the calculation across all pairs
	# http://stackoverflow.com/questions/6269526/r-applying-a-function-to-all-row-pairs-of-a-matrix-without-for-loop
	v <- 1:numtaxa
	allPairs <- combn(v, 2) # choose a pair from 1:numtaxa
	allPairs
	
	# For each pairwise comparison, run this function
	# to find which character pairs are invalid
	# due to one or both of them being "?" or "-"
	count_pairwise_invalids <- function(i, j, chardf)
		{
		# Count "?"
		TFq1 = chardf[,i] == "?"
		TFq2 = chardf[,j] == "?"
		TFq3 = chardf[,i] == "-"
		TFq4 = chardf[,j] == "-"
		TF = (TFq1 + TFq2 + TFq3 + TFq4) > 0
		count_of_invalid_comparisons = sum(TF)
		
		return(count_of_invalid_comparisons)
		}
	# Go across columns of allPairs with .margins=2
	counts_Qpairs = mapply(FUN=count_pairwise_invalids, i=allPairs[1,], j=allPairs[2,], MoreArgs=list(chardf=chardf))
	
	# Cut down the pairs to exclude invalid comparisons (all "?" to "?")
	allPairs2 = allPairs[, counts_Qpairs != numchars]
	
	grepl_w_length2_input <- function(inputs)	
		{
		TF = grepl(pattern=inputs[1], x=inputs[2])
		return(TF)
		}
	
	# For each pairwise comparison, run this function
	count_pairwise_matches_w_grepl <- function(index_pair, chardf, exclude="-")
		{
		i = index_pair[1]
		j = index_pair[2]
		
		# Count matches (with grepl)
		two_columns_values1 = cbind(chardf[,i], chardf[,j])
		two_columns_values2 = cbind(chardf[,j], chardf[,i])
		
		TF1 = apply(X=two_columns_values1, MARGIN=1, FUN=grepl_w_length2_input)
		TF2 = apply(X=two_columns_values2, MARGIN=1, FUN=grepl_w_length2_input)
		
		# If the comparison included one or more "-", switch the 
		# potential match to "FALSE"
		excludeTF1 = chardf[,i] == exclude
		excludeTF2 = chardf[,j] == exclude
		TF = (TF1 + TF2) > 0
		TF[excludeTF1] = FALSE 
		TF[excludeTF2] = FALSE 
		count_of_matches = sum(TF)
		#cat("\n", i, ",", j, "	=	", count_of_matches, sep="")
		
		return(count_of_matches)
		}
	
	
	# Can't use actual "?" in chardf, screws up counts with grepl
	chardf_noQs = chardf
	chardf_noQs[chardf_noQs == "?"] = "-"
	
	# Add brackets {}, to make the grepl search flexible
	# (e.g., {012} matches 123)
	if (grepl_w_brackets == TRUE)
		{
		add_brackets_TF = sapply(X=chardf_noQs, FUN=nchar) > 1
		numvals = sum(add_brackets_TF)
		newvals_mat = cbind("[", chardf_noQs[add_brackets_TF], "]")
		newvals = apply(X=newvals_mat, MARGIN=1, FUN=paste0, collapse="")
		newvals
		chardf_noQs[add_brackets_TF] = newvals
		} # END if (grepl_w_brackets == TRUE)
	chardf_noQs
	#matches_count = apply(X=allPairs2, MARGIN=2, FUN=count_pairwise_matches_w_grepl, chardf=chardf)
	#matches_count
	
	# This excludes matching "-" and "0" (formerly ? and 0), but
	# e.g. "-" and "-" will still match
	matches_count = apply(X=allPairs2, MARGIN=2, FUN=count_pairwise_matches_w_grepl, chardf=chardf_noQs)
	matches_count
	
	# Find the 
	expanded_matches_count = rep(0, length(counts_Qpairs))
	expanded_matches_count[counts_Qpairs != numchars] = matches_count
	
	rbind(allPairs, counts_Qpairs, expanded_matches_count)
	
	count_Qpairs_mat = matrix(data=NA, ncol=numtaxa, nrow=numtaxa)
	count_matches_mat = matrix(data=NA, ncol=numtaxa, nrow=numtaxa)
	for (i in 1:ncol(allPairs))
		{
		count_Qpairs_mat[allPairs[1,i], allPairs[2,i]] = counts_Qpairs[i]
		count_Qpairs_mat[allPairs[2,i], allPairs[1,i]] = counts_Qpairs[i]
		count_matches_mat[allPairs[1,i], allPairs[2,i]] = expanded_matches_count[i]
		count_matches_mat[allPairs[2,i], allPairs[1,i]] = expanded_matches_count[i]	
		}
	diag(count_Qpairs_mat) = numQs_by_taxon
	diag(count_matches_mat) = numchars_by_taxon
	
	count_Qpairs_mat
	count_matches_mat
	
	numchars_mat = matrix(data=numchars, nrow=numtaxa, ncol=numtaxa)

	# The Coefficient of Baraminic Distance
	# dij = mij / nij
	# mij = "number of mismatched characters between the ith and jth organism"
	# nij = "number of compared characters" between ith and jth organism
	# 
	# (dij is the "complementary form of the simple matching coefficient introduced by 
	#  Sokal and Michener (1958))
	# 
	# dij, "baraminic distance represents an estimator of a true biological distance 
	# (Dij)."
	# 
	# The authors recommend "that all characters be given equal weight. In other words, 
	# all characters contribute 0 (match) or 1 (mismatch) to the baraminic distance 
	# numerator.
	#
	# Source: Equation 3, page 198 of:
	# Robinson, D. Ashley; Cavanaugh, David P. (1998). 
	#  A Quantitative Approach to Baraminology With Examples from the 
	#  Catarrhine Primates. Creation Research Society Quarterly, 34(4), 
	#  196-208.
	bdist_dij = 1 - (count_matches_mat / (numchars_mat - count_Qpairs_mat))
	bdist_dij
	count_Qpairs_mat = as.data.frame(count_Qpairs_mat, stringsAsFactors=FALSE)
	count_matches_mat = as.data.frame(count_matches_mat, stringsAsFactors=FALSE)
	numchars_mat = as.data.frame(numchars_mat, stringsAsFactors=FALSE)
	bdist_dij = as.data.frame(bdist_dij, stringsAsFactors=FALSE)
	
	names(count_Qpairs_mat) = species_names
	names(count_matches_mat) = species_names
	names(numchars_mat) = species_names
	names(bdist_dij) = species_names

	rownames(count_Qpairs_mat) = species_names
	rownames(count_matches_mat) = species_names
	rownames(numchars_mat) = species_names
	rownames(bdist_dij) = species_names



	# "Baraminic units can be separately elucidated with a correlation
	#  analysis of the baraminic distances of each pair of
	#  organisms. Since this analysis is not a matrix-level comparison
	#  (as is the case for the criterial correlation analysis), the
	#  significance of Pearson product-moment correlation coefficients
	#  can be estimated in a classical manner with n – 2 degrees
	#  of freedom, where n equals the number of species."
	#  
	# "Positive correlations indicate the compared organisms share
	#  a similar set of distances that differ from other organisms by
	#  a similar magnitude. This condition may therefore be diagnostic
	#  of a monobaraminic relationship. In contrast, negative
	#  correlations indicate the compared species contain
	#  antithetical patterns of baraminic distances. Taxa similar to
	#  species X would be dissimilar to species Y and vice versa,
	#  which suggests species X and Y may be classified as apobaraminic.
	#  The baraminic relationship between organisms
	#  with uncorrelated baraminic distances is considered unresolved,
	#  but a monobaraminic affinity or apobaraminic discontinuity
	#  can be cautiously hypothesized with well
	#  supported trends."
	# 
	#  pp. 200-201 of: Robinson, D. Ashley; Cavanaugh, David P. (1998). 
	#  A Quantitative Approach to Baraminology With Examples from the 
	#  Catarrhine Primates. Creation Research Society Quarterly, 34(4), 
	#  196-208.
	# 

	# Pearson's product-moment correlation,
	# run on a matrix of distances
	cormat = cor(x=bdist_dij, method="pearson")
	cormat = as.data.frame(cormat, stringsAsFactors=FALSE)
	names(cormat) = names(chardf)
	row.names(cormat) = names(chardf)

	# Get p-values 
	# corr.test uses cor to find the correlations for either complete or pairwise data 
	# and reports the sample sizes and probability values as well. For symmetric 
	# matrices, raw probabilites are reported below the diagonal and correlations 
	# adjusted for multiple comparisons above the diagonal.
	corrtest_res = psych::corr.test(x=bdist_dij, y=NULL, use="pairwise", method="pearson", adjust="holm", alpha=0.05)
	cor_Ps_mat = corrtest_res$p
	cor_Ps_mat = as.data.frame(cor_Ps_mat, stringsAsFactors=FALSE)
	names(cor_Ps_mat) = names(chardf)
	row.names(cor_Ps_mat) = names(chardf)
		
	# Replace upper triangle with the no-corrections-for-sample-size P-values
	# (in order to match bdist / bdistmds)
	bdist_cor_Ps = cor_Ps_mat
	bdist_cor_Ps[upper.tri(bdist_cor_Ps)] = t(bdist_cor_Ps)[upper.tri(bdist_cor_Ps)]
	bdist_cor_Ps = as.data.frame(bdist_cor_Ps, stringsAsFactors=FALSE)
	names(bdist_cor_Ps) = names(chardf)
	row.names(bdist_cor_Ps) = names(chardf)
	

	# Replace lower triangle with the corrections-for-sample-size P-values
	# (in order to have a sample-size corrected version of bdist/bdistmds)
	cor_Ps_multiple_test_corrected = cor_Ps_mat
	cor_Ps_multiple_test_corrected[lower.tri(cor_Ps_multiple_test_corrected)] = t(cor_Ps_multiple_test_corrected)[lower.tri(cor_Ps_multiple_test_corrected)]
	cor_Ps_multiple_test_corrected = as.data.frame(cor_Ps_multiple_test_corrected, stringsAsFactors=FALSE)
	names(cor_Ps_multiple_test_corrected) = names(chardf)
	row.names(cor_Ps_multiple_test_corrected) = names(chardf)


	res = NULL
	res$numchars_mat = numchars_mat
	res$count_Qpairs_mat = count_Qpairs_mat
	res$count_matches_mat = count_matches_mat
	res$bdist_dij = bdist_dij
	res$cormat = cormat
	res$corrtest_res = corrtest_res
	res$cor_Ps_mat = cor_Ps_mat
	res$bdist_cor_Ps = bdist_cor_Ps
	res$cor_Ps_multiple_test_corrected = cor_Ps_multiple_test_corrected

	# Bootstrapping
	#######################################################
	# Bootstrap the chardf
	#######################################################
	#num_bootstraps = 100
	if ((num_bootstraps > 0) && (is_this_a_bootstrap_subrun == FALSE))
		{
		num_chars = nrow(chardf)
		list_of_bootstrapped_datasets = list()
		list_of_bdist_res = list()
		list_of_pval_matrices = list()
		array_of_pval_matrices = array(data=NA, dim=c(ncol(chardf), ncol(chardf), num_bootstraps))
		
		if (printflag)
			{
			cat("\nMaking ", num_bootstraps, " bootstrap replicates:\n")
			}
		for (iiii in 1:num_bootstraps)
			{
			if (printflag)
				{
				cat(iiii, " ", sep="")
				}
			sampled_chars = sample(x=1:num_chars, size=num_chars, replace=TRUE, prob=NULL)
			list_of_bootstrapped_datasets[[iiii]] = chardf[sampled_chars,]
			# chardf = list_of_bootstrapped_datasets[[iiii]]
			list_of_bdist_res[[iiii]] = fast_bdist(chardf=list_of_bootstrapped_datasets[[iiii]], alphabetize=FALSE, grepl_w_brackets=grepl_w_brackets, remove_brackets_from_chardf=remove_brackets_from_chardf, num_bootstraps=0, is_this_a_bootstrap_subrun=TRUE, printflag=printflag)
			
			#calc_bdist_stats(list_of_bootstrapped_datasets[[iiii]], min_char_relevance=min_char_relevance, min_taxic_relevance=min_taxic_relevance, parse_ambig=parse_ambig, do_frequencies_add_to_1=do_frequencies_add_to_1, use_ttl_numsp_in_charDiversity_calcs=use_ttl_numsp_in_charDiversity_calcs, use_num_valid_comparisons_in_bdist_calcs=use_num_valid_comparisons_in_bdist_calcs, split_on=split_on, split_multichars=split_multichars, dont_compare_ambig=dont_compare_ambig, div_matches_by_numchars=div_matches_by_numchars, max1_match=max1_match, alphabetical=alphabetical, num_bootstraps=0, is_this_a_bootstrap_subrun=TRUE, keep_taxa_TF=keep_taxa_TF)
			#cft(list_of_bdist_res[[iiii]]$bdist_cor_Ps)
			#print(dim(list_of_bdist_res[[iiii]]$bdist_cor_Ps))
			#print(dim(array_of_pval_matrices))
			list_of_pval_matrices[[iiii]] = list_of_bdist_res[[iiii]]$bdist_cor_Ps
			array_of_pval_matrices[,,iiii] = as.matrix(list_of_bdist_res[[iiii]]$bdist_cor_Ps)
			} # END for (i in 1:num_bootstraps)


		TFs_array = array_of_pval_matrices < 0.05
		bootstrap_freqs_pLT05 = apply(X=TFs_array, MARGIN=c(1,2), FUN=sum)
		bootstrap_freqs_pLT05 = as.data.frame(bootstrap_freqs_pLT05, stringsAsFactors=FALSE)
		names(bootstrap_freqs_pLT05) = names(chardf)
		row.names(bootstrap_freqs_pLT05) = names(chardf)
		bootstrap_freqs_pLT05
		
		res$bootstrap_freqs_pLT05 = bootstrap_freqs_pLT05
		} # END if (num_bootstraps > 0)

	extract='
	numchars_mat = res$numchars_mat
	count_Qpairs_mat = res$count_Qpairs_mat
	count_matches_mat = res$count_matches_mat
	bdist_dij = res$bdist_dij
	cormat = res$cormat
	corrtest_res = res$corrtest_res
	cor_Ps_mat = res$cor_Ps_mat
	bdist_cor_Ps = res$bdist_cor_Ps
	cor_Ps_multiple_test_corrected = res$cor_Ps_multiple_test_corrected
	bootstrap_freqs_pLT05 = res$bootstrap_freqs_pLT05
	'
	
	return(res)
	}









# Baraminic distance plot (as in bdistmds). 
# 
# This is a plot *not* of correlations between the characters of taxa (absolute similarity
# in character states). Instead, it is a
# plot of the correlations of the between-taxon distances (baraminic distances -- 
# basically a phenetic similarity measure).  The correlations of the *distances*
# give more of a relative scaling of similarity within the dataset.  Creationists think
# they can use this to identify baramins (God's specially created kinds), but I think 
# this measure has got to be nothing more than a statement of relative similarity, and 
# will depend entirely on the taxa sampled. If your dataset is hominins, then 
# Australopithecus and Homo species may have distances to other species in the data
# that are uncorrelated or negatively correlated. But if your dataset is primates,
# Australopithecus and Homo species will have distances to other species that will be
# pretty darn correlated.  It's all quite silly.
# 
# Black closed squares: BDS positive correlation, p<0.05, bootstrap>90%
# Grey closed squares: BDS positive correlation, p<0.05, bootstrap<=90%
# Black open circles: BDS negative correlation, p<0.05, bootstrap>90%
# Grey open circles: BDS negative correlation, p<0.05, bootstrap<=90%
# 
# taxon_to_order_on Order the output plot by baraminic distance to this taxon
# new_order Gives a user-specified order for the taxa. Overrides taxon_to_order_on
# 
bds_bootstrap_plot <- function(res, taxon_to_order_on=NULL, new_order=NULL, leftlabels=TRUE, toplabels=FALSE)
	{
	defaults='
	taxon_to_order_on = "H_sapiens"
	new_order=NULL
	'
	# 
	if (!is.null(taxon_to_order_on) && !is.null(new_order))
		{
		#cat("\nbds_bootstrap_plot() says: Warning: new_order is not NULL, this overrides taxon_to_order_on\n")
		taxon_to_order_on = NULL
		}

	tmp_bdist_dij = res$bdist_dij
	tmp_cormat = res$cormat
	tmp_bdist_cor_Ps = res$bdist_cor_Ps
	tmp_bootstrap_freqs_pLT05 = res$bootstrap_freqs_pLT05
	numrows = nrow(tmp_cormat)
	nums = 1:numrows

	
	if ( (length(new_order)==1) && (new_order == "heatmap") )
		{
		print("heatmapping")
		heatmap_res = heatmap(t(tmp_cormat))
		title("Heatmap of BDC matrix")
		new_order = heatmap_res$rowInd
		}

	# Sort BDC (Baraminic Distance Correlations) by distance
	# to a taxon
	if (is.null(taxon_to_order_on) == FALSE)
		{
		bdc_to_taxon = tmp_bdist_dij[rownames(tmp_cormat)==taxon_to_order_on,]
		new_order2 = order(bdc_to_taxon)

		names_to_plot = rownames(tmp_cormat)[new_order2]
		ordered_cormat = tmp_cormat[,new_order2]
		ordered_cormat = ordered_cormat[new_order2,]

		ordered_bdist_cor_Ps = tmp_bdist_cor_Ps[,new_order2]
		ordered_bdist_cor_Ps = ordered_bdist_cor_Ps[new_order2,]
	
		ordered_bootstrap_freqs_pLT05 = tmp_bootstrap_freqs_pLT05[,new_order2]
		ordered_bootstrap_freqs_pLT05 = ordered_bootstrap_freqs_pLT05[new_order2,]
		} 
		

	if (is.null(new_order) == FALSE)
		{
		# Now, apply reversal or manual order
		names_to_plot = rownames(tmp_cormat)
		ordered_cormat = tmp_cormat
		ordered_bdist_cor_Ps = tmp_bdist_cor_Ps
		ordered_bootstrap_freqs_pLT05 = tmp_bootstrap_freqs_pLT05
		
		names_to_plot = names_to_plot[new_order]
		ordered_cormat = ordered_cormat[,new_order]
		ordered_cormat = ordered_cormat[new_order,]
		ordered_bdist_cor_Ps = ordered_bdist_cor_Ps[,new_order]
		ordered_bdist_cor_Ps = ordered_bdist_cor_Ps[new_order,]
		ordered_bootstrap_freqs_pLT05 = ordered_bootstrap_freqs_pLT05[,new_order]
		ordered_bootstrap_freqs_pLT05 = ordered_bootstrap_freqs_pLT05[new_order,]
		}	

	# Reverse
	new_order = rev(nums)
	names_to_plot = names_to_plot[new_order]
	ordered_cormat = ordered_cormat[,new_order]
	ordered_cormat = ordered_cormat[new_order,]
	ordered_bdist_cor_Ps = ordered_bdist_cor_Ps[,new_order]
	ordered_bdist_cor_Ps = ordered_bdist_cor_Ps[new_order,]
	ordered_bootstrap_freqs_pLT05 = ordered_bootstrap_freqs_pLT05[,new_order]
	ordered_bootstrap_freqs_pLT05 = ordered_bootstrap_freqs_pLT05[new_order,]

	# Manual order, or just reverse
	if (is.null(new_order) == TRUE)
		{
		new_order = rev(nums)
		} else {
		new_order = new_order 
		}



	# Black closed squares: BDS positive correlation, p<0.05, bootstrap>90%
	# Grey closed squares: BDS positive correlation, p<0.05, bootstrap<=90%
	# Black open circles: BDS negative correlation, p<0.05, bootstrap>90%
	# Grey open circles: BDS negative correlation, p<0.05, bootstrap<=90%


	# Black closed squares: BDS positive correlation, p<0.05, bootstrap>90%
	black_closed_squares_TF = ((ordered_cormat > 0) + (ordered_bdist_cor_Ps < 0.05) + (ordered_bootstrap_freqs_pLT05 > 90)) == 3
	black_closed_squares_TF
	grey_closed_squares_TF = ((ordered_cormat > 0) + (ordered_bdist_cor_Ps < 0.05) + (ordered_bootstrap_freqs_pLT05 <= 90)) == 3
	black_open_circles_TF = ((ordered_cormat < 0) + (ordered_bdist_cor_Ps < 0.05) + (ordered_bootstrap_freqs_pLT05 > 90)) == 3
	grey_open_circles_TF = ((ordered_cormat < 0) + (ordered_bdist_cor_Ps < 0.05) + (ordered_bootstrap_freqs_pLT05 <= 90)) == 3

	xnums = matrix(data=nums, nrow=numrows, ncol=numrows, byrow=TRUE)
	ynums = matrix(data=nums, nrow=numrows, ncol=numrows, byrow=FALSE)
	
	if (leftlabels == TRUE)
		{
		leftmargin = 10
		} else {
		leftmargin = 4
		}
	if (toplabels == TRUE)
		{
		topmargin = 10
		} else {
		topmargin = 4
		}
	
	par(mar=c(4,leftmargin,topmargin,1))
	plot(x=nums, y=nums, pch=".", col="white",  xaxt="n", yaxt="n", xlab="", ylab="")
	if (leftlabels == TRUE)
		{
		mtext(text=names_to_plot, side=2, line=0.5, at=nums, padj=0, las=1)
		}
	if (toplabels == TRUE)
		{
		mtext(text=names_to_plot, side=3, line=0.5, at=nums, padj=0, las=3)
		}

	points(x=xnums[black_closed_squares_TF], y=ynums[black_closed_squares_TF], pch=15, col="black", cex=2)
	points(x=xnums[grey_closed_squares_TF], y=ynums[grey_closed_squares_TF], pch=15, col="grey", cex=2)
	points(x=xnums[black_open_circles_TF], y=ynums[black_open_circles_TF], pch=1, col="black", cex=2)
	points(x=xnums[grey_open_circles_TF], y=ynums[grey_open_circles_TF], pch=1, col="grey", cex=2)
	
	} # END bds_bootstrap_plot <- function(res, taxon_to_order_on=NULL, new_order=NULL)

