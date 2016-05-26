
bdist_subset <- function(chardf, relevance_res)
	{
	# Subset by 
	}


# Fix a NEXUS read_nexus_data2 chardf output for use in bdist
fix_NEXUS_charsdf_for_bdist <- function(chardf, remove_spaces=TRUE, alphabetize=TRUE)
	{
	if (alphabetize == TRUE)
		{
		chardf = chardf[, order(colnames(chardf))]
		}
	
	if (remove_spaces == TRUE)
		{
		chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern=" ", replacement="")
		}
	chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern="\\(", replacement="")
	chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern="\\)", replacement="")
	chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern="\\{", replacement="")
	chardf[!is.na(chardf)] = sapply(X=chardf, FUN=gsub, pattern="\\}", replacement="")
	return(chardf)
	}

# Fix a NEXUS read_nexus_data2 chardf output for use in bdist
fix_NEXUS_charslist_for_bdist <- function(chars_list, remove_spaces=TRUE, alphabetize=TRUE)
	{
	chardf = as.data.frame(x=chars_list, stringsAsFactors=FALSE)

	if (alphabetize == TRUE)
		{
		chardf = chardf[, order(colnames(chardf))]
		}

	chardf = fix_NEXUS_charsdf_for_bdist(chardf, remove_spaces=remove_spaces)
	return(chardf)
	}



# NOTE: I disagree with the whole idea of "relevance" -- everything is relevant!
# But, to replicate bdist / bdistmds, we have to reproduce these parts
# of the calculation
bdist_relevance <- function(chardf, min_char_relevance=0, min_taxic_relevance=0, alphabetize=TRUE, subset_data=TRUE)
	{
	defaults='
	min_char_relevance=0.75
	min_taxic_relevance=0
	alphabetize=TRUE
	subset_data=TRUE
	'
	
	# If alphabetizing:
	if (alphabetize == TRUE)
		{
		species_names = names(chardf)
		new_order = order(species_names)
		chardf = chardf[,new_order]
		} # END if (alphabetize == TRUE)
	# Be sure to re-initialize species names:
	species_names = names(chardf)


	# Taxic relevance before character cutoffs
	TF1 = chardf == "?"
	TF2 = chardf == "-"
	TF = TF1 + TF2
	taxic_relevance_rawdata = 1 - (colSums(TF) / (rep(nrow(chardf), times=ncol(chardf))))



	# Taxic relevance after character cutoffs
	TF1 = chardf == "?"
	TF2 = chardf == "-"
	TF = TF1 + TF2
	numQs_per_species = colSums(TF)
	numdata_per_species = nrow(chardf) - numQs_per_species
	ttl_nchar = numdata_per_species + numQs_per_species
	
	# Don't cut -- the user should make that decision deliberately
	#keep_taxa_TF = taxic_relevance > min_taxic_relevance
	#chardf = chardf[, keep_taxa_TF]
	
	# Make a table of taxic_relevance
	min_taxic_relevance_list = rep(min_taxic_relevance, times=length(species_names))
	taxic_rel_raw_above_cutoff_TF = taxic_relevance_rawdata > min_taxic_relevance_list
	taxic_relevance = cbind(1:length(species_names), species_names, numdata_per_species, numQs_per_species, ttl_nchar, min_taxic_relevance_list, taxic_relevance_rawdata, taxic_rel_raw_above_cutoff_TF)
	taxic_relevance = as.data.frame(taxic_relevance, stringsAsFactors=FALSE)
	names(taxic_relevance) = c("num", "name", "nchars", "nQs", "ttl_nchar", "minTaxicRel", "TaxicRel", "relTF")
	taxic_relevance$num = as.numeric(taxic_relevance$num)
	taxic_relevance$nchars = as.numeric(taxic_relevance$nchars)
	taxic_relevance$nQs = as.numeric(taxic_relevance$nQs)
	taxic_relevance$ttl_nchar = as.numeric(taxic_relevance$ttl_nchar)
	taxic_relevance$minTaxicRel = as.numeric(taxic_relevance$minTaxicRel)
	taxic_relevance$TaxicRel = as.numeric(taxic_relevance$TaxicRel)
	taxic_relevance$relTF = as.logical(taxic_relevance$relTF)



	# Make a table of character completeness
	#######################################################
	# Character relevance -- only "?" or "-" are irrelevant
	# (I guess)
	#######################################################
	TF1 = chardf == "?"
	TF2 = chardf == "-"
	TF = TF1 + TF2
	relevance = 1 - (rowSums(TF) / (rep(ncol(chardf), times=nrow(chardf))))
	relevance

	#min_char_relevance = 0.95
	include_TF = relevance >= min_char_relevance
	#chardf = chardf[include_TF,]
	
# 	if ( (sum(include_TF=FALSE)>0) && printflag )
# 		{
# 		cat("\read_bdist_to_chardf() says: ", sum(include_TF), "/", length(include_TF), " characters have completeness less than min_char_relevance=", min_char_relevance, " (a<", min_char_relevance, "). You should cut these (if you want to replicate bdist).")
# 		cat("\n")
# 		}

	# Number of Qs by character
	if (sum(include_TF) > 0)
		{
		TF1 = rowSums(chardf == "?")
		TF2 = rowSums(chardf == "-")
		} else {
		TF1 = rep(0, times=nrow(chardf))
		TF2 = rep(0, times=nrow(chardf))
		} # END if (sum(include_TF) > 0)

	numQs_per_char = TF1 + TF2

	num_species = rep(ncol(chardf), times=nrow(chardf))
	num_nonQ = ncol(chardf) - numQs_per_char
	#num_uniq_charstates
	#character_diversities


	# Repeat completeness calculations after taxa are cut
	keep_taxa_TF = taxic_relevance_rawdata > min_taxic_relevance
	chardf2 = chardf[, keep_taxa_TF]

	TF1 = chardf2 == "?"
	TF2 = chardf2 == "-"
	TF = TF1 + TF2
	relevance_sub = 1 - (rowSums(TF) / (rep(ncol(chardf2), times=nrow(chardf2))))

	include_TF_sub = relevance_sub >= min_char_relevance

	# Number of Qs by character
	if (sum(include_TF_sub) > 0)
		{
		TF1 = rowSums(chardf2 == "?")
		TF2 = rowSums(chardf2 == "-")
		} else {
		TF1 = rep(0, times=nrow(chardf2))
		TF2 = rep(0, times=nrow(chardf2))
		} # END if (sum(include_TF) > 0)

	numQs_per_char_sub = TF1 + TF2

	num_species_sub = rep(ncol(chardf2), times=nrow(chardf2))
	num_nonQ_sub = ncol(chardf2) - numQs_per_char_sub



	# Make an output table (BEFORE cutting taxa)
	charnum = 1:nrow(chardf)
	fract_nonQ = num_nonQ/num_species
	cutoffs = rep(min_char_relevance, times=nrow(chardf))
	stats_table = cbind(charnum, num_species, num_nonQ, fract_nonQ, relevance, cutoffs, include_TF)#, num_uniq_charstates, character_diversities)
	stats_table = as.data.frame(stats_table, stringsAsFactors=FALSE)
	names(stats_table) = c("charnum", "num_species", "num_nonQ", "fract_nonQ", "relevance", "cutoffs", "include_TF")#, "numUniqStates", "charDiversity")
	
	# Force the correct classes 
	stats_table$charnum = as.numeric(stats_table$charnum)
	stats_table$num_species = as.numeric(stats_table$num_species)
	stats_table$num_nonQ = as.numeric(stats_table$num_nonQ)
	stats_table$fract_nonQ = as.numeric(stats_table$fract_nonQ)
	stats_table$relevance = as.numeric(stats_table$relevance)
	stats_table$cutoffs = as.numeric(stats_table$cutoffs)
	stats_table$include_TF = as.logical(stats_table$include_TF)
	stats_table


	# Make an output table (AFTER cutting taxa)
	charnum_sub = 1:nrow(chardf2)
	num_species_sub = ncol(chardf2)
	fract_nonQ_sub = num_nonQ_sub/num_species_sub
	cutoffs_sub = rep(min_char_relevance, times=nrow(chardf2))
	stats_table_sub = cbind(charnum_sub, num_species_sub, num_nonQ_sub, fract_nonQ_sub, relevance_sub, cutoffs, include_TF_sub)#, num_uniq_charstates, character_diversities)
	stats_table_sub = as.data.frame(stats_table_sub, stringsAsFactors=FALSE)
	names(stats_table_sub) = c("charnum", "num_species", "num_nonQ", "fract_nonQ", "relevance", "cutoffs", "include_TF")#, "numUniqStates", "charDiversity")
	
	# Force the correct classes 
	stats_table_sub$charnum = as.numeric(stats_table_sub$charnum)
	stats_table_sub$num_species = as.numeric(stats_table_sub$num_species)
	stats_table_sub$num_nonQ = as.numeric(stats_table_sub$num_nonQ)
	stats_table_sub$fract_nonQ = as.numeric(stats_table_sub$fract_nonQ)
	stats_table_sub$relevance = as.numeric(stats_table_sub$relevance)
	stats_table_sub$cutoffs = as.numeric(stats_table_sub$cutoffs)
	stats_table_sub$include_TF = as.logical(stats_table_sub$include_TF)
	stats_table_sub
	
	# Save outputs
	relevance_res = NULL
	relevance_res$chardf = chardf
	
	# Subset the data?
	if (subset_data == TRUE)
		{
		subtaxa_chardf = chardf2
		subboth_chardf = chardf2[stats_table_sub$include_TF,]
		relevance_res$subtaxa_chardf = subtaxa_chardf
		relevance_res$subboth_chardf = subboth_chardf
		}
	
	
	relevance_res$species_names = species_names
	relevance_res$taxic_relevance = taxic_relevance
	relevance_res$stats_table = stats_table
	relevance_res$stats_table_sub = stats_table_sub
	
	# Extract
	extract='
	chardf = relevance_res$chardf
	species_names = relevance_res$species_names
	taxic_relevance = relevance_res$taxic_relevance
	stats_table = relevance_res$stats_table
	stats_table_sub = relevance_res$stats_table_sub

	subtaxa_chardf = relevance_res$subtaxa_chardf
	subboth_chardf = relevance_res$subboth_chardf
	'
	
	return(relevance_res)	
	} # END bdist_relevance <- function(chardf,  min_char_relevance=0, min_taxic_relevance=0)