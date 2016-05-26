# Read a bdist input file, like this:

# freddie 0 {10} 1 12    0 1 {01}
# velma   0 {01} 2 13 {12} 0    ?
# scooby  1    0 2 12    1 1    1
# shaggy  0    1 0 14    2 0 {01}
# daphne  0    0 2 12    0 1 {23}


read_bdist_to_chardf <- function(fn, min_char_relevance=0, min_taxic_relevance=0, correct_numeric_above_10=TRUE, correct_w_get_numstates_per_char=TRUE, return_list=TRUE, alphabetize=TRUE, printflag=TRUE, subset_data=TRUE)
	{
	defaults='
	fn = "/drives/GDrive/z_IDnew4/Wood_baraminology/2006_BDIST/testfile.txt"
	correct_numeric_above_10 = TRUE
	correct_w_get_numstates_per_char=TRUE
	'
	
	
	
	# Get the number of lines
	tmplines = scan(fn, what="character", skip=0, strip.white=TRUE, blank.lines.skip=TRUE, sep="\n")
	numlines = length(tmplines)
	
	# Read in to matrix
	tmplines = scan(fn, what="character", skip=0, strip.white=TRUE, blank.lines.skip=TRUE)
	
	# Convert ambiguities to standard format
	if ( any(grepl(pattern="\\{", x=tmplines)) == TRUE)
		{
		tmplines = gsub(pattern="\\{", replacement="(", x=tmplines)
		tmplines = gsub(pattern="\\}", replacement=")", x=tmplines)
		}
	tmplines
	
	# Rows are species, columns are characters
	tmpmat = matrix(data=tmplines, nrow=numlines, byrow=TRUE)
	tmpmat


	# If you have values above 10 as inputs (no ambiguities), convert to A, B, C, etc.
 	if (correct_numeric_above_10 == TRUE)
 		{
		# TNT can do 32 character states
		alphabet = LETTERS
		statecodes = c(0:9, LETTERS[1:22])
		
		# Find the numeric places:
		# Rows are species, columns are characters
		tmpmat2 = suppressWarnings(matrix(data=as.numeric(tmpmat[,-1]), nrow=numlines, byrow=FALSE))
		tmpmat2
		
		# Find the numerics (non-NAs) that are 10 or above
		TF = matrix(data=tmpmat2>9, nrow=numlines, byrow=FALSE)
		TF[is.na(TF)] = FALSE
		TF
		
		# Replace with statecodes
		tmpmat2[TF] = statecodes[tmpmat2[TF]]
		tmpmat2
		
		# Insert back in:
		tmpmat[,-1][TF] = tmpmat2[TF]
 		} # END if (correct_numeric_above_10 = TRUE)

	# Convert to data.frame
	species_names = tmpmat[,1]
	chardf = as.data.frame(x=t(tmpmat[,-1]), stringsAsFactors=FALSE)
	names(chardf) = species_names
	chardf
	
	# Save original chardf
	orig_chardf = chardf
	
	# Calculate character relevances, AND SUBSET
	# CHARACTER BUT NOT TAXA
	relevance_res = bdist_relevance(chardf=chardf,  min_char_relevance=min_char_relevance, min_taxic_relevance=min_taxic_relevance, alphabetize=alphabetize, subset_data=subset_data)
	# Extract
	chardf = relevance_res$chardf
	species_names = relevance_res$species_names
	taxic_relevance = relevance_res$taxic_relevance
	stats_table = relevance_res$stats_table
	
	if (printflag)
		{
		cat("\nread_bdist_to_chardf() says: Printing so-called 'taxic relevance' before and after cutting less complete characters. See also readres$taxic_relevance.")
		cat("\n\n")
		print(cft(taxic_relevance, numdigits_inbetween_have_fixed_digits=5))
		}

	# Correct, if desired
	if (correct_w_get_numstates_per_char == TRUE)
		{
		res = get_numstates_per_char(nexd=chardf, ambig_to_remove=c("\\(", "\\)", " ", ","), return_missing_chars="correct", printall="short", count_autapomorphies=TRUE, sort_ambigs=TRUE)
		missing_charstates_list = res$missing_charstates_list
		allQs_TF = res$allQs_TF
		chardf_corrected = res$nexdf
		main_chardf = chardf_corrected

		if (return_list == TRUE)
			{
			res$raw_chardf = chardf
			res$chardf_corrected = chardf_corrected
			res$chardf = main_chardf
			res$stats_table = stats_table
			res$taxic_relevance = taxic_relevance
			return(res)
			} else {
			return(main_chardf)
			} # END if (return_list == TRUE)
		} else {
		res = get_numstates_per_char(nexd=chardf, return_missing_chars="list", printall="short", count_autapomorphies=TRUE)
		missing_charstates_list = res$missing_charstates_list
		allQs_TF = res$allQs_TF

		if (return_list == TRUE)
			{
			res$raw_chardf = chardf
			res$chardf = chardf
			res$stats_table = stats_table
			res$taxic_relevance = taxic_relevance
			return(res)
			} else {
			return(main_chardf)
			} # END if (return_list == TRUE)
		} # END if (correct_w_get_numstates_per_char = TRUE)
	
	return(stop("STOP ERROR: read_bdist_to_chardf() says: this function should have exited earlier."))
	} # END read_bdist <- function()
