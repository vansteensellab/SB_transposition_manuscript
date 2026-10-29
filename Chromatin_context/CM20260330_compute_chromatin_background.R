#This is a script that takes the chromatin context of a large number of random genomic bins, filters out problematic bins (X and Y chromosome, encode blacklist regions, low (or too high) mappability regions) and computes mean and standard deviation.
#Conceptually based on the analysis in /home/c.moene/projects/hopping_technology/01_genomewide_chromatin_background/CM20260226_explore_background_location_set.Rmd

# input:
## sum_txt = output file from the chromatin context pipeline
## blacklist_bed = encode blacklist regions, tibble starting with columns "chr", "start", "end" 
## log_features = vector with the column names of columns containing features that should not be logtransformed
## non_cistrome_features = samples not processed by the cistrome pipeline (and not in 'log features'), for which the same psuedocounting should not be used
## input_cols = column names of the columns containing the processed ChIP input samples


chrom_background = function(sum_txt, blacklist_bed, 
                            log_features = c("LB1_pADam", "H3K9me3_pADam"), 
                            non_cistrome_features = c("TTseq_fw", "TTseq_rv"),
                            input_cols = c("input_1", "input_2")){

  #clean up the column names of the pipeline output file
  colnames_bw_cov = gsub("'", "",  colnames(sum_txt))
  colnames_bw_cov2 = gsub("-", "_",  colnames_bw_cov)
  colnames_bw_cov2[colnames_bw_cov2=="#chr"] <- "chr"
  
  context_tib = sum_txt
  colnames(context_tib) = colnames_bw_cov2
  
  if(all(c("TTseq_fw", "TTseq_rv") %in% colnames(context_tib))){
    context_tib = context_tib |>
      dplyr::mutate(TTseq_total = TTseq_fw + TTseq_rv)
  }
  ############
  #transformations
  ###########
  #pivot longer and annotate special features
  context_tib_long = tidyr::pivot_longer(context_tib, 
                                         cols = -c("chr", "start", "end"), 
                                         names_to = "feature", values_to = "signal") |>
    
    dplyr::mutate(signal_type = dplyr::case_when(feature %in% log_features ~ "logspace",
                                   feature %in% c(non_cistrome_features, "TTseq_total") ~ 'raw',
                                   .default = 'raw_cistrome'))
  
  #compute and add psuedocounts
  min_nonzero = context_tib_long |>
    dplyr::filter(signal > 0 & !is.na(signal) & !is.infinite(signal)) |>
    dplyr::group_by(feature) |>
    dplyr::summarize(minsignal = min(signal)) |>
    dplyr::mutate(psuedocount = 73*minsignal)
  
  
  # #add psuedocounts where necessary, then logtransform
  context_tib_ps = context_tib_long |>
    #we can remove NA and inf values anyways
    dplyr::filter( !is.na(signal) & !is.infinite(signal)) |>
    dplyr::left_join(min_nonzero) |>
    dplyr::mutate(signal_ps = signal + psuedocount,
           signal_log2 = dplyr::case_when(signal_type== "raw_cistrome" ~ log2(signal_ps),
                                   signal_type == "raw" & signal == 0 ~ NA,
                                   signal_type == "raw" & signal >0 ~ log2(signal),
                                   signal_type == "logspace" ~ signal)) #was already log transformed before the zscoring
  
  rm(context_tib_long)
  
  ##############
  #filtering
  ################
  
  #compute cut-offs for input samples
  sum_stats_input = context_tib_ps |>
    dplyr::filter(feature %in% input_cols) |>
    dplyr::filter(signal > 0 & !chr %in% c("chrX", "chrY")) |>
    dplyr::group_by(feature) |>
    dplyr::summarize(mean_input = mean(signal_log2),
                     sd_input = sd(signal_log2)) |>
    dplyr::mutate(cutoff_min = mean_input - 2*sd_input,
                  cutoff_max = mean_input + 2*sd_input)
  
  #mark bins with too high or low input signal
  bins_to_remove_lowmap = context_tib_ps |>
    dplyr::filter(feature %in% input_cols ) |>
    dplyr::left_join(sum_stats_input) |>
    dplyr::filter(signal_log2 < cutoff_min | signal_log2 > cutoff_max) |>
    dplyr::mutate(bin_remove_lowmap = 'remove') |>
    dplyr::select(chr, start, end, bin_remove_lowmap)
  
  
  #overlap blacklist region
  blacklist_gr = GenomicRanges::GRanges(blacklist_bed)
  random_bins_gr = GenomicRanges::GRanges(context_tib[,c("chr", "start", "end")])
  
  bins_in_blacklist = GenomicRanges::findOverlaps(random_bins_gr,blacklist_gr, ignore.strand = T )
  
  #find all bins to remove
  bins_to_remove = context_tib[,c("chr", "start", "end")] |> #original wide tibble, one row per bin
    dplyr::mutate(rownr = dplyr::row_number()) |>
    #mark bins overlapping blacklist or XY chrom
    dplyr::mutate(bin_remove_blacklist = dplyr::case_when(rownr %in% S4Vectors::queryHits(bins_in_blacklist) ~ "remove",
                                                          chr %in% c("chrX", "chrY") ~ "remove",
                                                          .default = 'keep')) |>
    #mark bins with mappability issue
    dplyr::left_join(bins_to_remove_lowmap) |>
    dplyr::filter(bin_remove_lowmap == "remove" | bin_remove_blacklist == "remove") |>
    dplyr::mutate(removal = "remove") |>
    #collect all bins to remove, deduplicate
    dplyr::select(chr, start, end, removal) |>
    dplyr::distinct()
     
  #remove the bins 
  context_ps_filt = context_tib_ps |>
    dplyr::left_join(bins_to_remove) |>
    dplyr::filter(is.na(removal)) #if the bin should be removed, it contains 'remove', so keep all NAs
  
 
  #calculate summary statistics
  sum_values_log2 = context_ps_filt |>
    dplyr::group_by(feature) |>
    dplyr::filter(!is.na(signal_log2))|> #remove the NA log2 values (are 0 signal values) from TTseq, so we can also calculate there
    dplyr::summarize(mean_log2 = mean(signal_log2),
              sd_log2 = sd(signal_log2),
              N_bins = dplyr::n()) |>
    #add the psuedocounting data (to process locations of interest later)
    dplyr::left_join(min_nonzero) |>
    #add input cutoffs (to process locations of interest later)
    dplyr::left_join(sum_stats_input) |>
    dplyr::select(-c('mean_input', 'sd_input', 'minsignal'))
  
  sum_values_log2
  
}

# Now write a script that takes the chromatin context at a specific set of locations and performs the normalizations on this data


# input:
## sum_txt = output file from the chromatin context pipeline
## background_sum_tib = tibble, output from the function above
## blacklist_bed = encode blacklist regions, tibble starting with columns "chr", "start", "end" 
## log_features = vector with the column names of columns containing features that should not be logtransformed
## non_cistrome_features = samples not processed by the cistrome pipeline (and not in 'log features'), for which the same psuedocounting should not be used
## input_cols = column names of the columns containing the processed ChIP input samples


chrom_insert_sites = function(sum_txt, background_sum_tib, blacklist_bed, 
                            log_features = c("LB1_pADam", "H3K9me3_pADam"), 
                            non_cistrome_features = c("TTseq_fw", "TTseq_rv"),
                            input_cols = c("input_1", "input_2")){
  
  #clean up the column names of the pipeline output file
  colnames_bw_cov = gsub("'", "",  colnames(sum_txt))
  colnames_bw_cov2 = gsub("-", "_",  colnames_bw_cov)
  colnames_bw_cov2[colnames_bw_cov2=="#chr"] <- "chr"
  
  context_tib = sum_txt
  colnames(context_tib) = colnames_bw_cov2
  
  if(all(c("TTseq_fw", "TTseq_rv") %in% colnames(context_tib))){
    context_tib = context_tib |>
      dplyr::mutate(TTseq_total = TTseq_fw + TTseq_rv)
  }
  
  ############
  #transformations
  ###########
  #pivot longer and annotate special features
  context_tib_long = tidyr::pivot_longer(context_tib, 
                                         cols = -c("chr", "start", "end"), 
                                         names_to = "feature", values_to = "signal") |>
    
    dplyr::mutate(signal_type = dplyr::case_when(feature %in% log_features ~ "logspace",
                                                 feature %in% c(non_cistrome_features, "TTseq_total") ~ 'raw',
                                                 .default = 'raw_cistrome'))
  
  # add psuedocount and log transform
  context_tib_transf = context_tib_long|>
    # remove NA and inf values 
    dplyr::filter( !is.na(signal) & !is.infinite(signal)) %>%
    dplyr::left_join(background_sum_tib) %>%
    dplyr::mutate(
      #psuedocount
      signal_ps = signal + psuedocount,
      #log transform when needed
      signal_log2 = dplyr::case_when(signal_type== "raw_cistrome" ~ log2(signal_ps),
                                     signal_type == "raw" & signal == 0 ~ NA,
                                     signal_type == "raw" & signal >0 ~ log2(signal),
                                     signal_type == "logspace" ~ signal), #was already log transformed before the zscoring
      # z-score transformation
      zscore = (signal_log2 - mean_log2)/sd_log2)
  
  ###############
  #mark problematic bins (dont remove them)
  ###########
  #mark bins with too high or low input signal
  bins_to_remove_lowmap = context_tib_transf |>
    dplyr::filter(feature %in% input_cols ) |>
    dplyr::filter(signal_log2 < cutoff_min | signal_log2 > cutoff_max) |>
    dplyr::mutate(bin_remove_lowmap = 'remove') |>
    dplyr::select(chr, start, end, bin_remove_lowmap)
  
  #overlap blacklist region
  blacklist_gr = GenomicRanges::GRanges(blacklist_bed)
  random_bins_gr = GenomicRanges::GRanges(context_tib[,c("chr", "start", "end")])
  
  bins_in_blacklist = GenomicRanges::findOverlaps(random_bins_gr,blacklist_gr, ignore.strand = T )
  
  bins_to_remove = context_tib[,c("chr", "start", "end")] |> #original wide tibble, one row per bin
    dplyr::mutate(rownr = dplyr::row_number()) |>
    #mark bins overlapping blacklist or XY chrom
    dplyr::mutate(bin_remove_blacklist = dplyr::case_when(rownr %in% S4Vectors::queryHits(bins_in_blacklist) ~ "remove",
                                                          chr %in% c("chrX", "chrY") ~ "remove",
                                                          .default = 'keep')) |>
    #mark bins with mappability issue
    dplyr::left_join(bins_to_remove_lowmap) |>
    dplyr::filter(bin_remove_lowmap == "remove" | bin_remove_blacklist == "remove") |>
    dplyr::mutate(removal = "remove") |>
    #collect all bins to remove, deduplicate
    dplyr::select(chr, start, end, removal) |>
    dplyr::distinct()
  
  context_tib_transf |>
    dplyr::left_join(bins_to_remove) |>
    tidyr::replace_na(list(removal = 'keep'))
}

