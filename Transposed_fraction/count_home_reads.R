#load packages silently
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))

#collect snakemake variables
home_seq <- snakemake@params[["home_seq"]]
fq_file <- snakemake@params[["fastq_path"]]
sample_name <- snakemake@params[["sample_name"]]
LP <- snakemake@params[["lp"]]
orientation <- snakemake@params[["orientation"]]

output_file <- snakemake@output[[1]]

#define the home location sequence (16bp, because the parsed files contain sequences of min 16bp)
home_location_seq = DNAString(home_seq)

#load the parsed fastq file
parsed_fq <- readFastq(fq_file)  # ShortReadQ object

# Get first 16bp from each read
start_of_read <- subseq(sread(parsed_fq), 1, 16)

# find the reads that match the home location (max 2 mismatches, is 10% error rate rounded up)
# home_location_seq is 16bp long and start_of_read also, so if identical there will be exactly 1 match
read_match_home <- vcountPattern(home_location_seq, start_of_read, max.mismatch=2, with.indels=FALSE)

# count the matches
N_home_reads <- sum(read_match_home == 1)
N_total_reads <- length(start_of_read)

#make output tibble (one row)
result_tib = tibble(sample_name = sample_name, 
                    launch_pad = LP,
                    orientation = orientation,
                    total_reads = N_total_reads,
                    home_reads = N_home_reads,
                    fract_home_reads = N_home_reads/N_total_reads)

#write output as tsv
write_tsv(result_tib, file = output_file)
