log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(polyester)

args <- c(snakemake@input, snakemake@params, snakemake@output)
fun_args <- c(
    formalArgs(simulate_experiment), 
    "readlen", 
    "lib_sizes", 
    "distr", 
    "error_model", 
    "bias", 
    "gcbias", 
    "frag_GC_bias", 
    "strand_specific", 
    "meanmodel", 
    "write_info", 
    "seed", 
    "transcriptid", 
    "gzip", 
    "exononly", 
    "idfield", 
    "attrsep")
args <- args[names(args) %in% fun_args]

# Reshape fold_changes to matrix
groups <- length(args[["num_reps"]])
args[["fold_changes"]] <- matrix(args[["fold_changes"]], ncol = groups)

do.call(simulate_experiment, args)
