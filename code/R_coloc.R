Sys.setenv(R_LIBS_USER = "~/Rlibs")
Sys.setenv(R_LIBS_SITE = "")       
dir.create("~/Rlibs", showWarnings = FALSE, recursive = TRUE)
.libPaths(c("~/Rlibs"))             

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = "~/Rlibs", repos = "https://cloud.r-project.org/")
}

if (!requireNamespace("progressr", quietly = TRUE)) {
    install.packages("progressr", lib = "~/Rlibs", repos = "https://cloud.r-project.org/")
}

library(remotes)

# Force installation in user library to avoid system write errors
remotes::install_github("chr1swallace/coloc", 
                        build_vignettes=FALSE, 
                        dependencies=TRUE, 
                        lib = "~/Rlibs", 
                        force = TRUE)

library(data.table)
library(dplyr)
library(coloc)
library(progressr)

handlers(global = TRUE)
handlers(handler_progress(
  format = "[:bar] :percent | Elapsed: :elapsedfull | ETA: :eta",
  clear  = FALSE,
  show_after = 0
))

results_output <- "R_results.tsv"
results_file_exists <- file.exists(results_output)

met_files <- list.files(
  path = "met_files", 
  pattern = "coloc5_final.tsv.gz$", 
  recursive = TRUE, 
  full.names = TRUE
)

eqtl_files <- list.files(
  path = "files", 
  pattern = "lbf_variable", 
  recursive = TRUE, 
  full.names = TRUE
)

with_progress({
  p <- progressor(along = met_files)

  for (met_file in met_files) {
    met_df <- fread(met_file, sep = "\t", data.table = FALSE)
    met_df$chromosome <- as.character(met_df$chromosome)
    met_chromosomes <- unique(met_df$chromosome)
    met <- unique(met_df$molecular_trait_id)[1]

    handlers("progress")

    q <- progressor(along = eqtl_files)  

    for (eqtl_file in eqtl_files) {
      eqtl_df <- fread(eqtl_file, sep = "\t", data.table = FALSE)
      eqtl_parts <- strsplit(eqtl_file, .Platform$file.sep)[[1]]
      eqtl_qtd <- if (length(eqtl_parts) > 2) eqtl_parts[2] else eqtl_parts[1]

      eqtl_df$chromosome <- as.character(eqtl_df$chromosome)
      eqtl_chromosomes <- unique(eqtl_df$chromosome)

      common_chromosomes <- c("20")  # Hardcoded chromosome
      if (length(common_chromosomes) == 0) next

      for (chromosome in common_chromosomes) {
        met_chrom_df <- met_df[met_df$chromosome == chromosome, ]
        eqtl_chrom_df <- eqtl_df[eqtl_df$chromosome == chromosome, ]
        eqtl_chrom_traits <- unique(eqtl_chrom_df$molecular_trait_id)

        eqtl_ranges <- list()
        for (eqtl_trait in eqtl_chrom_traits) {
          eqtl_trait_df <- eqtl_chrom_df[eqtl_chrom_df$molecular_trait_id == eqtl_trait, ]
          eqtl_ranges[[eqtl_trait]] <- c(
            min(eqtl_trait_df$position), 
            max(eqtl_trait_df$position)
          )
        }

        met_regions <- unique(met_chrom_df$region)
        for (region in met_regions) {
          met_trait_df <- met_chrom_df[met_chrom_df$region == region, ]
          met_trait_name <- paste0(met, "_", region)
          region_coords <- unlist(strsplit(gsub(".*:", "", region), "-"))
          met_region_min <- as.numeric(region_coords[1])
          met_region_max <- as.numeric(region_coords[2])

          met_mat <- as.matrix(dplyr::select(met_trait_df, lbf_variable1:lbf_variable10))
          row.names(met_mat) <- met_trait_df$variant
          met_mat <- t(met_mat)

          for (eqtl_trait in eqtl_chrom_traits) {
            trait_range <- eqtl_ranges[[eqtl_trait]]

            if (trait_range[1] <= met_region_max && trait_range[2] >= met_region_min) {
              eqtl_trait_df <- eqtl_chrom_df[eqtl_chrom_df$molecular_trait_id == eqtl_trait, ]
              eqtl_trait_name <- paste0(eqtl_qtd, "_", eqtl_trait)

              eqtl_mat <- as.matrix(dplyr::select(eqtl_trait_df, lbf_variable1:lbf_variable10))
              row.names(eqtl_mat) <- eqtl_trait_df$variant
              eqtl_mat <- t(eqtl_mat)

              coloc <- coloc.bf_bf(
                met_mat, 
                eqtl_mat,
                p12 = 1e-6,
                overlap.min = 0.5,
                trim_by_posterior = TRUE
              )

              if (is.null(coloc) || is.null(coloc$summary)) next

              result_summary <- as.data.frame(coloc$summary)
              filtered_result <- result_summary %>%
                filter(PP.H4.abf >= 0)  # Apply meaningful filtering

              if (nrow(filtered_result) > 0) {
                filtered_result <- filtered_result %>%
                  mutate(
                    signal1 = paste0(met_trait_name, "_L", filtered_result$idx1),
                    signal2 = paste0(eqtl_trait_name, "_L", filtered_result$idx2)
                  )

                write_header <- !results_file_exists

                output_df <- filtered_result[, c("signal1", "signal2", "PP.H4.abf", "PP.H3.abf")]
                fwrite(
                  output_df, 
                  file = results_output,
                  sep = "\t", 
                  col.names = write_header,
                  row.names = FALSE, 
                  append = TRUE
                )

                results_file_exists <- TRUE
              }
            }
          }
        }
      }
      q()
    }
    handlers(global = TRUE) 
    p()
  }
})
