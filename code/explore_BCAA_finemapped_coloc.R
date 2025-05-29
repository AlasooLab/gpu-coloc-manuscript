library("data.table")
library("dplyr")
ge_data = readr::read_tsv("data/coloc/ge_vs_met56_coloc_results.tsv")
lc_data = readr::read_tsv("data/coloc/leafcutter_vs_met56_coloc_results.tsv")

coloc_data = dplyr::bind_rows(ge_data, lc_data)

dplyr::filter(coloc_data, gene_name == "BCAT2") %>% View()

new_coloc = readr::read_tsv("data/met_56_lbfs.tsv")

dplyr::filter(new_coloc, signal1 %like% "Total_BCAA_chr19:47301513-50300951") %>% dplyr::filter(PP.H4 > 0.8) %>% View()

dplyr::filter(new_coloc, signal1 %like% "Total_BCAA_chr4:86804123-89803718") %>% dplyr::filter(PP.H4 > 0.8) %>% View()



dplyr::filter(new_coloc, signal1 %like% "HDL_C_chr16:55457585-58457455_L3") %>% dplyr::filter(PP.H4 > 0.8) %>% View()


split_signal = tidyr::separate(new_coloc, signal2, c("dataset_id", "rest"), sep = "_", extra = "merge") %>%
  tidyr::separate(rest, c("molecular_trait_id", "cs_index"), sep = "_L") %>%
  dplyr::mutate(cs_index = paste0("L", cs_index))

#Import leafcutter meta
file_paths = list.files("data/big_data/leafcutter_meta/", full.names = T)

dataset_labels = list.files("data/big_data/leafcutter_meta/") %>%
  stringr::str_remove("leafcutter_") %>%
  stringr::str_remove("_Ensembl_105_phenotype_metadata.tsv.gz")

file_list = setNames(as.list(file_paths), dataset_labels)
lc_meta = purrr::map_df(file_list, readr::read_tsv, .id = "dataset_id") %>%
  dplyr::transmute(dataset_id, molecular_trait_id = phenotype_id, gene_id, gene_name)

lc_meta2 = dplyr::mutate(lc_meta, molecular_trait_id = paste(dataset_id, molecular_trait_id, sep = "_"))

joined = dplyr::left_join(split_signal, lc_meta)

#LBFs
met_coloc = readr::read_tsv("data/met_lbfs_with_ensemble.tsv")
fingen_coloc = readr::read_tsv("data/finngen_lbf_with_ensemble.tsv")

#ABFs
met_abf_coloc = readr::read_tsv("data/met_abfs_with_ensemble.tsv")
fingen_abf_coloc = readr::read_tsv("data/finngen_abf_with_ensemble.tsv")


finemapped_results = readr::read_tsv("~/Downloads/56_metabolites_finemapping_credible_sets.tsv")


#BCAT2
dplyr::filter(met_coloc, signal1 %like% "Total_BCAA_chr19:47301513-50300951", PP.H4 > 0.9) %>% View()

#PPM1K
dplyr::filter(met_coloc, signal1 %like% "Total_BCAA_chr4:86804123-89803718") %>% dplyr::filter(PP.H4 > 0.9) %>% View()


