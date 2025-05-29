library("dplyr")
library("locuszoomr")
library("data.table")
library("EnsDb.Hsapiens.v86")

#BCAT2
met_coloc = readr::read_tsv("data/met_lbfs_with_ensemble.tsv")
dplyr::filter(met_coloc, signal1 %like% "Total_BCAA_chr19:47301513-50300951", PP.H4 > 0.9) %>% View()


#Import fine-mapped GWAS
bcat2_finemap = readr::read_tsv("data/Total_BCAA_coloc5_final.tsv.gz")
bcat2_region = dplyr::filter(bcat2_finemap, region == "chr19:47301513-50300951")

#sQTL singal
gwas_df = dplyr::transmute(bcat2_region, chrom = chromosome, pos = position, rsid = variant, other_allele = "A", effect_allele = "T", lbf_variable2, variant = variant) %>%
  as.data.frame()
loc <- locus(data = gwas_df, gene = 'BCAT2', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", yvar = 'lbf_variable2')

pdf("figures/BCAT2_met_sQTL_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()


#eQTL singal
gwas_df = dplyr::transmute(bcat2_region, chrom = chromosome, pos = position, rsid = variant, other_allele = "A", effect_allele = "T", lbf_variable4, variant = variant) %>%
  as.data.frame()
loc <- locus(data = gwas_df, gene = 'BCAT2', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", yvar = 'lbf_variable4')

pdf("figures/BCAT2_met_eQTL_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

#Import GEUVADIS sQTL
geuv_lbf = readr::read_tsv("data/QTD000114.lbf_variable.txt.gz")
bcat2_sqtl_lbf = dplyr::filter(geuv_lbf, molecular_trait_id == "19:48800297:48807000:clu_9896_-")
write.table(bcat2_sqtl_lbf, "data/bcat2_sqtl_lbf.tsv", sep = "\t", row.names = F, quote = F)

qtl_df = dplyr::transmute(bcat2_sqtl_lbf, chrom = chromosome, pos = position, rsid = variant, other_allele = "A", effect_allele = "T", lbf_variable1, variant = variant) %>%
  as.data.frame()
loc <- locus(data = qtl_df, gene = 'BCAT2', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", yvar = 'lbf_variable1', index_snp = "chr19_48806519_G_C")

pdf("figures/BCAT2_GEUV_sQTL_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

#Import FUSION eQTL
fusion_lbf = readr::read_tsv("data/QTD000090.lbf_variable.txt.gz")
bcat2_eqtl_lbf = dplyr::filter(fusion_lbf, molecular_trait_id == "ENSG00000105552")
write.table(bcat2_eqtl_lbf, "data/bcat2_eqtl_lbf.tsv", sep = "\t", row.names = F, quote = F)

qtl_df = dplyr::transmute(bcat2_eqtl_lbf, chrom = chromosome, pos = position, rsid = variant, other_allele = "A", effect_allele = "T", lbf_variable1, variant = variant) %>%
  as.data.frame()
loc <- locus(data = qtl_df, gene = 'BCAT2', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", yvar = 'lbf_variable1', index_snp = "chr19_48797174_G_A")

pdf("figures/BCAT2_FUSION_eQTL_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()
