library("dplyr")
library("locuszoomr")
library("EnsDb.Hsapiens.v86")

#Extract coloc
met_coloc = readr::read_tsv("data/met_lbfs_with_ensemble.tsv")

#Import LDL
ds = arrow::open_dataset("../EstBB-UKBB-metaanalysis/data/gwas_sumstats/LDL_C_EstBB_UKBB_EUR_metaanalysis.parquet")
chr19_data = dplyr::filter(ds, CHROM == 19) %>% collect()

metabolite_df = dplyr::transmute(chr19_data, chrom = CHROM, pos = GENPOS, 
                                 other_allele = ALLELE0, effect_allele = ALLELE1, nlogp = LOG10P, beta = Effect, se = StdErr, variant = paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = "_")) %>%
  as.data.frame()
metabolite_df = dplyr::mutate(metabolite_df, rsid = variant)

loc <- locus(data = metabolite_df, gene = 'LDLR', flank = 1.5e5, ens_db = "EnsDb.Hsapiens.v86", yvar = "nlogp")
summary(loc)

pdf("figures/LDLR_marginal_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

#Import fine mapping results
ldlc_finemap = readr::read_tsv("../metabolite_GWAS/big_data/UKBB_sumstats_for_sharing/coloc5/LDL_C_coloc5_final.tsv.gz")
ldlr = dplyr::filter(ldlc_finemap, region == "chr19:9587925-12587668")

#Missense signal
gwas_df = dplyr::transmute(ldlr, chrom = chromosome, pos = position, rsid = variant, other_allele = "A", effect_allele = "T", lbf_variable1, variant = variant) %>%
  as.data.frame()
loc <- locus(data = gwas_df, gene = 'LDLR', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", yvar = 'lbf_variable1')

pdf("figures/LDLR_Rahu_missense_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

#Intron retention signal
gwas_df = dplyr::transmute(ldlr, chrom = chromosome, pos = position, rsid = variant, other_allele = "A", effect_allele = "T", lbf_variable2, variant = variant) %>%
  as.data.frame()
loc <- locus(data = gwas_df, gene = 'LDLR', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", yvar = 'lbf_variable2')

pdf("figures/LDLR_Rahu_IR_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

#Read LDLR eQTL data
bpT = readr::read_tsv("data/QTD000031.lbf_variable.txt.gz")
ldlr_eqtl = dplyr::filter(bpT, molecular_trait_id == "ENSG00000130164")

qtl_df = dplyr::transmute(ldlr_eqtl, chrom = chromosome, pos = position, rsid = variant, other_allele = "A", effect_allele = "T", lbf_variable1, variant = variant) %>%
  as.data.frame()

loc <- locus(data = qtl_df, gene = 'LDLR', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", yvar = 'lbf_variable1')
#loc$data = loc$data[!is.na(loc$data$variant),]

pdf("figures/LDLR_BLUEPRINT_eQTL_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

majiq = readr::read_tsv("data/majiq_irQTL.tsv.gz", col_names = F)
colnames(majiq) = colnames(bpT)
majiq_locus = dplyr::filter(majiq, molecular_trait_id == "ENSG00000130164:t:11123174-11123344:11120523-11123173")

qtl_df = dplyr::transmute(majiq_locus, chrom = chromosome, pos = position, rsid = variant, other_allele = "A", effect_allele = "T", lbf_variable1, variant = variant) %>%
  as.data.frame()

loc <- locus(data = qtl_df, gene = 'LDLR', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", yvar = 'lbf_variable1')
#loc$data = loc$data[!is.na(loc$data$variant),]

pdf("figures/LDLR_GTEx_LCL_irQTL_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()


