library("dplyr")
library("locuszoomr")
library("EnsDb.Hsapiens.v86")


fingen_coloc = readr::read_tsv("data/finngen_lbf_with_ensemble.tsv")
dplyr::filter(fingen_coloc, gene_name == "CRP") %>% View()

#Read CRP sQTL data
liver = readr::read_tsv("data/QTD000270.cc.tsv.gz")
crp_sqtl = dplyr::filter(liver, molecular_trait_id == "1:159713640:159714007:clu_11689_-")

qtl_df = dplyr::transmute(crp_sqtl, chrom = chromosome, pos = position, rsid, other_allele = ref, effect_allele = alt, p = pvalue, beta, se, variant = rsid) %>%
  as.data.frame()

loc <- locus(data = qtl_df, gene = 'CRP', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs1800947")
#loc$data = loc$data[!is.na(loc$data$variant),]

pdf("figures/CRP_GTEx_liver_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()

#Read pneumonia data
pneumo = readr::read_tsv("data/finngen_R12_J10_PNEUMOBACT.gz")
colnames(pneumo)[1] = "chrom"
pneumo_region = dplyr::filter(pneumo, chrom == 1, pos > 158213648 & pos < 161213648)

gwas_df = dplyr::transmute(pneumo_region, chrom = chrom, pos = pos, rsid = rsids, other_allele = ref, effect_allele = alt, pvalue = pval, beta = beta, se = sebeta, variant = rsids) %>%
  as.data.frame()

loc <- locus(data = gwas_df, gene = 'CRP', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs1800947")
summary(loc)

pdf("figures/CRP_pneumonia_LocusZoom.pdf", width = 5.5, height = 4)
locus_plot(loc)
dev.off()
