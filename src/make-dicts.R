## make dict/named list of gene:snp

dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/EGEUV1/quant/refbias2/Btrecase/SpikeMixV3_2/GT"
files <- list.files(dir, full.names = TRUE)
files <- grep("rbias", files, value=TRUE)
stan_files <- grep("GT.stan1.input.rds", files, value = TRUE, fixed = TRUE)
genes <- unique(gsub(".*(ENSG\\d+)\\..*", "\\1", stan_files))

snp_list <- lapply(stan_files, function(x) names(readRDS(x)))
names(snp_list) <- genes

write(rjson::toJSON(snp_list), "GT_dict.json")




dir <- "/home/abo27/rds/rds-mrc-bsu/ev250/psoriasis/refbias/Btrecase/SpikePrior/fisher001/rna/"
files <- list.files(dir, full.names = TRUE)

ns <- files[grep("refbias.ENSG\\d+\\.normal_skin.noGT.stan.input.rds", files)]
ps <- files[grep("refbias.ENSG\\d+\\.Psoriasis_skin.noGT.stan.input.rds", files)]

ngenes <- gsub(".*(ENSG\\d+).*", "\\1", ns)
pgenes <- gsub(".*(ENSG\\d+).*", "\\1", ps)
genes <- intersect(ngenes, pgenes)

snp_list <- lapply(genes, function(gene) {
    nsnps <- names(readRDS(grep(gene, ns, value = TRUE)))
    psnps <- names(readRDS(grep(gene, ps, value = TRUE)))
    intersect(nsnps, psnps)
})
names(snp_list) <- genes
write(rjson::toJSON(snp_list), "noGT_dict.json")
