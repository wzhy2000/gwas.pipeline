library(gwas.pipeline);

setwd("/ycga-ba/home/zw224/fr/bmi4-test");

plink            <- "/home/bioinfo/software/plink/plink";
mainDir          <- "/ycga-ba/home/zw224/fr/bmi4-test";
plink.bfile      <- "testdata"
file.pheno.long  <- "bmi4-phe-c1c2-unrelated.csv"
phe.cols.sex     <- 2;
phe.cols.Y       <- c(25:46); 
phe.cols.Z       <- c(3:24);

p0 <- new.proj( plink, mainDir, plink.bfile, "", file.pheno.long, phe.cols.Y,phe.cols.Z, phe.cols.sex );

p0 <- do.plink.qc1(p0, options=list(p.mind=0.1, p.geno=0.1));

p0 <- do.post.qc1(p0, options=NULL);

p0 <- do.plink.qc2(p0, options=NULL);

p0 <- do.post.qc2(p0, options=list(mds.cov=5));

p0 <- do.plink.gwas(p0, options=NULL );

p0$gene$file.gene.hg19 <- "/ycga-ba/home/zw224/fr/bmi4/bmi-c1c2-hg19.SetID"
p0 <- do.SKAT(p0, options=list(covariate.count=7) );

p0 <- do.LSKAT(p0, options=list(covariate.count=7) );

p0 <- do.blasso(p0, options=list(covariate.count=7));

p0 <- do.glasso(p0, options=list(covariate.count=7));

save(p0, file="bmi4.rdata");
