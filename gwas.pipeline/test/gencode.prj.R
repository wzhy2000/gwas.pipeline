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

p0 <- do.gencode( p0, "gencode.v19.annotation.gtf", "bmi-c1c2-hg19-1.5k.SetID", ext_bp_size=1500, n.cores=7);


# the fowllowing results are ext_bp_size=50K
# 
# NROW(gencode_transcript_ext)
# 47972

# tb<-read.table( "bmi-c1c2-hg19.SetID" );
# dim(tb)
#[1] 1008886       2
# length(unique(tb[,1]))
#[1] 44012
# length(unique(tb[,2]))
#[1] 414317


# q<-sqldf("select V1, count(V1) as cnt from tb group by V1 order by cnt desc")
# head(q, n=20)
#             V1  cnt
#1         CSMD1 1023
#2        RBFOX1  839
#3         PTPRD  688
#4  RP11-420N3.2  670
#5         CDH13  649
#6          FHIT  596
#7       CNTNAP2  546
#8          WWOX  535
#9       MACROD2  530
#10         DAB1  520
#11        LSAMP  477
#12        LRP1B  418
#13        ASIC2  416
#14       CTNNA2  406
#15         DLG2  402
#16       CTNNA3  396
#17        MAGI2  386
#18        PRKG1  373
#19       PCDH15  361
#20         NRG1  358
