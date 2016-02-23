library(snowfall);
library(sqldf);
library(snpStats);

# upstream;downstream
# UTR5;UTR3
# exonic;splicing
filter_comma<-function( newtb )
{
	newtbm_nes<-sqldf("select * from newtb where antype!='exonic;splicing' and antype!='UTR5;UTR3' and antype!='upstream;downstream'");
	newtbm_es<-sqldf("select * from newtb where antype='exonic;splicing' or antype='UTR5;UTR3' or antype='upstream;downstream'");
	newtb0<-sqldf("select * from newtbm_es where genes not like '%;%'");
	newtb1<-sqldf("select * from newtbm_es where genes like '%;%'");
	newtbx <- c();
	for(i in 1:dim(newtb1)[1])
	{
		genes<-strsplit( as.character(newtb1[i,2]), ";");
		temp <- data.frame(rep(newtb1[i,],1), V7=c(1:length(genes[[1]])) )
		temp[,2] <- genes[[1]]
		genes<-strsplit( as.character(newtb1[i,1]), ";");
		temp[,1] <- genes[[1]]
		newtbx <- rbind( newtbx, temp);
	}


	newtbm_es <- rbind(newtb0, newtbx[,-9] );
	tbr <- rbind(newtbm_nes, newtbm_es);
	return(tbr);
}

# exonic GENEA,GENEB
# ncRNA_exonic
filter_exonic<-function( newtb )
{
	newtbm_exonic<-sqldf("select * from newtb where antype='exonic' or antype='ncRNA_exonic'");
	newtb0<-sqldf("select * from newtbm_exonic where genes not like '%,%'");
	newtb1<-sqldf("select * from newtbm_exonic where genes like '%,%'");
	newtbx <- c();
	for(i in 1:dim(newtb1)[1])
	{
		#cat("Row", i, "\n");
		genes<-strsplit( as.character(newtb1[i,2]), ",");
		temp <- data.frame(rep(newtb1[i,],1), V7=c(1:length(genes[[1]])) )
		temp[,2] <- genes[[1]]
		newtbx <- rbind( newtbx, temp);
	}

	newtbm_exonic <- rbind(newtb0, newtbx[,-9]);
	return(newtbm_exonic);
}


# intronic
# ncRNA_intronic
filter_intronic<-function( newtb )
{
	newtbm_intronic<-sqldf("select * from newtb where antype='intronic' or antype='ncRNA_intronic'");
	newtb0<-sqldf("select * from newtbm_intronic where genes not like '%,%'");
	newtb1<-sqldf("select * from newtbm_intronic where genes like '%,%'");
	newtbx <- c();
	for(i in 1:dim(newtb1)[1])
	{
		#cat("Row", i, "\n");
		genes<-strsplit( as.character(newtb1[i,2]), ",");
		temp <- data.frame(rep(newtb1[i,],1), V7=c(1:length(genes[[1]])) )
		temp[,2] <- genes[[1]]
		newtbx <- rbind( newtbx, temp);
	}

	newtbm_intronic <- rbind(newtb0, newtbx[,-9]);
	return(newtbm_intronic)
}


# splicing
# ncRNA_splicing
filter_splicing<-function( newtb )
{
	newtbm_splicing<-sqldf("select * from newtb where antype='splicing' or antype='ncRNA_splicing'");
	newtb0<-sqldf("select * from newtbm_splicing where genes not like '%(%'");
	newtb1<-sqldf("select * from newtbm_splicing where genes like '%(%'");

	temp.s1<-strsplit(as.character(newtb1$genes), "\\(");
	temp <- c();
	for(i in 1:dim(newtb1)[1])
	{
		temp <- c(temp, temp.s1[[i]][1]);
	}
	newtb1$genes <-temp

	newtbm_splicing <- rbind(newtb0, newtb1 );
	return(newtbm_splicing)
}

# ncRNA_UTR3
# ncRNA_UTR5
# UTR3
# UTR5
filter_utr35<-function( newtb )
{
	newtbm_utr<-sqldf("select * from newtb where antype='ncRNA_UTR3' or antype='ncRNA_UTR5' or antype='UTR3' or antype='UTR5'");

	newtb0<-sqldf("select * from newtbm_utr where genes not like '%(%'");
	newtb1<-sqldf("select * from newtbm_utr where genes like '%(%'");
	temp.s1<-strsplit(as.character(newtb1$genes), "\\(");
	temp <- c();
	for(i in 1:dim(newtb1)[1])
	{
		temp <- c(temp, temp.s1[[i]][1]);
	}
	newtb1$genes <-temp
	newtbm_utr <- rbind(newtb0, newtb1 );
	
	return(newtbm_utr);
}

# downstream
# upstream
filter_updown<-function( newtb )
{
	newtbm_updown<-sqldf("select * from newtb where antype='downstream' or antype='upstream'");

	newtb0<-sqldf("select * from newtbm_updown where genes not like '%,%'");
	newtb1<-sqldf("select * from newtbm_updown where genes like '%,%'");
	newtbx <- c();
	for(i in 1:dim(newtb1)[1])
	{
		#cat("Row", i, "\n");
		genes<-strsplit( as.character(newtb1[i,2]), ",");
		temp <- data.frame(rep(newtb1[i,],1), V7=c(1:length(genes[[1]])) )
		temp[,2] <- genes[[1]]
		newtbx <- rbind( newtbx, temp);
	}

	newtbm_updown <- rbind(newtb0, newtbx[,-9]);
	
	return(newtbm_updown);
}

gene_sorting <- function( file.annovar.variant, file.gen.set )
{
	tb<-read.table( file.annovar.variant, sep="\t", header=F)

	write.table(tb$V3, "temp_file_annovar_variant.txt", quote=F, row.names=F, col.names=F)
	tbv3<-read.table("temp_file_annovar_variant.txt", sep=" ", header=F)
	newtb <- cbind(antype=tb$V1, genes=tb$V2, tbv3)

	#remove intergenic
	newtb.intergenic<-newtb[which(newtb$antype=="intergenic"),]
	newtb <- newtb[-which(newtb$antype=="intergenic"),]

cat("Search genes including COMMA...\n");
	newtb           <- filter_comma( newtb);
cat("Search genes including EXONIC...\n");	
	newtbm_exonic   <- filter_exonic( newtb);
cat("Search genes including INTRONIC...\n");	
	newtbm_intronic <- filter_intronic( newtb);
cat("Search genes including SPLICING...\n");	
	newtbm_splicing <- filter_splicing( newtb);
cat("Search genes including UTR35...\n");	
	newtbm_utr35    <- filter_utr35( newtb);
cat("Search genes including UPDOWN...\n");	
	newtbm_updown   <- filter_updown( newtb);

cat("Merge UPDOWN, UTR35, SPLICING, INTRONIC, EXONIC...\n");	
	newtbn <- rbind( newtbm_updown, newtbm_utr35, newtbm_splicing, newtbm_intronic, newtbm_exonic);
	colnames(newtbn) <- c("antype", "genes", "chr", "start", "stop", "refer", "alter", "snp");

	if(1)
	{
		idx <- which(newtbn[,2]=="LIMS3,LIMS3L");
		if (length(idx)>0)
			newtbn<- newtbn[-idx,]

		idx <- which(newtbn[,2]=="C1QTNF5,MFRP");
		if (length(idx)>0)
			newtbn <- rbind(newtbn, newtbn[idx[1],])
		idx <- which(newtbn[,2]=="C1QTNF5,MFRP");
		if (length(idx)>0)
		{
			newtbn[idx[1],2]<-"C1QTNF5";
			newtbn[idx[2],2]<-"MFRP";
		}
		idx <- which(newtbn[,2]=="SNRPN,SNURF");
		if (length(idx)>0)
			newtbn <- rbind(newtbn, newtbn[idx[1],])
		idx <- which(newtbn[,2]=="SNRPN,SNURF");
		if (length(idx)>0)
		{
			newtbn[idx[1],2]<-"SNRPN";
			newtbn[idx[2],2]<-"SNURF";
		}	
	}


 	newtb.intergenic[,2]<-"NONE"
	colnames(newtb.intergenic) <- c("antype", "genes", "chr", "start", "stop", "refer", "alter", "snp");
	snp_all <- rbind( newtbn, newtb.intergenic);

cat("Write GENE Set into ", file.gen.set, "...\n");	
	snp_gene <- sqldf("select * from newtbn order by chr, start");
	write.table(newtbn[,c(2,8)], file=file.gen.set, sep=" ", row.names=F, quote=F, col.names=F) 

	tb.gene   <- sqldf("select genes,count(genes) from newtbn group by genes");
	tb.antype <- sqldf("select antype,count(antype) from newtbn group by antype");

cat("Write summary into ", file.gen.set, "XXX.gen.sum.csv and XXX.antype.sum.csv...\n");	
	write.csv(tb.gene,   file=paste(file.gen.set, ".gen.sum.csv", sep=""), quote=F ) 
	write.csv(tb.antype, file=paste(file.gen.set, ".antype.sum.csv", sep=""), quote=F ) 
	
	return( list(snp_gene=snp_gene, snp_all=snp_all) );
}

gene_range_sorting <- function(snp_all, n.cpu=7 )
{
	tb_genes <- sqldf("select genes, count(genes) as snps, min(chr) as chr, max(chr), min(start) as minp, max(start) as maxp from snp_all group by genes");
	tb_genes2 <- sqldf("select * from tb_genes order by chr desc, minp");
	
	n.times <- floor( NROW(tb_genes)/1000 ) +1;
cat("ROW(gene)=", NROW(tb_genes2), "Times=", n.times, "\n");

	boot.fun <- function(sect)
	{
		library(sqldf);

		n.start <- 1+sect*1000;
		n.stop <- (1+sect)*1000;
		if (n.stop>NROW(tb_genes2)) n.stop <- NROW(tb_genes2);

		snp_ext<-c();
		for(i in n.start:n.stop)
		{
			chr <- tb_genes2[i,3];
			minp<- tb_genes2[i,5];
			maxp<- tb_genes2[i,6];

			sql <- paste("select * from snp_all where chr='", chr, "' and start>=", minp-100*1000, " and start<=", maxp+100*1000, sep="" );
			genes.tb <- sqldf(sql);
			genes.new <- cbind( biggene=tb_genes2[i,1], genes.tb);
			snp_ext <- rbind( snp_ext, genes.new );

			if( i %% 100 == 0) 
				cat("row=", i, "\n");
		}
		
		return(snp_ext);
	}
	
	sfInit(parallel=TRUE, cpus=n.cpu, type="SOCK")

	sfExport("tb_genes2", "snp_all");

	res <- sfClusterApplyLB(0:(n.times-1), boot.fun)
	
	snp_ext <- c();
	for(i in 1:length(res))
	{
		cat("SECT=", i, "NROW=", NROW(res[[i]]));
		snp_ext <- rbind(snp_ext, res[[i]]);
	}
	
	cat("Completed!")
	
	return(snp_ext);
}


lskat_pl_annovar<-function( path.annovar, file.plink.bed, file.plink.bim, file.plink.fam, file.gene.tab, file.gen.setid )
{
	snp.mat <- read.plink( file.plink.bed,  file.plink.bim, file.plink.fam);
	snp.map <- cbind(snp.mat$map[,1], snp.mat$map[,4], snp.mat$map[,4], snp.mat$map[,5], snp.mat$map[,6], snp.mat$map[,2]);
	write.table(snp.map, file.gene.tab, sep=" ", quote=F, col.names=F, row.names=F);

	cmd.annovar <- c("cd ", path.annovar, "\n", "./annotate_variation.pl", "-out", file.gene.tab, "-build hg19 ", file.gene.tab, "humandb/", "\n");
	t1 <- try(system(paste(cmd.annovar, collapse=" "), intern = TRUE))
	
	show(t1);
	
	annotate_variation.tab<- paste(file.gene.tab, "variant_function", sep=".");
	if (file.exists( annotate_variation.tab ))
	{
		t2 <- gene_sorting( annotate_variation.tab , file.gen.setid);
		save( t2, file=paste( file.gen.setid, ".sort.rdata", sep="") )
		
		t3 <- gene_range_sorting( t2$snp_all);
		save( t3, file=paste( file.gen.setid, ".geninfo.rdata", sep="") );
	}
}

# !!! http://www.openbioinformatics.org/annovar/annovar_gene.html
# annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene humandb/
#
# cat FHS-bmi-gene.tab
# chr,pos0,pos1, A ,G, comment:  
# 1 948921 948921 T C comments: rs15842, a SNP in 5' UTR of ISG15
# 1 1404001 1404001 G T comments: rs149123833, a SNP in 3' UTR of ATAD3C
# 1 5935162 5935162 A T comments: rs1287637, a splice site variant in NPHP4
# 1 162736463 162736463 C T comments: rs1000050, a SNP in Illumina SNP arrays
#
# annotate_variation.pl -out FHS-bmi-gene.tab -build hg19 FHS-bmi-gene.tab humandb/
# NOTICE: The --geneanno operation is set to ON by default
# NOTICE: Reading gene annotation from humandb/hg19_refGene.txt ... Done with 48660 transcripts (including 10375 without coding sequence annotation) for 25588 unique genes
# NOTICE: Reading FASTA sequences from humandb/hg19_refGeneMrna.fa ... Done with 14 sequences
# WARNING: A total of 333 sequences will be ignored due to lack of correct ORF annotation
# NOTICE: Finished gene-based annotation on 15 genetic variants in example/ex1.avinput
# NOTICE: Output files were written to FHS-bmi-gene.tab.variant_function, FHS-bmi-gene.tab.exonic_variant_function

path.annovar <-"/home/zw224/g/annovar"
lskat_pl_annovar( path.annovar, 
		"/home/zw224/f/bmi/FHS-bmi-v1.bed" , 
		"/home/zw224/f/bmi/FHS-bmi-v1.bim" , 
		"/home/zw224/f/bmi/FHS-bmi-v1.fam" , 
		"/home/zw224/f/bmi/FHS-bmi-gene.tab", 
		"/home/zw224/f/bmi/FHS-bmi-gen-hg19.SetID");
