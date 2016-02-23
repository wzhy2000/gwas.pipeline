import_gencode <-function( file.gencode.gtf, V3.type=NULL )
{
	if (!is.null(V3.type))
		# Get level 1 & 2 annotation (manually annotated) only: 
		# Get all "gene" lines: 
		awk.cmd <- paste( "awk '($3==\"", V3.type, "\" && $0~\"level (1|2);\" ){gsub( /\\\";?/, \"\", $10);gsub( /\\\";?/, \"\", $18);print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$18}' ", file.gencode.gtf, sep="")
	else
		awk.cmd <- paste( "awk '{gsub( /\\\";?/, \"\", $10);gsub( /\\\";?/, \"\", $18);print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$18}' ", file.gencode.gtf, sep="");
	
cat(awk.cmd, "\n");

	# for gzipped GTF file
	if(tools::file_ext(file.gencode.gtf)=="gz" )
		awk.cmd <- paste( "zcat ",file.gencode.gtf," | awk '($3==\"", V3.type, "\"){gsub( /\\\";?/, \"\", $10);gsub( /\\\";?/, \"\", $18);print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$18}' ", sep="");
		
	bigdf <- read.table( pipe(awk.cmd), header = F );
	colnames(bigdf) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "same", "gene_id", "gene_name");
	
	return(bigdf);
}

lskat_pl_gencode<-function( file.gencode.gtf, file.plink.bim, file.gen.setid, ext_bp_size=50000, n.cores=1 )
{
	gencode_transcript_ext <- try( import_gencode( file.gencode.gtf, "gene" ) );
   	if(is.null(gencode_transcript_ext) || class(gencode_transcript_ext)=="try-error")
   		stop("Gencode data can not be found in the GTF file specified by the parameter of file.gencode.gtf.");

   	cat("  Gencode file (", file.gencode.gtf, ") has been loaded.\n")
    cat("  For exon", NROW(gencode_transcript_ext), "items are selected from GENCODE dataset.\n");

	#gene_name <- unlist(lapply( strsplit(as.character(gencode_transcript_ext$gene_name), "[.]"), function(x) { return(x[1]) } ) );
	#gencode_transcript_ext$gene_name <- gene_name;
	
	gene_grp1 <- sqldf("select gene_name, V1, count(V1) as sects from gencode_transcript_ext where V2='HAVANA' group by gene_name, V1 order by sects desc");
	gene_grp <- sqldf("select gene_name, count(V1) as chrs from gene_grp1 group by gene_name order by chrs desc");
	
	havana.grp <- sqldf("select * from gencode_transcript_ext where V2='HAVANA'");

#> head(gene_grp, n=20)
#     gene_name chrs
#1          7SK    4
#2        Y_RNA    3
#3  AJ271736.10    2
#4      AKAP17A    2
#5        AMDP1    2
#6         ASMT    2
#7        ASMTL    2
#8    ASMTL-AS1    2
#9         CD99    2
#10      CD99P1    2
#11       CRLF2    2
#12      CSF2RA    2

	tb.bim <- read.table( file.plink.bim );

	gene.list <- mclapply(1:NROW(havana.grp), function(i){
	#gene.list <- lapply(1:100, function(i){
		chr.i <- substring( as.character(havana.grp$V1[i]), 4 )

		idx <- which( tb.bim$V1 == chr.i & tb.bim[,4] > havana.grp$V4[i]-ext_bp_size & tb.bim[,4]<havana.grp$V5[i] + ext_bp_size );
		
		if(length(idx)>0)
			return( cbind( as.character(havana.grp$gene_name[i]), as.character(tb.bim[idx,2])))
		else
			return(NULL);
	}, mc.cores=n.cores);

	df.setid <- do.call("rbind", gene.list);

	df <- df.setid[ !duplicated(df.setid),] ;
	
	write.table(df, file=file.gen.setid, quote=F, col.names=F, row.names=F);
	
	return(file.gen.setid);
}

do.gencode<-function( plink.obj, file.gtf, file.gene.set, options=list(ext_bp_size=50000,n.cores=7) )
{
	cat("[ Gencode ...]\n");

	library(parallel);
	
	file.plink.bim <- plink.obj$genotype$file.plink.bim;

	file.gene.hg19 <- lskat_pl_gencode( file.gtf, file.plink.bim, file.gene.set, options$ext_bp_size, options$n.cores);
	
	plink.obj$gene <- list( file.gene.hg19 = file.gene.hg19 );
	
	return(plink.obj);
}
