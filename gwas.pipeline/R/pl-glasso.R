do.glasso<-function( plink.obj, options=list(qc.method="qc2")  )
{
	library(gwas.lasso)
	
	subDir  <- "glasso"
	mainDir <- plink.obj$main.path;
	dir.create( file.path(mainDir, subDir), showWarnings = FALSE );

	file.phe.glasso <- "glasso/bmi-phe-glasso.csv"

	n.cov   <- options$covariate.count;
	tb <- read.csv( plink.obj$phenotype$file.phe.long, header=T);
	tbcov <- read.csv( plink.obj$qc2$file.pca.cov, header=T);
	tbcov <- tbcov[, c(2:(n.cov+2)), drop=F];
	colnames( tbcov ) <- c("shareid", paste("COV", 1:n.cov, sep="") );

	tb <- tb [ match( tbcov$shareid, tb[,1]), ];
	tb.names <- colnames(tb)
	tb.names[1] <- "shareid";
	tb.names[ plink.obj$phenotype$cols.pheno ] <- paste("Y_", 1:length(plink.obj$phenotype$cols.pheno), sep="");
	tb.names[ plink.obj$phenotype$cols.time ]  <- paste("Z_", 1:length(plink.obj$phenotype$cols.time), sep="");

	colnames(tb) <- tb.names;
	newphe <- cbind(tb, tbcov[,-1, drop=F])
	write.csv( newphe, file=file.phe.glasso , row.names=F, quote=F);

	show(head(newphe));

	plink.bfile <- plink.obj$genotype$plink.bfile.nobed;
	if(options$qc.method=="qc2") plink.bfile <- plink.obj$qc2$plink.out.bfile
	if(options$qc.method=="impute") plink.bfile <- plink.obj$impute$plink.out.bfile
		
	file.plink.bed <- paste( plink.bfile, "bed", sep="." );
	file.plink.bim <- paste( plink.bfile, "bim", sep="." );
	file.plink.fam <- paste( plink.bfile, "fam", sep="." );
	file.ret.rdata <- "glasso/glasso-ret.rdata";

	ret <- gls.plink( file.phe.glasso, 
				file.plink.bed, 
				file.plink.bim, 
				file.plink.fam, 
				"Y", 
				"Z", 
				paste("COV", 1:n.cov, sep=""), 
				refit = TRUE, 
				add.used = T, 
				dom.used = T, 
				fgwas.filter = T,
				force.split=T,
				plink.command = plink.obj$plink.path,
				options=options )

	save(ret, file=file.ret.rdata);
	summary(ret);
	plot(ret);

	plink.obj$glasso <- list(rdata=file.ret.rdata);
	return( plink.obj );
}