do.plink.gwas <- function( plink.obj, options=list())
{
	subDir  <- "plink.gwas/"
	mainDir <- plink.obj$main.path;
	dir.create( file.path(mainDir, subDir), showWarnings = FALSE );
	
	options0 <- list();
	if(!is.null(options)) options0[names(options)] <- options;

	plink.path        <- plink.obj$plink.path;
	plink.qc2.bfile   <- plink.obj$qc2$plink.out.bfile;
	file.pca.cov      <- plink.obj$qc2$file.pca.cov;
	plink.out.log     <- paste( subDir, plink.obj$genotype$plink.bfile.nobed, "-R.log", sep=""); 

	cat("[ PLINK GWAS ...]\n");
	cat("[ PLINK GWAS ...]\n", file=plink.out.log, append = FALSE);
	
	plink.out.log10  <- paste( subDir, plink.obj$genotype$plink.bfile.nobed, ".log10", sep="");
	plink.out.line   <- paste( subDir, plink.obj$genotype$plink.bfile.nobed, ".line",  sep="");
	file.logis.ret   <- paste( plink.out.log10, ".assoc.ret", sep="");
	file.linear.ret  <- paste( plink.out.line,  ".assoc.ret", sep="");

	t1 <- plink_command(plink.path, c("--noweb", 
					"--bfile", plink.qc2.bfile, 
					"--logistic",
					"--covar",  file.pca.cov, 
					"--out", plink.out.log10 ));
	cat(t1, file=plink.out.log, append = TRUE);	

	draw_logistic_fig2 ( file.logis.ret, plink.out.log10, "FHS BMI analysis");
	
	r.top <- NULL;
	r.sig <- show_plink_sign( file.logis.ret, tag="ADD" );
	if( is.null( r.sig$sig ) )
	{
		r.top <- show_plink_topn( file.logis.ret, tag="ADD", top=10 );
		cat("The top 10 SNPs of Y = b0 + b1.ADD + b2.COV1 + b3.COV2 + ... + e\n");
		show( r.top );
	}
	else
	{
		cat("Thr significant SNPs of Y = b0 + b1.ADD + b2.COV1 + b3.COV2 + ... + e\n");
		show( r.sig );
	}
	
	t2 <- plink_command(plink.path, c("--noweb", 
					"--bfile", plink.qc2.bfile, 
					"--linear --genotypic ",
					"--covar",  file.pca.cov, 
					"--out", plink.out.line ));
					
	cat(t2, file=pl.file.log, append = TRUE);	

	draw_linear_QQ_fig ( file.linear.ret, plink.out.line, " ");
	draw_linear_man_fig( file.linear.ret, plink.out.line, " ");
	
	r.top.add <- NULL;
	r.top.dom <- NULL;
	r.top.geno <- NULL;
	r.sig.add <- show_plink_sign( file.linear.ret, tag="ADD" );
	r.sig.dom <- show_plink_sign( file.linear.ret, tag="DOMDEV" );
	r.sig.geno <- show_plink_sign( file.linear.ret, tag="GENO_2DF" );
	if( is.null( r.sig.add$sig ) || is.null( r.sig.dom$sig ) || is.null( r.sig.geno$sig ))
	{
		r.top.add  <- show_plink_topn( file.linear.ret, tag="ADD", top=10 );
		r.top.dom  <- show_plink_topn( file.linear.ret, tag="DOMDEV", top=10 );
		r.top.geno <- show_plink_topn( file.linear.ret, tag="GENO_2DF", top=10 );

		cat("The top 10 SNPs of Y = b0 + b1.ADD + b2.DOM + b3.COV1 + b4.COV2 + ... + e\n");
		cat("\n[ ADO ]\n");
		show( r.top.add );
		cat("\n[ DOM ]\n");
		show( r.top.dom );
		cat("\n[ GENO ]\n");
		show( r.top.geno );
	}
	else
	{
		cat("The significant SNPs of Y = b0 + b1.ADD + b2.DOM + b3.COV1 + b4.COV2 + ... + e\n");
		cat("\n[ ADO ]\n");
		show( r.sig.add );
		cat("\n[ DOM ]\n");
		show( r.sig.dom );
		cat("\n[ GENO ]\n");
		show( r.sig.geno );
	}

	plink.obj$plink.gwas <- list( 
			    r.sig       = r.sig,
				r.sig.add   = r.sig.add,
	 			r.sig.dom   = r.sig.dom, 
				r.sig.geno  = r.sig.geno,
				r.top       = r.top,
				r.top.add   = r.top.add,
	 			r.top.dom   = r.top.dom, 
				r.top.geno  = r.top.geno);
				
	return(plink.obj);
}	