new.proj <- function( plink.path, main.path, plink.bfile.nobed, file.pheno.mean, file.pheno.long, cols.pheno, cols.time, cols.cov )
{
	plink.obj <- list();
	
	plink.obj$plink.path <- plink.path;
	plink.obj$main.path  <- main.path;

	plink.obj$genotype <- list();
	plink.obj$genotype$plink.bfile.nobed <- plink.bfile.nobed;
	plink.obj$genotype$file.plink.bed <- paste( plink.bfile.nobed, "bed", sep=".")
	plink.obj$genotype$file.plink.bim <- paste( plink.bfile.nobed, "bim", sep=".")
	plink.obj$genotype$file.plink.fam <- paste( plink.bfile.nobed, "fam", sep=".")

	plink.obj$phenotype <- list();
	plink.obj$phenotype$file.phe.mean <- file.pheno.mean;
	plink.obj$phenotype$file.phe.long <- file.pheno.long;
	plink.obj$phenotype$cols.pheno    <- cols.pheno;
	plink.obj$phenotype$cols.time     <- cols.time; 
	plink.obj$phenotype$cols.cov      <- cols.cov; 
	
	return(plink.obj);
}

