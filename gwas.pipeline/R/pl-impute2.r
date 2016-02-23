get_con_param<-function(parm.id)
{
        for (e in commandArgs())
        {
                ta = strsplit(e,"=", fixed=TRUE);
                if(! is.na( ta[[1]][2]))
                {
                        temp = ta[[1]][2];
                        if( ta[[1]][1] == parm.id) {
                                temp = as.character(temp);
                                return (temp);
                        }
                }
        }

        return(NA);
}

cpu.fun<-function( chr.i )
{
	setwd( path.setwd );
	
	file.log <- paste(file.bfile.name, "-R-chr", chr.i, ".log", sep="");
	
	file.ped.out <- paste(file.bfile.name, "-chr", chr.i, sep="");
	file.ped <- paste(file.ped.out, ".ped", sep="");
	file.map <- paste(file.ped.out, ".map", sep="");

	cmd.plink <- paste(path.plink, "--noweb --bfile ", file.bfile, "--chr", chr.i, "--recode", "--out", file.ped.out, sep=" ");

	if(!file.exists(file.ped))
	{
		t1 <- system( cmd.plink,  intern = TRUE );
		cat(t1, file=file.log, append=T, sep="\n");
	}

	file.gtool.gen    <- paste(file.bfile.name, "-chr", chr.i, ".ped.gen", sep="");
	file.gtool.sample <- paste(file.bfile.name, "-chr", chr.i, ".ped.sample", sep="");
	
	cmd.gtool <- paste(path.gtool, "-P  --ped", file.ped, "--map", file.map, sep=" ");
	if(!file.exists(file.gtool.gen))
	{
		t2 <- system( cmd.gtool,   intern = TRUE  );
		cat(t2, file=file.log, append=T, sep="\n");
	}
	
	tb.map <- read.table(file.map);
	min.pos <- min(tb.map$V4);
	max.pos <- max(tb.map$V4);
	options(scipen=30);
	range <- unique( seq(1, NROW(tb.map), 1000 ), NROW(tb.map) );
	range[length(range)] <- range[length(range)] + 1;

	file.map2  <- paste(path.impute2.data, "genetic_map_chr", chr.i, "_combined_b37.txt", sep="");
	file.haps  <- paste(path.impute2.data, "ALL_1000G_phase1integrated_v3_chr", chr.i, "_impute.hap.gz", sep="");
	file.legend<- paste(path.impute2.data, "ALL_1000G_phase1integrated_v3_chr", chr.i, "_impute.legend.gz", sep="");
	file.impute2 <- paste( file.gtool.gen, ".impute2", sep="");
	
	for(i in 1:(length(range)-1) )
	{
		cat(">>>>>>>>>>>>>>>>>>>>>>>>>--", i, "--<<<<<<<<<<<<<<<<<<<<<<<\n");

		file.seq.prephase<- paste( file.gtool.gen, ".seq", i, ".prephasing", sep="");
		file.seq.prephase.haps <- paste( file.gtool.gen, ".seq", i, ".prephasing_haps", sep="");
		file.seq.impute2 <- paste( file.gtool.gen, ".seq", i, ".impute2", sep="");
		
		keep.snp<-tb.map$V2[ range[i] : (range[i+1]-1) ];
		keep.file <- tempfile()
		write.table(keep.snp, file=keep.file, quote=F, col.names=F, row.names=F)
		
		cmd.prephasing <- paste(path.impute2, 
					 " -prephase_g ", 
					 " -m ", file.map2,
					 " -g ", file.gtool.gen,
					 " -int ", tb.map$V4[ range[i] ] , tb.map$V4[ range[i+1] ] - 1,
					 " -Ne 20000 ",		 
					 " -allow_large_regions",
					 " -include_snps", keep.file, 
					 " -o ", file.seq.prephase, sep=" ");

		if(!file.exists(file.seq.impute2))
		{
			cat("CALL:", cmd.prephasing, file=file.log, append=T, sep="\n");
			t3 <- system(cmd.prephasing, intern = TRUE  );
			cat(t3, file=file.log, append=T, sep="\n");
		}
		
		cmd.impute2 <- paste(path.impute2, 
				 " -use_prephased_g", 
				 " -m ", file.map2,
				 " -h ", file.haps,
				 " -l ", file.legend,
				 " -known_haps_g ", file.seq.prephase.haps,
				 " -int ", tb.map$V4[ range[i] ] , tb.map$V4[ range[i+1] ] - 1,
				 " -allow_large_regions",
				 " -Ne 20000 ",	
				 " -include_snps", keep.file, 
				 " -o ",file.seq.impute2, sep=" "  );

		if(!file.exists(file.seq.impute2))
		{
			cat("CALL:", cmd.impute2, file=file.log, append=T, sep="\n");
			t4<- system(cmd.impute2, intern = TRUE );
			cat(t4, file=file.log, append=T, sep="\n");

			file.imputed.ped <- paste(file.gtool.gen, ".seq", i, ".impute2.ped", sep="");
			file.imputed.map <- paste(file.gtool.gen, ".seq", i, ".impute2.map", sep="");


			cmd.gtool2 <- paste(path.gtool, "-G --g", file.seq.impute2, "--s", file.gtool.sample, "--ped", file.imputed.ped, "--map", file.imputed.map, sep=" ");
			cat(cmd.gtool2, file=file.log, append=T, sep="\n");

			if(!file.exists(file.imputed.ped))
			{
				t5 <- system( cmd.gtool2, intern = TRUE  );
				cat(t5, file=file.log, append=T, sep="\n");
			}
		}
	}

	cat(">>>>>>>>>>>>>>>>>>>>>>>>> D O N E <<<<<<<<<<<<<<<<<<<<<<<\n");

	file.seq.plink <- c();
	if(!file.exists(file.impute2))
	{
		i <- 1;

		file.imputed.ped  <- paste(file.gtool.gen, ".seq", i, ".impute2.ped", sep="");
		file.imputed.map  <- paste(file.gtool.gen, ".seq", i, ".impute2.map", sep="");
		file.imputed.new1 <- paste(file.gtool.gen, ".seq", i, ".impute2", sep="");

		if(file.exists(file.imputed.ped))
		{
		}
		else
		{
			file.imputed.new1 <- paste(file.gtool.gen, ".seq", i, sep="");
			cmd.plink <- paste(path.plink, "--noweb --file ", file.ped.out, "--chr", chr.i, "--from-bp", tb.map$V4[ range[1] ], "--to-bp", tb.map$V4[ range[2] ] - 1, "--recode --out", file.imputed.new1, sep=" ");
			system( cmd.plink );
		}
		
		for(i in 2:(length(range)-1) )
		{
			file.imputed.ped <- paste(file.gtool.gen, ".seq", i, ".impute2.ped", sep="");
			file.imputed.map <- paste(file.gtool.gen, ".seq", i, ".impute2.map", sep="");
			
			if(file.exists(file.imputed.ped))
			{
				file.seq.plink <- c( file.seq.plink, paste( file.imputed.ped, file.imputed.map, sep=" ") )
			}
			else
			{
				file.imputed.new <- paste(file.gtool.gen, ".seq", i, sep="");
				
				cmd.plink <- paste(path.plink, "--noweb --file ", file.ped.out, "--chr", chr.i, "--from-bp", tb.map$V4[ range[i] ], "--to-bp", tb.map$V4[ range[i+1] ]-1, "--recode --out", file.imputed.new, sep=" ");
				system( cmd.plink );
			
				file.imputed.ped <- paste(file.imputed.new, ".ped", sep="");
				file.imputed.map <- paste(file.imputed.new, ".map", sep="");

				if(file.exists(file.imputed.ped))
					file.seq.plink <- c( file.seq.plink, paste( file.imputed.ped, file.imputed.map, sep=" ") )
				
			}
		}
	}

	merge.txt <- paste("merge-", chr.x, ".txt", sep="");
	write.table(file.seq.plink, file=merge.txt, quote=F, row.names=F, col.names=F);
	
	file.impute2.ped <- paste(file.ped.out, ".impute2", sep="");
	cmd.plink <- paste(path.plink, "--noweb --file ", file.imputed.new1, "--merge-list", merge.txt, "--make-bed", "--out", file.impute2.ped, sep=" ");
	system( cmd.plink );
	
	return(1);
}

merge.chrs<-function(file.bfile.name, chr.groups, file.out.bfile)
{
	merge.txt <- "merge-all-chrs.lst"
	file.imputed.new1 <- paste(file.bfile.name, "-chr", chr.groups[1], sep="")
	
	cat("", merge.txt, append=FALSE);
	for(chr.i in 2:length(chr.groups) )
	{
		file.nobed <- paste(file.bfile.name, "-chr", chr.i, sep="");
		cat( paste( file.nobed, ".bed ", sep=""), merge.txt, append=TRUE);
		cat( paste( file.nobed, ".bim ", sep=""), merge.txt, append=TRUE);
		cat( paste( file.nobed, ".fam ", sep=""), merge.txt, append=TRUE);
		cat("\n",  merge.txt, append=TRUE);
	}
	
	cmd.plink <- paste(path.plink, "--noweb --file ", file.imputed.new1, "--merge-list", merge.txt, "--make-bed", "--out", file.out.bfile, sep=" ");
	system( cmd.plink );
	
	return( file.out.bfile );

}

do.impute2<-function( plink.obj, path.impute2, path.gtool, path.impute2.data, options=list(n.cores=7) )
{
	require("snowfall");

	cat("[ Impute2 ...]\n");

	subDir  <- "impute2/"
	mainDir <- plink.obj$main.path;
	dir.create( file.path(mainDir, subDir), showWarnings = FALSE );

	path.setwd <- file.path(mainDir, subDir);
	
	op.cpu            <- options$n.cores;
	path.plink        <- plink.obj$plink.path
	file.bfile        <- plink.obj$qc2$plink.out.bfile;
	file.bfile.name   <- "plink-imp2";

	chr.groups <- c(1:22);
	cat("Starting parallel computing, snowfall/snow......\n"); 
	snowfall::sfInit(parallel = TRUE, cpus = op.cpu, type = "SOCK")

	snowfall::sfExport("path.impute2.data", "path.plink", "path.impute2", "path.gtool", "file.bfile", "file.bfile.name", "path.setwd" );

	gls.cluster <- snowfall::sfClusterApplyLB( chr.groups, cpu.fun);

	snowfall::sfStop();

	cat("Stopping parallel computing......\n");
	
	file.out.bfile <- merge.chrs( file.bfile.name, chr.groups, paste(file.bfile, "impute2", sep=".") );
	
	plink.obj$genotype$impute <- list( file.out.bfile = file.out.bfile );
	
	return(plink.obj);
}

