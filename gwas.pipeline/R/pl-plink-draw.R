draw_manhattan<-function( res, map.title="", sig.thres, dot.cex, y.max=NA)
{
	#par( xaxs="r",  yaxs="r");  # Extend axis limits by 4%
	#par( tck = -0.008 );
	#par( mgp = c(2,0.2,0) );
	#par( mar=c( 4, 4, 1.5, 1.5)+0.1);

	pvalues <- -log10(res[,3]);
	nrow    <- dim(res)[1];
	log10.max <- round(max(pvalues, na.rm=T))+1;

	if (is.na(y.max))
		y.max <- log10.max;

	if (length(which(pvalues>y.max))>0)
		pvalues[which(pvalues>y.max)] <- y.max;

	plot( 1,1, type="n", xlab="Chromosome", ylab=expression(-log[10](italic(p))),
		  cex.axis=0.7, xlim=c(1, nrow), ylim=c(0,  y.max), xaxt="n", main=map.title );

	p.lab <-  - c( log10( sig.thres) );

	abline( h= c(p.lab), col="gray", lwd=1, lty="dashed");

	#text( x=0, p.lab[1] + 0.1, "p=0.01", cex=0.8, srt=90, adj=c(0.5, -1)); 

	cols <- c( "darkgreen","black",  "orange",  "red", "blue", "purple");
	
	p.cex <- rep(0.5*dot.cex, length(pvalues));
	p.cex[which(pvalues>0.4*y.max)] <- 0.5*dot.cex + (pvalues[which(pvalues>y.max*0.4)]-0.4*y.max)/y.max*dot.cex;
	points( pvalues, pch=20, col=cols[ (res[,1]%%6+1)], cex=p.cex);

	x.off <- 0;
	x <- c();
	x.ps <- c();
	for(i in 1:18)
	{
		x.ps<-c( x.ps, length(which(res[,1]==i) ) );
		x <- c(x, x.ps[i]/2 + x.off);
		x.off <- x.off + x.ps[i];
	}

	x.ps<-c( x.ps, length(which(res[,1]>18 ) )  )
	x <- c(x, x.ps[19]*2/3 + x.off);
	
	axis(side=1, at=x, col="black", labels=c(paste("", 1:18, sep=""),"..."), col.axis="black", col.lab="black", cex.axis=0.4, padj=-0.2 )
}	


draw_linear_man_fig<-function( file.linear.ret, fig.prefix, fig.title)
{
	tb<-read.table( file.linear.ret, header=T)

	pdf(paste(fig.prefix, "-manh-add.pdf", sep=""), width=6, height=3.5)
	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)
	t1.len <-length(which(tb[,5]=="ADD"))
	draw_manhattan( tb[which(tb[,5]=="ADD"),c(1,2,9)], map.title=paste("Additive/", fig.title, sep=""), 0.05/t1.len, 0.7 );
	dev.off();

	pdf(paste(fig.prefix, "-manh-dom.pdf", sep=""), width=6, height=3.5)
	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)
	t2.len <-length(which(tb[,5]=="DOMDEV"))
	draw_manhattan( tb[which(tb[,5]=="DOMDEV"),c(1,2,9)], map.title=paste("Dominant/", fig.title, sep=""), 0.05/t2.len, 0.7 );
	dev.off();

	pdf(paste(fig.prefix, "-manh-geno.pdf", sep=""), width=6, height=3.5)
	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)
	t5.len <-length(which(tb[,5]=="GENO_2DF"))
	draw_manhattan( tb[which(tb[,5]=="GENO_2DF"),c(1,2,9)], map.title=paste("Geno/", fig.title, sep=""), 0.05/t5.len, 0.7 );
	dev.off();
}

draw_linear_QQ_fig<-function( file.linear.ret, fig.prefix, fig.title)
{
	tb<-read.table( file.linear.ret, header=T)

	pdf(paste(fig.prefix, "-qq-add.pdf", sep=""), width=4, height=4)
	par(cex=0.6);
	par(cex.axis=1.2);
	y8<-tb[which(tb[,5]=="ADD"),8];
	y8 <- y8[!is.na(y8)]
	qqnorm(y8)
	lines(seq(-4,4), seq(-4,4), lwd=2, col="green")
	dev.off()


	pdf(paste(fig.prefix, "-qq-dom.pdf", sep=""), width=4, height=4)
	par(cex=0.6);
	par(cex.axis=1.2);
	y8<-tb[which(tb[,5]=="DOMDEV"),8];
	y8 <- y8[!is.na(y8)]
	qqnorm(y8)
	lines(seq(-4,4), seq(-4,4), lwd=2, col="green")
	dev.off()

	pdf(paste(fig.prefix, "-qq-geno.pdf", sep=""), width=4, height=4)
	par(cex=0.6);
	par(cex.axis=1.2);
	y8<-tb[which(tb[,5]=="GENO_2DF"),8];
	y8 <- y8[!is.na(y8)]
	qqplot(y8, rchisq(length(y8), df=2) )
	lines(seq(0,35), seq(0,35), lwd=2, col="green")
	dev.off()
}

draw_logistic_fig2<-function( file.log.ret, fig.prefix, fig.title)
{
 	tb<-read.table( file.log.ret, header=T)

 	tb0<-tb[which(tb$TEST=="ADD"),]
	pdf(paste(fig.prefix, "-man.pdf", sep=""), width=6, height=3.5)
	par(mar=c(3.2, 3.5,  2.5, 1))
	par(mgp=c(1.4, 0.3, 0), tck=-0.03)
	draw_manhattan( tb0[,c(1,2,9)], map.title=fig.title, 0.05/dim(tb0)[1], 0.7 );
	dev.off();
	
	pdf(paste(fig.prefix, "-qq.pdf", sep=""), width=4, height=4)
	par(cex=0.6);
	par(cex.axis=1.2);

	y8<-tb[which(tb[,5]=="ADD"),8];
	y8 <- y8[!is.na(y8)]
	qqnorm(y8)
	
	lines(seq(-4,4), seq(-4,4), lwd=2, col="green")
	dev.off()	
}

show_plink_topn<-function( file.linear, tag="ADD", top=10 )
{
	tb<-read.table(file.linear, header=T)
	tb.add<-tb[which((tb$TEST==tag) & !is.na(tb[,9])),]
 
 	tb.add.st <- sort.int( tb.add[,9], decreasing=F, index.return=T);
 	
 	tb.add.top <- tb.add[ tb.add.st$ix[1:top], ];
 	
 	ret<-c();
 	for(i in 1:top)
 	{
	 	t <- which( tb[,1] == tb.add.top[i, 1] & tb[,2] == tb.add.top[i, 2] );
	 	ret <- rbind(ret, tb[t,]);
 	}
 	
	return(list(full=ret, top=tb.add.top));
}

show_plink_sign<-function( file.linear, tag="ADD", sig=0.05 )
{
	tb<-read.table(file.linear, header=T)
	tb.add<-tb[which((tb$TEST==tag) & !is.na(tb[,9])),]
 
 	tb.add.st <- sort.int( tb.add[,9], decreasing=F, index.return=T);
 	tb.add.st<- tb.add[ tb.add.st$ix, ];
 
 	sign <- which( tb.add.st[,9] <= sig/dim(tb.add.st)[1] )
 	if (length(sign)>0)
 	{	
		tb.add.sig <- tb.add.st[sign,]
		
		ret<-c();
		for(i in 1:length(sign))
		{
			t <- which( tb[,1] == tb.add.sig[i, 1] & tb[,2] == tb.add.sig[i, 2] );
			ret <- rbind(ret, tb[t,]);
		}
		return(list(full=ret, sig=tb.add.sig));
	}
	else
		return(list(full=NULL, sig=NULL));

}
