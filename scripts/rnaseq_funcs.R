
# calculate how many replicates in each group are above 1 cpm expression
# you might have to alter that threshold depending on sequencing depth and other considerations, for a library size of 12 million reads, that means around 12 counts
# see suggestions in Filtering to remove low counts in Chen et al., F1000Research 2016, 5:1438, http://dx.doi.org/10.12688/f1000research.8987.2

importFromEdgeR <- function(obj, res, org_name, useF= T, ...)
{
  #import from Edger QLF DGEList object and a TopTags results table into a summarized experiment object for enrichment browser
  # allows to add leading edge and stuff
  # uses the F statistic by default for ranking
  stopifnot(is(obj, "DGEList") && is(res, "TopTags"))
  res <- res$table
  padj <- grep("^F[WED]{1,2}R$", colnames(res), value = TRUE)
  rnames <- c("logFC", "F", "PValue", padj)
  stopifnot(all(rnames %in% colnames(res)))
  
  se <- SummarizedExperiment(assays = list(counts = obj$counts),
                             colData = obj$samples)
  
  assay(se, "norm") <- limma::voom(obj)$E
  names(assays(se))[1] <- "raw" 
  
  # (1) rowData
  rrnames <- c("FC.COL", "PVAL.COL", "ADJP.COL")
  rrnames <- vapply(rrnames, configEBrowser, character(1), USE.NAMES = FALSE)
  ind <- match(rnames, colnames(res))
  colnames(res)[ind] <- c(rrnames[1], "edgeR.STAT", rrnames[2:3])
  if (useF){ # sing of FC * sqrt(F)
    res$edgeR.STAT <- sign(res[, rrnames[1]]) * sqrt(res[, "edgeR.STAT"])
  } else { # -1*sign(FC)*log10(P)
    res$edgeR.STAT <- -1 * sign(res[, rrnames[1]]) * log10(res[, rrnames[2]])
  }
  
  rowData(se) <- res
  
  # (2) colData
  vars <- names(attr(obj$design, "contrasts"))
  design <- obj$design[, -1, drop = FALSE]
  gcol <- configEBrowser("GRP.COL")
  grp.col <- design[,ncol(design)]
  ngrps <- length(unique(grp.col))
  if(length(vars) > 2 || ngrps > 2) 
    stop("Supported experimental designs include binary",
         " group comparisons with an optional blocking",
         " variable for paired samples / sample batches")
  # (2i) group
  se[[gcol]] <- unname(grp.col)
  # (3) metadata 
  metadata(se)$dataType <- "rseq"
  metadata(se)$annotation <- org_name
  metadata(se)$common.dispersion <- obj$common.dispersion
  return(se)
}



check_expression_threshold <- function(exobj,n=NULL,min_cpm=1.0,min_count=NULL,by_group=TRUE,groups=NULL){
  # this is similar to filterByExpr in the edgeR package
  # check ?filterByExpr and maybe use that if you feel more comfortable
  # returns a logical vector with length of genes with TRUE if expression in n replicates of a group is higher than min_cpm and min_count
  # n defaults to int(3/4) of the smallest group size. if by_group=FALSE no grouping is used, groups can be a factor vector indicating group membership, if a DGEList object is used and group is NULL, groups=dgeobject$samples$group
  # if min_count or min_cpm are set to NULL, either is ignored
  if (by_group){
    if (is.null(groups)) {
      if ("samples" %in% names(exobj) && ! is.null(exobj$samples$group)) {
        groups <- exobj$samples$group
      }
      else {
        groups=as.factor(rep(1,ncol(exobj)))
      }
    }
    else {
      groups = as.factor(groups)
    }
  } 
  else {
    groups=as.factor(rep(1,ncol(exobj)))
  }
  if (is.null(n)){
    n <- as.integer(min(tabulate(groups))*3/4)
  }
  keep=rep(TRUE,nrow(exobj))
  if (! is.null(min_cpm)) {
    if (by_group){
      suf_by_group=tapply(colnames(exobj),groups, function (x) rowSums(cpm(exobj,normalized.lib.sizes = FALSE)[,x] >= min_cpm ))
      suf_by_group=sapply(names(suf_by_group), function (x) suf_by_group[[x]]) 
      keep= keep & (rowSums(suf_by_group >= n) > 0)
    } else {
      suf_by_group = rowSums(cpm(exobj,normalized.lib.sizes = FALSE) >= min_cpm )
      keep = keep & (suf_by_group >= n)
    }
    
  }
  if (! is.null(min_count)) {
    if (by_group){
    suf_by_group=tapply(colnames(exobj),groups, function (x) rowSums(as.matrix(exobj[,x]) >= min_count ))
    suf_by_group=sapply(names(suf_by_group), function (x) suf_by_group[[x]])
    keep = keep & (rowSums(suf_by_group >= n) > 0)
    } else {
      suf_by_group = rowSums(as.matrix(exobj) >= min_count )
      keep = keep & (suf_by_group >= n)
    }
    
  }
  return(keep)
}

#set up colors for samples (similar colors for replicates, can be used later for plots)
get_group_sample_cols <- function(groups=mmdge$samples$group,replicates=mmdge$samples$replicates,palette="Set3"){
  #takes two factor vectors for group and replicates and returns a vector for group and sample colors according to group and replicate vectors as a list
  # if the first element is a DGElist object, it extracts groups and replicates from that
  if(class(groups) == "DGEList") {
    replicates = groups$samples$replicates
    groups = groups$samples$group
  }
  if (length(levels(groups)) <=12){
    group_cols=brewer.pal(n=length(levels(groups)),name=palette)
  } else {
    group_cols=colorRampPalette(brewer.pal(n=12,name=palette))(length(levels(groups)))
  }
  
  names(group_cols)=levels(groups)
  sample_cols=group_cols[groups]
  # change the colors a tiny bit for each replicate
  sample_cols=sapply(1:length(sample_cols), function (x) 
    adjustcolor(sample_cols[x],transform=diag(c(1, 1, 1, 1) - (replicates[x]-mean(replicates))*c(0.20,0.20,0.20,0) ) ))
  names(sample_cols)=rownames(replicates)
  return(list("group_cols"=group_cols,"sample_cols"=sample_cols))
}

plot_MA_matrix <- function(mmdge, averageover=c("all","group")){
  # calculates average A and M values per all samples, group and ind (if defined)
  # MA plots between all replicates of a group
  # each point corresponds to a gene with the log expression averaged over a group of samples (A) on the x axis and the log2 fold expression change of a sample against the average as the y axis (M)
  # this plot can reveal systematic biases (especially comparing a sample to other replicates) or changes in expression of large number of genes (against the overall average)
  # need to first create matrix of A (average log expression) and M (log fold expression difference) values
  # all is plotted in one big matrix to plot, use ggplot with hexagonal binning to reduce the number of points
  AMm=get_A_M_values(mmdge,averageover)
  # plot MA of each sample against overall mean of all other samples
  # the red trend line indicates the mean M value and should be close to y=0
  if ("all" %in% averageover) {
    p = ggplot(AMm,aes(x=A,y=M)) + stat_binhex() + geom_smooth(color="red",se=FALSE) + 
      facet_grid(rows=vars(group),cols=vars(replicates)) + 
      labs(y = "log-ratio (sample against average)", x = "average log CPM (all samples)" )
    plot(p)
  }
  # plot MA of each replicate against mean of its sample group, these should be more similar and can reveal problems between replicates
  if ("group" %in% averageover) {
    p =  ggplot(AMm,aes(x=Ag,y=Mg)) + stat_binhex() + geom_smooth(color="red",se=FALSE) + 
      facet_grid(rows=vars(group),cols=vars(replicates)) + 
      labs(y = "log-ratio (sample against group average)", x = "average log CPM (replicate samples in group)" )  
    plot(p)
  }
}

plot_filtered_genes <- function(mmdge,minval=c(0.5,1,2.5,5,7.5,10,15),n=3,use.cpm=TRUE,palette="Blues",
                                per_lib=FALSE,gen_above_threshold=FALSE,add_kept_gene_frac=FALSE) {
  # n: number of samples that have to have expression above threshold
  # per_lib: only look at genes filtered per library
  
  grp_samp_cols=get_group_sample_cols(mmdge)
  counts=if(use.cpm) cpm(mmdge) else mmdge$counts
  if (per_lib){
    fracs= sapply(minval,function (x) colSums(counts < x) )/nrow(counts)
  } else {
    fracs = sapply(minval,function (x) { keep=rowSums(counts >= x) > n;  y = if(gen_above_threshold) x else 0;  colSums(counts[keep,] <= y)/sum(keep) }  )
    if (add_kept_gene_frac) {
      fracs=rbind(fracs,sapply(minval,function (x) { keep=rowSums(counts >= x) > n;  sum(keep)/nrow(counts) }) )
      rownames(fracs)[nrow(fracs)] = "kept overall"
    }
  }
  
  colnames(fracs)=minval
  old.par <-par(no.readonly = T)
  par(mar=c(10,3,4,5))
  leg_title=if(use.cpm) "min. cpm" else "min. count"
  main_tit_part=if(gen_above_threshold) "Frac. o. genes below threshold after filt. above threshold in at least" else "Frac. o. unexpressed genes after filt. above threshold in at least"
  main_title=if(per_lib) "Frac. o. genes filt. with threshold per lib." else paste(main_tit_part,as.character(n),"samples") 
  barplot(t(fracs),beside=TRUE,legend.text=TRUE,args.legend=list(title=leg_title,x ="right",bty="n",inset=c(-0.10,0),xpd = TRUE),
          col=brewer.pal(n=ncol(fracs),name=palette),cex.names=0.5,las=2,main=main_title)
  par(old.par)
}


plotBCV_groups <- function(dgelist){
  # plot the biological coefficient of variation per group
  groups <- levels(dgelist$samples$group)
  if (length(groups) == 3) {
    par(mfrow = c(2,2))
  } else {
    par(mfrow = n2mfrow(length(groups)))
  }
  for(i in groups) {
    x <- estimateDisp(dgelist[,dgelist$samples$group == i])
    plotBCV(x,main=i)
  }
  # fill up page so not to leave plots empty
  while ( ! par()$page) { plot.new() }
  par(mfrow=c(1,1))
}


plot_overall_counts <- function(mmdge,normalized=TRUE,prior.count=0.5,which_plots="both") {
  ## prior.count : added for getting logcpm not to be bothered by 0 vals
  ## produces two figures, which_plots: first, second, both
  grp_samp_cols=get_group_sample_cols(group=mmdge)
  pcpm=prior.count/10^6
  norm=if(normalized) "TMM normalized" else ("library size normalized")
  if(which_plots == "both" | which_plots == "first") {
    # boxplot of log(CPM) with CPM: counts per million reads in library
    par(mfrow=c(2,2))
    names=paste(mmdge$samples$group,mmdge$samples$replicates,sep="_")
    boxplot(log10(cpm(mmdge,normalized.lib.sizes = normalized)+pcpm)[1:100,],col=c(grp_samp_cols$sample_cols),las=2,cex.axis=0.75, ylab="log10(CPM)",names=names,main=paste(norm,"CPM values of genes"))
    
    # plot library sizes (reads mapped to filtered genes)
    barplot(mmdge$samples$lib.size/10^6,col=grp_samp_cols$sample_cols,main="Library sizes",names.arg=names, las=2, cex.names =0.75, ylab="million reads")
    ## plot fraction of zero counts per sample
    barplot(apply(mmdge$counts,2,function (x) length(x[ x == 0])/length(x)),col=grp_samp_cols$sample_cols,
            names.arg=names, las=2, cex.names =0.75, ylab="fraction 0 count genes",main="genes with 0 counts")
    plot.new()
    par(xpd=TRUE)
    legend("topleft",legend=paste(names,rownames(mmdge$samples),sep=" : "),cex=0.75,fill=grp_samp_cols$sample_cols)
    par(xpd=FALSE)
    par(mfrow=c(1,1))
    
  }
  if(which_plots == "both" | which_plots == "second") {
  # plot density of log cpm
  # change data frame to be used with ggplot
  lcpm_m=melt(log10(cpm(mmdge,normalized.lib.sizes = normalized)+pcpm),value.name = "log10CPM")
  colnames(lcpm_m)[1:2]=c("gene","sample")
  # density plot using ggplot2
  p = ggplot(lcpm_m, aes(x=log10CPM,color=sample)) + stat_density(geom="line",position = "identity") + 
    scale_color_manual(values = grp_samp_cols$sample_cols) + guides(col=guide_legend(ncol=1)) + labs(x="log10(CPM)", title=paste(norm,"CPM values of genes") ) +
    theme(legend.title = element_text(), legend.text = element_text(size=7.0))
  plot(p)
  }
}


get_A_M_values <- function(mmdge,avagainst=c("all","group")){
  # this function returns M and A values averaged over all samples, by group, and by ind
  # gives table with columns:
  # A, M : averaged over all samples
  # Ag,Mg : averaged by group
  # averages over all samples
  lcpm_n=cpm(mmdge,log=TRUE)
  As=rowMeans(lcpm_n) 
  # Ms over all sample averages
  # melt data frame for use in ggplot2
  Mm=melt(lcpm_n)
  colnames(Mm)=c("gene","sample","lCPM")
  AMm=cbind(Mm,As[Mm$gene])
  colnames(AMm)=c("gene","sample","lCPM","A")
  AMm$M=AMm$lCPM-AMm$A
  AMm=cbind(AMm,mmdge$samples[Mm$sample,c("group","replicates")])
  # calculate A values per group (not overall)
  # this might be more appropriate, to see whether replicates deviate from each other
  if ("group" %in% avagainst){
    Ag=tapply(AMm$lCPM,AMm[,c("gene","group")],mean,simplify = TRUE)
    AMm$Ag=sapply(seq(1,nrow(AMm)), function (x) Ag[AMm$gene[x],AMm$group[x]] ) 
    AMm$Mg=AMm$lCPM-AMm$Ag
  }
  return(AMm)
}





niceMDS <- function(exobj,annos=NULL,color_by="group", shape_by="replicates",
                    lables=NULL,ndims=5,xdim="Dim 1",ydim="Dim 2", title=NULL,top=500,gene.selection="pairwise", filename = NULL){
  ### takes an object which plotMDS from the limma package can take and an annotations sheet and plots the dimensions given by xdim and ydim ("Dim n")
  ### for coloring, shapes and labelling give the name of the factor column in the annotation sheet
  ### returns a data frame with the MDS dimensions and the annotation sheet concatenated for easy plotting with ggplot2
  # plot Multidimensional scaling (MDS) plots for samples using plotMDS from the limma package
  # Like Principle component analysis, MDS allows to represent high dimensional data in lower dimensions
  # from Chen et al., F1000Research 2016, 5:1438, http://dx.doi.org/10.12688/f1000research.8987.2 :
  # In the MDS plot, the distance between each pair of samples can be interpreted as the leading log-fold change between the samples for the genes that best distinguish that pair of samples. By default, leading fold-change is defined as the root-mean-square of the largest 500 log2-fold changes between that pair of samples.
  # in this plots, replicates should be closer to each other, than unrelated samples
  # you can color/label samples by different annotations, to see whether other influences lead to clustering of samples 
  # the MDS dimensions/components are ranked by the amount of difference the explain between samples. you can also look at higher change dimensions, apart from the first 2
  # the direct plots produced are alright, but do not look as nice and are a bit harder to configure, so we just extract the result data and plot it in ggplot
  mdsdata=plotMDS(exobj, top=top, gene.selection=gene.selection, ndim=ndims,plot=FALSE)
  if (is.null(annos) & is(exobj,"DGEList")) {
    annos = exobj$samples
  }
  mds_gg=cbind(mdsdata$cmdscale.out,annos[rownames(mdsdata$cmdscale.out),])
  mds_gg[c(color_by,shape_by)] <- lapply(mds_gg[c(color_by,shape_by)] , factor)
  #mds_gg$replicates=as.factor(mds_gg$replicates)
  dim_names = paste("Dim",1:ndims)
  colnames(mds_gg)[1:ndims]=dim_names
  #xdim=dim_names[1]; ydim=dim_names[2];
  gg = ggplot(mds_gg,aes(x = !!ensym(xdim), y = !!ensym(ydim),colour=!!ensym(color_by),shape=!!ensym(shape_by))) + geom_point() +
    labs(x = paste("leading logFC",xdim), y = paste("leading logFC",ydim) ) + ggtitle(title) +
    geom_hline(color="grey", linetype="dashed", yintercept = 0) + geom_vline(color="grey", linetype="dashed", xintercept=0)
  if (!is.null(lables)) {
    gg = gg + geom_label_repel(aes(label = mds_gg[,lables]))
  }
  if (is.null(filename)){
    plot(gg)
  } else {
    ggsave(filename, plot=gg)
  }
  return(mds_gg)
}

## helper function to only get the topvar rows of a data frame with highest variation
get_top_var<- function(exprs_mat,topvar=1000){
  ## only returns the topvar rows with the highest variation from a dataframe/matrix exprs_mat
  row_vars <- apply(exprs_mat,1,var)
  if ( ! is.finite(topvar) || topvar < 1 || topvar > length(row_vars)) {
    topvar = length(row_vars)
  }
  topvar<-as.integer(topvar)
  return(exprs_mat[order(row_vars,decreasing = TRUE)[1:topvar],])
}

get_topvar_by_group <- function(exprs_mat,groups,topvar=100) {
  # performs variance estimation per level in factor vector group for all rows in exprs_mat
  # returns the topvar genes for each group with the highest within group variance
  row_vars <- as.data.frame(sapply(levels(groups), function(x) rowVars(exprs_mat[],cols=which(groups == x))))
  rownames(row_vars) <- rownames(exprs_mat)
  if ( ! is.finite(topvar) || topvar < 1 || topvar > nrow(row_vars)) {
    topvar = nrow(row_vars)
  }
  topvar<-as.integer(topvar)
  top_row_var <- lapply(colnames(row_vars),function(x) { y<-row_vars[order(row_vars[,x],decreasing = TRUE)[1:topvar],x,drop=FALSE]; return(setNames(y[[1]],rownames(y))) } )
  names(top_row_var) <- colnames(row_vars)
  top_row_var[["COMBINED"]]=unique(as.vector(sapply(top_row_var,names)))
  names(top_row_var[["COMBINED"]])=top_row_var[["COMBINED"]]
  return(top_row_var)
  
}


# Performing Principle component analysis (PCA) and plotting Principle Components for the n genes varying most over all samples
# Here we only take the 1000 most varying genes for doing a PCA, similar to what DESeq2 does with their plotPCA function
nicePCA <- function(exobj=cpm(mmdge,log=TRUE),annos=mmdge$samples,topvar=1000,color_by="group", shape_by="replicates",lables=NULL,xdim="PC1",ydim="PC2",center=TRUE,scale=FALSE,title=NULL,plot_skree=TRUE, filename = NULL ){
  ### takes a numerical matrix with expression values of genes x sample (eg. cpm(mmdge,log=TRUE)) and an annotations sheet, performes centering over genes and PCA on the transposed matrix and plots a skree plot and the dimensions given by xdim and ydim ("Dim n")
  ### the PCA is performed on the topvar most varying genes (default 5000), if topvar = 0,NA,NULL or similar all genes are considered
  ### for coloring, shapes and labelling give the name of the factor column in the annotation sheet
  ### returns a prcomp output with a data frame added under pca$tab with the principle components and the annotation sheet concatenated for easy plotting with ggplot2
  ### if given a prcomp object, only does the plotting
  if (! class(exobj) == "prcomp") {
    if (is(exobj,"DGEList")) {
      annos=exobj$samples;
      exobj=cpm(exobj,log=TRUE);
    }
    exobj_top_var=get_top_var(exobj,topvar)
    # perform PCA (transposed exobj matrix to do it in sample space and to center over each gene (columns), not samples)
    pca <- prcomp(t(exobj_top_var), center = center, scale. = scale) 
    #skree plot (Variance explained by the different components)
    if (plot_skree) {
      barplot(pca$sdev^2/sum(pca$sdev^2),ylab="% variance explained",names.arg = colnames(pca$x),las=2)
    }
  } else {
    pca <- exobj
  }
  if (is.data.frame(annos)){
    pca$tab=cbind(pca$x,annos[rownames(pca$x),])
    pca$tab[c(color_by,shape_by)] <- lapply(pca$tab[c(color_by,shape_by)] , factor)
  }
  # get nice strings for variance explained for PCs
  var_expl=paste0("(",round(100*(pca$sdev^2/sum(pca$sdev^2)),digits=1),"%)")
  names(var_expl)=colnames(pca$x)
  # do actual plotting
  gg = ggplot(pca$tab,aes(x = !!ensym(xdim), y = !!ensym(ydim),colour=!!ensym(color_by),shape=!!ensym(shape_by))) + 
    geom_point(size=3) + labs(x = paste(xdim,var_expl[xdim]), y = paste(ydim,var_expl[ydim]) ) + ggtitle(title) +
    geom_hline(color="grey", linetype="dashed", yintercept = 0) + geom_vline(color="grey", linetype="dashed", xintercept=0)
  if (!is.null(lables)) {
    gg = gg + geom_label_repel(aes(label = pca$tab[,lables]))
  }
  if(is.null(filename)){
    plot(gg)
  } else {
    ggsave(filename,gg)
  }
  
  return(pca)
}


get_r2_corPC <- function(pca,annos,dims=c(1:4)){
  res=list()
  for (x in names(annos)){
    res[[x]]=data.frame(p_lm=rep(NA,length(dims)), adj.r2=rep(NA,length(dims)),p_krusk=rep(NA,length(dims)))
    rownames(res[[x]])=as.character(dims)
    for(i in dims){
      res[[x]][as.character(i),"adj.r2"] <- summary(lm(pca$x[,i] ~ annos[,x]))$adj.r.squared
      #f <- a$fstatistic
      #res[[x]][as.character(i),"p_lm"] <- pf(f[1],f[2],f[3],lower.tail=F)
      res[[x]][as.character(i),"p_lm"] <- anova(lm(pca$x[,i] ~ annos[,x]))$`Pr(>F)`[1]
      res[[x]][as.character(i),"p_krusk"] <-kruskal.test(pca$x[,i]~annos[,x])$p.value
    }
  }
  return(res)
}

perform_gPCA <- function(dgeobj,batch="group",topvar=1000){
  require("gPCA")
  gPCA.out <- gPCA.batchdetect(t(get_top_var(cpm(dgeobj,log=TRUE),topvar)),batch=as.numeric(dgeobj$samples[,batch]))
  old.par=par()
  gDist(gPCA.out)
  PCplot(gPCA.out,ug="guided",type="1v2")
  PCplot(gPCA.out,ug="unguided",type="1v2")
  #PCplot(gPCA.out,ug="guided",type="comp",npcs=3)
  par(old.par)
  CumulativeVarPlot(gPCA.out,ug="unguided",col="blue")
  CumulativeVarPlot(gPCA.out,ug="guided",col="blue")
  return(gPCA.out)
}

vulcano.plot <- function(fit2, contrast, method="BH", pch=16,cex=0.45, cols=c("red","green","orange"), xlab="log2FC", ylab="-log10P", ...) {
  plot(fit2[["coefficients"]][,contrast],-log10(fit2[["p.value"]][,contrast]), pch=pch,cex=cex, xlab=xlab, ylab=ylab, ...)
  sig= p.adjust(fit2[["p.value"]][,contrast],method=method) <= 0.05
  lfgt2= abs(fit2[["coefficients"]][,contrast]) > 1.0
  points(fit2[["coefficients"]][sig & ! lfgt2,contrast],-log10(fit2[["p.value"]][sig & ! lfgt2,contrast]), pch=pch,cex=cex,col=cols[1])
  points(fit2[["coefficients"]][lfgt2 & ! sig,contrast],-log10(fit2[["p.value"]][lfgt2 & ! sig,contrast]), pch=pch,cex=cex,col=cols[2])
  points(fit2[["coefficients"]][sig & lfgt2,contrast],-log10(fit2[["p.value"]][sig & lfgt2,contrast]), pch=pch,cex=cex,col=cols[3])
}

plotCompareP_alt <- function (p1, p2, vpDonor = NULL, dupcorvalue = 0, fraction = 0.2, psize=2, 
                              xlabel = "duplicateCorrelation", ylabel = "dream", adjust = FALSE, title = NULL, base_size=10) 
{
  xlabel = bquote(.(xlabel) ~ (-log[10] ~ p))
  ylabel = bquote(.(ylabel) ~ (-log[10] ~ p))
  df2 = if (adjust) data.frame(p1 = -log10(p.adjust(p1,method="BH")), p2 = -log10(p.adjust(p2,method="BH"))) else data.frame(p1 = -log10(p1), p2 = -log10(p2))
  
  df2$vpdonor=0
  if (! is.null(vpDonor)){
    if (length(unique(c(length(p1), length(p2), length(vpDonor)))) != 
        1) {
      stop("p1, p2 and vpDonor must have the same number of entries")
    }
    if (length(dupcorvalue) != 1) {
      stop("dupcorvalue must be a scalar")
    }
    df2$vpDonor = vpDonor 
    df2$delta = vpDonor - dupcorvalue
    N = nrow(df2)
    c1 = sort(df2$delta)[N * fraction]
    c2 = sort(df2$delta)[length(df2$delta) - N * fraction]
    l1 = lm(p2 ~ p1, df2[df2$delta >= c2, ])
    l2 = lm(p2 ~ p1, df2[df2$delta <= c1, ])
    df_line = data.frame(rbind(coef(l1), coef(l2)))
    colnames(df_line) = c("a", "b")
    df_line$type = c("darkred", "navy")
  }
  lim = c(0, max(max(df2$p1[df2$p1 < Inf]), max(df2$p2[df2$p2 < Inf])))
  title = if (adjust) paste(title,"BH adjusted P") else title
  gp = ggplot(df2, aes(p1, p2, color = vpDonor)) + geom_point(size = psize) + geom_abline() + 
    geom_hline(yintercept = -log10(0.05), colour="black",linetype="dotted") + geom_vline(xintercept = -log10(0.05), colour="black",linetype="dotted") + 
    theme_bw(base_size = base_size) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + xlim(lim) + ylim(lim) +
    labs(x=xlabel, y=ylabel, title=title) 
  if (! is.null(vpDonor)) {
    gp = gp + geom_abline(intercept = df_line$a, slope = df_line$b, color = df_line$type, linetype = 2) + 
    scale_color_gradientn(name = "Ind", colours = c("blue", "green", "red"), values = rescale(c(0, dupcorvalue, 1)), 
                          guide = "colorbar", limits = c(0, 1))
  }
  gp
}

glimmaMA.DESeqResults  <- function(x,dds,counts=DESeq2::counts(dds),groups=colData(dds)$group, status=NULL, anno=NULL,
                                   display.columns = NULL, status.cols=c("dodgerblue", "silver", "firebrick"), sample.cols=NULL, transform.counts = c("logcpm", "cpm", "rpkm", "none"),
                                   main="MA Plot", xlab="logCPM", ylab="logFC", html=NULL, width = 920, height = 920, ...)
{
  transform.counts <- match.arg(transform.counts)
  res.df <- as.data.frame(x)
  # filter out genes that have missing data
  complete_genes <- complete.cases(res.df)
  res.df <- res.df[complete_genes, ]
  x <- x[complete_genes, ]
  
  total_genes <- length(complete_genes)
  filtered_genes <- sum(!complete_genes)
  message(filtered_genes, " of ", total_genes, " genes were filtered out in DESeq2 tests")
  
  # extract status if it is not given
  if (is.null(status))
  {
    status <- ifelse(
      res.df$padj < 0.05,
      ifelse(res.df$log2FoldChange < 0, -1, 1),
      0
    )
  }
  else
  {
    if (length(status)!=length(complete_genes)) stop("Status vector
      must have the same number of genes as the main arguments.")
    status <- status[complete_genes]
  }
  
  # create initial table with logCPM and logFC features
  table <- data.frame(signif(log(res.df$baseMean + 0.5), digits=4),
                      signif(res.df$log2FoldChange, digits=4))
  colnames(table) <- c(xlab, ylab)
  table <- cbind(table, PValue=signif(res.df$pvalue, digits=4), AdjPValue=signif(res.df$padj, digits=4))
  table <- cbind(gene=rownames(x), table)
  xData <- Glimma:::buildXYData(table, status, main, display.columns, anno[complete_genes, ], counts, xlab, ylab, status.cols, sample.cols, groups, transform.counts)
  return(Glimma:::glimmaXYWidget(xData, width, height, html))
}




glimmaVolcano.DeSeqResults  <- function(x,dds,counts=DESeq2::counts(dds),groups=colData(dds)$group, status=NULL, anno=NULL,
                                   display.columns = NULL, status.cols=c("dodgerblue", "silver", "firebrick"), sample.cols=NULL, 
                                   transform.counts = c("logcpm", "cpm", "rpkm", "none"),
                                   main="Volcano plot", xlab="logFC", ylab="negLog10P", html=NULL, width = 920, height = 920, ...)
{
  transform.counts <- match.arg(transform.counts)
  res.df <- as.data.frame(x)
  # filter out genes that have missing data
  complete_genes <- complete.cases(res.df)
  res.df <- res.df[complete_genes, ]
  x <- x[complete_genes, ]
  
  total_genes <- length(complete_genes)
  filtered_genes <- sum(!complete_genes)
  message(filtered_genes, " of ", total_genes, " genes were filtered out in DESeq2 tests")
  
  # extract status if it is not given
  if (is.null(status))
  {
    status <- ifelse(
      res.df$padj < 0.05,
      ifelse(res.df$log2FoldChange < 0, -1, 1),
      0
    )
  }
  else
  {
    if (length(status)!=length(complete_genes)) stop("Status vector
      must have the same number of genes as the main arguments.")
    status <- status[complete_genes]
  }
  
  # create initial table with logFC and -log10(pvalue) features
  table <- data.frame(signif(res.df$log2FoldChange, digits=4),
                      signif(-log10(res.df$pvalue), digits=4))
  colnames(table) <- c(xlab, ylab)
  table <- cbind(table, PValue=signif(res.df$pvalue, digits=4), AdjPValue=signif(res.df$padj, digits=4))
  table <- cbind(gene=rownames(x), table)
  xData <- Glimma:::buildXYData(table, status, main, display.columns, anno[complete_genes, ], counts, xlab, ylab, status.cols, sample.cols, groups, transform.counts)
  return(Glimma:::glimmaXYWidget(xData, width, height, html))
}

# for enrichment browser:
# use fsge multilevel as suggested, hopefully alright
use.fgse <-
  function(se ,
           gs,
           alpha = 0.05,
           #perm = 10000,
           full_res = FALSE,
           stat_name = "DESeq2.STAT",
           minSize = 10,
           maxSize = 500) {
    #eset = as(se,"ExpressionSet")
    genes = rownames(rowData(se))
    if ( ! stat_name %in% colnames(rowData(se)) ){
      stat_name = grep("STAT$", colnames(rowData(se)),  value = TRUE)[1]
    }
    # if( stat_name == "edgeR.STAT") {
    #   ## from https://support.bioconductor.org/p/120780/ should also do z value, but cannot easily get the df
    #   t.values <- sign(rowData(se)[, "FC"]) * sqrt(rowData(se)[, "edgeR.STAT"])
    #   if (grep("^df.total$", colnames(rowData(se)))){
    #     t.values <- zscoreT(t.values, df=rowData(se)[, "df.total"]) 
    # } else {
      t.values <- rowData(se)[, stat_name]
#      }
#    }
    names(t.values) = genes
    fgseaRes <-
      fgsea(gs,
            t.values,
            #nperm = perm,
            minSize = minSize,
            maxSize = maxSize)
    fgseaRes[is.na(fgseaRes$pval), c("pval","padj")] = 1.0
    #gStestsBH=p.adjust(gStests,method = "BH")
    if (full_res) {
      return(fgseaRes)
    } else {
      return(setNames(fgseaRes$pval, fgseaRes$pathway))
    }
  }

get.nsig <- function(res) {
  p.values = rowData(res$se)$ADJ.PVAL
  sig.gn = rownames(rowData(res$se))[p.values < 0.05]
  gns = res$res.tbl$GENE.SET
  NR.SIG.GENES = sapply(gns, function (x)
    sum(sig.gn  %in% res$gs[[x]]))
  return(NR.SIG.GENES)
}

get.nsig_up_down <- function(res,stat_name = "DESeq2.STAT") {
  #eset = as(res$se,"ExpressionSet")
  p.values = rowData(res$se)$ADJ.PVAL
  if ( ! stat_name %in% colnames(rowData(res$se)) ){
    stat_name = grep("STAT$", colnames(rowData(res$se)),  value = TRUE)[1]
  }
  stat = rowData(res$se)[, stat_name]
  sig.gn.up = rownames(rowData(res$se))[p.values < 0.05 & stat > 0]
  sig.gn.dwn = rownames(rowData(res$se))[p.values < 0.05 &
                                           stat < 0]
  gns = res$res.tbl$GENE.SET
  NR.SIG.GENES.UP = sapply(gns, function (x)
    sum(sig.gn.up  %in% res$gs[[x]]))
  NR.SIG.GENES.DOWN = sapply(gns, function (x)
    sum(sig.gn.dwn  %in% res$gs[[x]]))
  return(data.frame(N.SIG.UP = NR.SIG.GENES.UP, N.SIG.DOWN = NR.SIG.GENES.DOWN))
}

get.ngenes <- function(res) {
  #eset = as(res$se,"ExpressionSet")
  NR.GENES = sapply(res$res.tbl$GENE.SET, function (x)
    sum(rownames(rowData(res$se))  %in% res$gs[[x]]))
  return(NR.GENES)
}

geneReportLE <- function(s, gt, out.dir, sall, se, stat_name = "DESeq2.STAT") {
  # create enrichment plot
  if ( ! stat_name %in% colnames(rowData(se)) ){
    stat_name = grep("STAT$", colnames(rowData(se)),  value = TRUE)[1]
  }
  p1 = plotEnrichment(unlist(GSEABase::geneIds(sall)),
                      setNames(rowData(se)[, stat_name], rownames(rowData(se))) ) + labs(title = GSEABase::setName(sall))
  ggsave(
    paste0(out.dir, "/", GSEABase::setName(s), "_enrichmentplot.png"),
    plot = p1,
    width = 15,
    height = 8,
    units = "cm"
  )
  # (0) extract gene information
  gt <- gt[GSEABase::geneIds(s), ]
  gt <- EnrichmentBrowser:::.sortGeneTable(gt)
  
  # (1) html table
  htmlRep <- ReportingTools::HTMLReport(
    basePath = dirname(out.dir),
    reportDirectory = basename(out.dir),
    shortName = GSEABase::setName(s),
    title = paste(GSEABase::setName(s), "Gene Report", sep = ": ")
  )
  ReportingTools::publish(gt, htmlRep, )#.modifyDF = list(EnrichmentBrowser:::.ncbiGeneLink))
  # (2) flat file
  fname <- paste0(GSEABase::setName(s), "_genes.txt")
  ofname <- file.path(out.dir, fname)
  if (!file.exists(ofname))
    write.table(
      gt,
      file = ofname,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  dlink <- ReportingTools::Link("Download .txt", fname)
  # (3) add enrichment plot
  ReportingTools::publish(dlink, htmlRep)
  himg <-
    hwriter::hwriteImage(paste0(GSEABase::setName(s), "_enrichmentplot.png"), link =
                           "Enrichment Plot")
  ReportingTools::publish(hwriter::hwrite(himg, br = TRUE), htmlRep)
  
  rep <- ReportingTools::finish(htmlRep)
  return(rep)
}

#eaBrowse_alt(fgse_reactome, html.only = TRUE,out.dir = paste0(getwd(),"/EBrowser2/REACTOME_", i, "/"), report.name =  paste("REACTOME",i, sep="_"))
eaBrowse_alt <-function (res,
            nr.show = -1,
            graph.view = NULL,
            html.only = FALSE,
            out.dir = NULL,
            report.name = NULL) {
    isAvailable("ReportingTools", type = "software")
    method <- ifelse(is(res$method, "character"), res$method,
                     NA)
    se <- res$se
    #changed:
    se_all <- se
    alpha <- res$alpha
    gs <- res$gs
    if (is.null(out.dir)) {
      out.dir <- configEBrowser("OUTDIR.DEFAULT")
      stamp <- format(Sys.time(), "%a_%b%d_%Y_%H%M%S")
      out.dir <- file.path(out.dir, stamp)
    } else { 
      out.dir <- path.expand(out.dir) 
      }
      
    if (!file.exists(out.dir))
      dir.create(out.dir, recursive = TRUE)
    if (is.null(report.name))
      report.name <- method
    if (nr.show < 1)
      nr.show <- res$nr.sigs
    stopifnot(nr.show > 0)
    if (nr.show > nrow(res$res.tbl))
      nr.show <- nrow(res$res.tbl)
    res <- res$res.tbl[seq_len(nr.show),]
    gsc <- EnrichmentBrowser:::.gsList2Collect(gs[res[, 1]])
    if(grepl("^GO", names(gsc)[1]) & ! grepl("^GO:", names(gsc)[1])){
      names(gsc) <- gsub("^GO","GO:",names(gsc))
      for(i in names(gsc)){
        new_i <- gsub("^GO","GO:",i)
        #GSEABase::setName(gsc[[i]]) <- new_i
        GSEABase::description(gsc[[i]]) <- AnnotationDbi::Term(new_i)
        GSEABase::collectionType(gsc[[i]])@ontology <- AnnotationDbi::Ontology(new_i)
      }
    }
    res[, 1] <- vapply(res[, 1], function(s)
      unlist(strsplit(s,
                      "_"))[1], character(1))
    is.kegg <-
      is(GSEABase:::collectionType(gsc[[1]]), "KEGGCollection")
    is.go <- is(GSEABase:::collectionType(gsc[[1]]), "GOCollection")
    gs.title <- vapply(gsc, description, character(1))
    nr.genes <-
      vapply(gsc, function(g)
        length(GSEABase::geneIds(g)), integer(1))
    cnames <- c(colnames(res)[1], "TITLE")
    resn <- DataFrame(res[, 1], gs.title)
    if (!("NR.GENES" %in% colnames(res))) {
      cnames <- c(cnames, "NR.GENES")
      resn <- DataFrame(resn, nr.genes)
    }
    cnames <- c(cnames, colnames(res)[2:ncol(res)])
    resn <- DataFrame(resn, res[, 2:ncol(res)])
    colnames(resn) <- cnames
    res <- resn
    im <- GSEABase::incidence(gsc)
    org <- organism(gsc[[1]])
    if (is.na(org))
      org <- metadata(se)$annotation
    if (!length(org))
      stop(
        "Organism annotation not found!\n",
        "Organism under study must be annotated via metadata(se)$annotation"
      )
    message("Creating gene report ...")
    se <- se[colnames(im),]
    rcols <- sapply(c("FC.COL", "ADJP.COL"), configEBrowser)
    fDat <- rowData(se, use.names = TRUE)[, rcols]
    fDat <- as.data.frame(fDat)
    if( ! configEBrowser("EZ.COL") %in% names(rowData(se))) {
      rowData(se)[, configEBrowser("EZ.COL")] <- rownames(rowData(se))
    }
    if (org == "none") {
      #configEBrowser("EZ.COL") <- "GeneID"
      EZ.COL <- configEBrowser("EZ.COL")
      gn.cols <- sapply(c("SYM.COL", "GN.COL"), configEBrowser)
      gt <- rowData(se)[, c(EZ.COL, gn.cols)]
      gt <- gt[!duplicated(gt[, EZ.COL]), ]
      rownames(gt) <- gt[[EZ.COL]]
      if (!is.null(fDat)) {
        fDat[, 1] <- round(fDat[, 1], digits = 2)
        fDat[, 2] <- signif(fDat[, 2], digits = 2)
        gt <- cbind(gt, fDat)
      }
    } else {
      gt <- suppressMessages(EnrichmentBrowser:::.geneTable(im, org, fcs = fDat))
    }
    gn.cols <- sapply(c("SYM.COL", "GN.COL"), configEBrowser)
    nonames <- is.na(gt[, gn.cols[1]])
    gt[ nonames, gn.cols[1]] <- rownames(gt)[ nonames]
    gt[ ! nonames, gn.cols[1]] <- paste( rownames(gt)[ ! nonames] , gt[ ! nonames, gn.cols[1]], sep = "/")
    rowData(se)[, gn.cols] <- DataFrame(gt[, gn.cols])
    
    gt.reps <-
      sapply(gsc, function(s)
        EnrichmentBrowser:::.geneReport(s, gt, out.dir))
    link <- paste0(names(gsc), ".html")
    res[, "NR.GENES"] <-
      hwriter::hwrite(res[, "NR.GENES"], link = link,
                      table = FALSE)
    if ("LEADING.EDGE" %in% colnames(res)) {
      gsc2 <- gsc
      le_genes <- setNames(res$LEADING.EDGE, res$GENE.SET)
      for (x in names(gsc2)) {
        GSEABase::geneIds(gsc2[[x]]) <- le_genes[[GSEABase::collectionType(gsc2[[x]])@ids]]
        GSEABase::setName(gsc2[[x]]) <-
          paste0(GSEABase::setName(gsc2[[x]]), "_Leading.Edge")
      }
      res$LEADING.EDGE = sapply(res$LEADING.EDGE, length)
      gt2.reps <-
        sapply(gsc2, function(s) {
          gset = sub("_Leading.Edge", "", GSEABase::setName(s))
          
          geneReportLE(s, gt, out.dir, gsc[[gset]] , se_all)
        })  #EnrichmentBrowser:::.geneReport(s, gt, out.dir))
      link2  <- paste0(names(gsc2), ".html")
      res[, "LEADING.EDGE"] <-
        hwriter::hwrite(res[, "LEADING.EDGE"], link = link2, table = FALSE)
    }
    message("Creating set view ...")
    out.prefix <- file.path(out.dir, names(gsc))
    names(out.prefix) <- names(gsc)
    vcol <-
      sapply(gsc, function(s)
        EnrichmentBrowser:::.viewSet(se[GSEABase::geneIds(s),],
                                     out.prefix[GSEABase::setName(s)]))
    vcol <- hwriter::hwriteImage(
      sub("sview.html", "volc.png",
          vcol),
      link = vcol,
      table = FALSE,
      height = 50,
      width = 50,
      target = "_blank"
    )
    res <- DataFrame(res, vcol)
    colnames(res)[ncol(res)] <- "SET.VIEW"
    if (is.kegg) {
      message("Creating kegg view ...")
      isAvailable("pathview", type = "software")
      vcol <-
        sapply(gsc, function(s)
          EnrichmentBrowser:::.viewPath(GSEABase::setName(s),
                                        se[GSEABase::geneIds(s),], out.prefix[GSEABase::setName(s)]))
      vcol <- hwriter::hwriteImage(
        sub("kview.html", "kpath.png",
            vcol),
        link = vcol,
        table = FALSE,
        height = 50,
        width = 50,
        target = "_blank"
      )
      res <- DataFrame(res, vcol)
      colnames(res)[ncol(res)] <- "PATH.VIEW"
    }
    if (!is.null(graph.view)) {
      message("Creating graph view ...")
      vcol <-
        sapply(gsc, function(s)
          EnrichmentBrowser:::.viewGraph(
            se[GSEABase::geneIds(s),],
            EnrichmentBrowser:::.queryGRN(GSEABase::geneIds(s), graph.view, index = FALSE),
            alpha,
            out.prefix[GSEABase::setName(s)]
          ))
      vcol <- hwriter::hwriteImage(
        sub("html$", "png", vcol),
        link = vcol,
        table = FALSE,
        height = 50,
        width = 50,
        target = "_blank"
      )
      res <- DataFrame(res, vcol)
      colnames(res)[ncol(res)] <- "GRAPH.VIEW"
    }
    link <- NULL
    GS.COL <- configEBrowser("GS.COL")
    if (is.kegg) {
      link <-
      sapply(gsc, function(s)
        EnrichmentBrowser:::.getHTMLOfMarkedPathway(GSEABase::setName(s),
                                                    GSEABase::geneIds(s)[fDat[GSEABase::geneIds(s), 2] < alpha])) 
      } else if (is.go) {
      link <- paste0(configEBrowser("GO.SHOW.URL"), res[, GS.COL])
      }
    if (!is.null(link)) {
      res[, GS.COL] <- hwriter::hwrite(res[, GS.COL], link = link,
                                       table = FALSE)
    }
    htmlRep <- ReportingTools::HTMLReport(
      shortName = report.name,
      title = paste(toupper(method), configEBrowser("RESULT.TITLE"),
                    sep = " - "),
      basePath = dirname(out.dir),
      reportDirectory = basename(out.dir)
    )
    res <- as.data.frame(res)
    ReportingTools::publish(res, htmlRep)
    rep <- ReportingTools::finish(htmlRep)
    if (!html.only && interactive()) {
      message(paste("Your output files are in", out.dir))
      message(paste0("HTML report: ", report.name, ".html"))
      if (Sys.getenv("RSTUDIO") == "1")
        rep <- URLencode(rep)
      browseURL(rep)
    }
  }

#environment(eaBrowse_alt) <- as.environment("package:EnrichmentBrowser")

sbea_fgse <-
  function(seset,
           gs,
           alpha = 0.05,
           #perm = 10000,
           padj.method = "BH", 
           stat_name = "DESeq2.STAT",
           minSize=10,
           maxSize=500) {
    sbea.test <-
      sbea(
        method = use.fgse,
        se = seset,
        gs = gs ,
        alpha = alpha,
        #perm = perm,
        padj.method = padj.method,
        minSize=minSize,
        maxSize=maxSize
      )
    sbea.test$method <- "fgse"
    sbea.test$se <- seset
    sbea.test$res.tbl$NR.SIG.GENES <- get.nsig(sbea.test)
    sbea.test$res.tbl$NR.GENES <- get.ngenes(sbea.test)
    sbea.test$res.tbl <-
      cbind(sbea.test$res.tbl, get.nsig_up_down(sbea.test))
    fgse_res <- use.fgse(sbea.test$se, gs, full_res = TRUE)
    sbea.test$res.tbl$PVAL <-
      fgse_res$pval[match(sbea.test$res.tbl$GENE.SET, fgse_res$pathway)]
    sbea.test$res.tbl$ADJ.PVAL <-
      fgse_res$padj[match(sbea.test$res.tbl$GENE.SET, fgse_res$pathway)]
    sbea.test$res.tbl$LEADING.EDGE <-
      fgse_res$leadingEdge[match(sbea.test$res.tbl$GENE.SET, fgse_res$pathway)]
    sbea.test$res.tbl$ES <-
      fgse_res$ES[match(sbea.test$res.tbl$GENE.SET, fgse_res$pathway)]
    sbea.test$res.tbl$NES <-
      fgse_res$NES[match(sbea.test$res.tbl$GENE.SET, fgse_res$pathway)]
    sbea.test$res.tbl <-
      sbea.test$res.tbl[, c(
        "GENE.SET",
        "PVAL",
        "ADJ.PVAL",
        "NR.GENES",
        "LEADING.EDGE",
        "ES",
        "NES",
        "NR.SIG.GENES",
        "N.SIG.UP",
        "N.SIG.DOWN"
      )]
    
    return(sbea.test)
  }

exprsHeatmap_alt <-
  function(expr,
           grp,
           scale.rows = TRUE,
           log.thresh = 100)
  {
    # log-transform?
    dd <- diff(range(expr, na.rm = TRUE))
    if (dd > log.thresh)
      expr <- log(expr + 1, base = 2)
    
    # scale?
    if (scale.rows)
      expr <- t(scale(t(expr)))
    
    # group colors
    grp <- as.factor(grp)
    if (length(levels(grp) > 2)) {
      coll <- rainbow(length(levels(grp)))
      cols_ord = order(grp)
    } else {
      coll <- c("#B62A84", "#2AB68C")
      cols_ord = NULL
    }
    
    names(coll) <- levels(grp)
    coll <- list(Group = coll)
    
    # annotation
    df <- data.frame(Group = grp)
    ha <- ComplexHeatmap::HeatmapAnnotation(df = df, col = coll)
    # remove rows that are all 0 or all equal
    no_change = apply(expr, 1, function(x) length(unique(x))==1 )
    expr = expr[! no_change,]
    # plot
    print(
      ComplexHeatmap::Heatmap(
        expr,
        name = "Expression",
        top_annotation = ha,
        show_row_names = (nrow(expr) < 41),
        show_column_names = (ncol(expr) < 41),
        column_title = "Samples",
        row_title = "Features",
        column_order = cols_ord,
      )
    )
  }



unlockBinding("exprsHeatmap", env = as.environment("package:EnrichmentBrowser"))
assignInNamespace("exprsHeatmap",
                  exprsHeatmap_alt,
                  envir = as.environment("package:EnrichmentBrowser"))
lockBinding("exprsHeatmap", env = as.environment("package:EnrichmentBrowser"))

write.gmt <- function(gs, filename) {
  write(sapply(names(gs), function(x)
    paste(x, x, paste(gs[[x]], collapse = "\t"), sep = "\t")), file = filename)
}


create_heatmap = function(res_tab, dds, main = "", ngenes = 0, alpha = 0.05, scale_rows =TRUE, 
                          column_names_rot = 90, gradient_col = c("blue", "white", "red"),
                          show_annotation_name = TRUE){
  sig_genes = rownames(res_tab)[! is.na(res_tab$FDR) & res_tab$FDR <= alpha]
  if(length(sig_genes) <= 1){return()}
  if (ngenes > 0){
    sig_genes = sig_genes[1:min(length(sig_genes), ngenes)]
  }
  
  exprs = cpm(dds[sig_genes,], 
              log = TRUE, prior.count = 1)
  #exprs = log2(exprs + 1)
  gene_symbols = as.character(dds$genes[sig_genes, ] %>% 
                                mutate(Name = if_else(is.na(Name), GeneID, Name)) %>%
                                pull(Name))
  # gradient cols
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), gradient_col)
  # group colors
  grp <- droplevels(dds$samples$group)
  if (length(levels(grp)) == 2) {
    coll <- c("turquoise","salmon"  )
  } else if (length(levels(grp)) <=12){
    coll <- brewer.pal(n=length(levels(grp)),name="Set3")
  } else {
    coll <- colorRampPalette(brewer.pal(n=12,name="Set3"))(length(levels(grp)))
  }
  #coll <- brewer.pal(length(levels(grp)), "Set1")
  cols_ord = order(grp)
  names(coll) <- levels(grp)
  coll <- list(Group = coll)
  # annotation
  df <- data.frame(Group = grp)
  if(scale_rows){
    exprs = t(scale(t(exprs)))
  }
  ha <- ComplexHeatmap::HeatmapAnnotation(Group = grp, col = coll, 
                                          show_annotation_name = show_annotation_name)
  # remove rows that are all 0 or all equal
  #no_change = apply(exprs, 1, function(x) length(unique(x))==1 )
  #exprs = exprs[! no_change,]
  #exprs
  #return(gene_symbols)
  # plot
  return(ComplexHeatmap::Heatmap(
    exprs,
    col = col_fun,
    column_names_rot = column_names_rot,
    name = if (scale_rows) "scaled log(CPM)" else "log2(CPM)",
    top_annotation = ha,
    show_row_names = (nrow(exprs) <= 50),
    show_column_names = (ncol(exprs) <= 50),
    row_title = NULL,
    column_title = main,
    column_order = cols_ord,
    cluster_rows = TRUE,
    row_labels = gene_symbols,
    row_names_gp = grid::gpar(fontsize = 8),
    clustering_distance_rows = "pearson"
  )
  )
  
}

createMAplot = function(data, main = NULL, alpha = 0.1, gnames = "Name", colors = c("blue", "gray", "red"), 
                        ntop = 45, size = 1, font.size = 8, xlab = "Log2 mean CPM", 
                        ylab = "Log2 fold change", fc = 2, maxoverlaps = 100,
                        ggtheme = theme_classic(),...){
  data[is.na(data$FDR)]$FDR = 1.0
  #data$log2baseMean = log2(data$baseMean + 1)
  data$sig = as.factor(sign(data$logFC)*( data$FDR <= alpha))
  data <- data[order(data$FDR), ]
  nsig_up = sum(data$FDR <= alpha & data$logFC > 0)
  nsig_down = sum(data$FDR <= alpha & data$logFC < 0)
  labs_data <- data %>% filter(logFC > 0) %>% slice_head(n = min(nsig_up,ntop))
  labs_data <- bind_rows(labs_data, 
                         data %>% filter(logFC < 0) %>% slice_head(n = min(nsig_down,ntop)))
  p = ggplot(data, aes(x = logCPM, y = logFC)) +
    geom_point(aes(color = sig), size = size) +  
    scale_color_manual(values=colors, guide = guide_legend(reverse=TRUE), labels = c("down", "-", "up"))
  p = p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = get(gnames)), 
                                   box.padding = unit(0.35, "lines"),
                                   point.padding = unit(0.3, "lines"),
                                   force = 1, max.overlaps = maxoverlaps,
                                   size = font.size/3)
  p = p + scale_x_continuous(breaks=seq(0, max(data$logCPM), 2))+
    labs(x = xlab, y = ylab, title = main, color = "signif.") + # to remove legend title use color = ""
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2),
               color = c("black", "black", "black")) #+ theme(legend.position="bottom", legend.direction = "horizontal")
  p <- ggpar(p, ggtheme = ggtheme, ...)
  return(p)
}

volcano.plot.cust <- function(data, pch=16, cex=0.45,  main = NULL, alpha = 0.05, gnames = "Name", 
                              colors = c("blue", "gray", "red"), 
                              ntop = 45, size = 1, font.size = 8, fc = 2, maxoverlaps = 100,cols=c("blue","grey", "red"), 
                              xlab="Log2(FC)", ylab=expression("-Log10(Padj)"), adjP = "padj", lgfc = "logFC", 
                              ggtheme = theme_bw(), ...){
  data <- data %>% rename(FDR = {{adjP}}, logFC = {{lgfc}}) %>% arrange(FDR)
  data[is.na(data$FDR)]$FDR = 1.0
  #data$log2baseMean = log2(data$baseMean + 1)
  data$sig = as.factor(sign(data$logFC)*( data$FDR <= alpha))
  data <- data[order(data$FDR), ]
  nsig_up = sum(data$FDR <= alpha & data$logFC > 0)
  nsig_down = sum(data$FDR <= alpha & data$logFC < 0)
  labs_data <- data %>% filter(logFC > 0) %>% slice_head(n = min(nsig_up,ntop))
  labs_data <- bind_rows(labs_data, 
                         data %>% filter(logFC < 0) %>% slice_head(n = min(nsig_down,ntop)))
  p = ggplot(data, aes(x = logFC, y = -log10(FDR))) +
    geom_point(aes(color = sig), size = size) +  
    scale_color_manual(values=colors, guide = guide_legend(reverse=TRUE), labels = c("down", "-", "up"), ) +
    scale_y_continuous(expand = expansion(mult = c(0, .05)))
  
  p = p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = get(gnames)), 
                                   box.padding = unit(0.25, "lines"),
                                   point.padding = unit(0.1, "lines"),
                                   force = 1, force_pull = 4, max.overlaps = maxoverlaps,
                                   size = font.size/4)
  p = p +
    labs(x = xlab, y = ylab, title = main, color = "signif.") + # to remove legend title use color = ""
    geom_vline(xintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2),
               color = c("black", "darkgrey", "darkgrey")) 
  p <- ggpar(p, ggtheme = ggtheme, ...)
  return(p + theme(legend.position="right", legend.title = element_blank()))
}

