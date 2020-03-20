library("knitr")
library(BiocStyle)
opts_chunk$set(tidy=FALSE,dev="pdf",fig.show="hide",
               fig.width=4,fig.height=4.5,
               message=FALSE, warning=FALSE)
BiocStyle::latex(width=90)
options(digits=3, width=80, prompt=" ", continue=" ")
## ----simResult, cache=TRUE, results='hide', fig.height=4, fig.width=4---------
# load library
library('variancePartition')

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
data(varPartData)

# Specify variables to consider
# Age is continuous so model it as a fixed effect
# Individual and Tissue are both categorical, 
# so model them as random effects
# Note the syntax used to specify random effects
form <- ~ Age + (1|Individual) + (1|Tissue) + (1|Batch) #???

# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified, 
#     a linear mixed model is used
# If all variables are modeled as fixed effects, 
#		a linear model is used
# each entry in results is a regression model fit on a single gene
# 2) extract variance fractions from each model fit
# for each gene, returns fraction of variation attributable 
#		to each variable 
# Interpretation: the variance explained by each variables 
# after correcting for all other variables
# Note that geneExpr can either be a matrix, 
# and EList output by voom() in the limma package, 
# or an ExpressionSet
varPart <- fitExtractVarPartModel( geneExpr, form, info )

# sort variables (i.e. columns) by median fraction 
#		of variance explained
vp <- sortCols( varPart )

# Figure 1a
# Bar plot of variance fractions for the first 10 genes
plotPercentBars( vp[1:10,] )

#
# Figure 1b
# violin plot of contribution of each variable to total variance
plotVarPart( vp )

## ----accessResults, cache=TRUE, warning=FALSE---------------------------------
# Access first entries
head(varPart)

# Access first entries for Individual
head(varPart$Individual)

# sort genes based on variance explained by Individual
head(varPart[order(varPart$Individual, decreasing=TRUE),])

## ----savePlot, cache=TRUE, eval=FALSE-----------------------------------------
#  fig <- plotVarPart( vp )
#  ggsave(file, fig)

## ----plotStratify, cache=TRUE, warning=FALSE, fig.height=4, fig.width=4-------
# get gene with the highest variation across Tissues
# create data.frame with expression of gene i and Tissue 
#		type for each sample
i <- which.max( varPart$Tissue )
GE <- data.frame( Expression = geneExpr[i,], Tissue = info$Tissue)

# Figure 2a 
# plot expression stratified by Tissue 
plotStratify( Expression ~ Tissue, GE, main=rownames(geneExpr)[i])
#
# get gene with the highest variation across Individuals
# create data.frame with expression of gene i and Tissue 
#		type for each sample
i <- which.max( varPart$Individual )
GE <- data.frame( Expression = geneExpr[i,], 
                  Individual = info$Individual)

# Figure 2b
# plot expression stratified by Tissue 
label <- paste("Individual:", format(varPart$Individual[i]*100, 
                                     digits=3), "%")
main <- rownames(geneExpr)[i]
plotStratify(  Expression ~ Individual, GE, colorBy=NULL, 
               text=label, main=main)

## ----cache=TRUE---------------------------------------------------------------
library('lme4')

# fit regression model for the first gene 
form_test <- geneExpr[1,] ~ Age + (1|Individual) + (1|Tissue)
fit <- lmer(form_test, info, REML=FALSE )

# extract variance statistics
calcVarPart(fit)

## ----canCorPairs, cache=TRUE, results='hide', fig.width=5, fig.height=5-------
form <- ~ Individual + Tissue + Batch + Age + Height

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value  
C = canCorPairs( form, info)

# Plot correlation matrix
plotCorrMatrix( C )

## ----simResult-omit, cache=TRUE, results='hide'-------------------------------

form <- ~ Age + (1|Individual) + (1|Tissue) + (1|Batch)

# Fit model
results <- fitVarPartModel( geneExpr, form, info )

# Extract results
varPart <- extractVarPart( results )

## ----simResult-fast-two-step, cache=TRUE, results='hide'----------------------
# Fit model and run summary() function on each model fit 
vpSummaries <- fitVarPartModel( geneExpr, form, info, fxn=summary )

## ----simResult-fast-two-step1, cache=TRUE-------------------------------------
# Show results of summary() for the first gene
vpSummaries[[1]]

## ----echo=TRUE, cache=TRUE, results='hide'------------------------------------
form <- ~ (1|Tissue) + (1|Individual) + (1|Batch) + Age
varPart <- fitExtractVarPartModel( geneExpr, form, info )

## ----echo=TRUE, cache=TRUE, results='hide'------------------------------------
library('limma')
# subtract out effect of Batch
fit <- lmFit( geneExpr, model.matrix(~ Batch, info))
res <- residuals( fit, geneExpr)

# fit model on residuals
form <- ~ (1|Tissue) + (1|Individual) + Age

varPartResid <- fitExtractVarPartModel( res, form, info )

## ----echo=TRUE, cache=TRUE, results='hide'------------------------------------
# subtract out effect of Batch with linear mixed model
modelFit <- fitVarPartModel( geneExpr, ~ (1|Batch), info )
res <- residuals( modelFit )

# fit model on residuals
form <- ~ (1|Tissue) + (1|Individual) + Age

varPartResid <- fitExtractVarPartModel( res, form, info )

## ----echo=TRUE, cache=TRUE, results='hide'------------------------------------
# extract residuals directly without storing intermediate results
residList <- fitVarPartModel( geneExpr, ~ (1|Batch), info, 
                              fxn=residuals )

# convert list to matrix
residMatrix = do.call(rbind, residList)

## ----withinTissue, echo=TRUE, cache=TRUE, results='hide', fig.height=5, fig.width=4----
# specify formula to model within/between individual variance 
# separately for each tissue
# Note that including +0 ensures each tissue is modeled explicitly
# Otherwise, the first tissue would be used as baseline
form <- ~ (Tissue+0|Individual) + Age + (1|Tissue) + (1|Batch)

# fit model and extract variance percents
varPart <- fitExtractVarPartModel( geneExpr, form, info, showWarnings=FALSE )

# violin plot
plotVarPart( sortCols(varPart), label.angle=60 )

## ----echo=TRUE, cache=TRUE, results='hide'------------------------------------
form <- ~ (1|Individual) + (1|Tissue) + Age + Height

# fit model
res <- fitVarPartModel( geneExpr[1:4,], form, info )

## ----echo=TRUE, cache=TRUE----------------------------------------------------
# evaluate the collinearity score on the first model fit
# this reports the correlation matrix between coefficient estimates
# for fixed effects
# the collinearity score is the maximum absolute correlation value
# If the collinearity score > .99 then the variance partition 
# estimates may be problematic
# In that case, a least one variable should be omitted
colinearityScore( res[[1]] )

## ----echo=TRUE, cache=TRUE, results='hide'------------------------------------

form <- ~ (1|Individual) + (1|Tissue) + Age + Height

# Specify custom weights
# In this example the weights are simulated from a 
# uniform distribution and are not meaningful.
weights <- matrix(runif(length(geneExpr)), nrow=nrow(geneExpr))

# Specify custom weights
res <- fitExtractVarPartModel( geneExpr[1:4,], form, info, 
                               weightsMatrix=weights[1:4,] )

## ----vpInteraction, echo=TRUE, cache=TRUE, results='hide', fig.width=4, fig.height=4----
form <- ~ (1|Individual) + Age + Height + (1|Tissue) + (1|Batch) + 
  (1|Batch:Tissue)

# fit model
vpInteraction <- fitExtractVarPartModel( geneExpr, form, info )

plotVarPart( sortCols( vpInteraction ) )

## ----echo=TRUE, cache=TRUE, results='hide'------------------------------------
library('limma')
library('edgeR')

# identify genes that pass expression cutoff
isexpr <- rowSums(cpm(geneCounts)>1) >= 0.5 * ncol(geneCounts)

# create data structure with only expressed genes
gExpr <- DGEList(counts=geneCounts[isexpr,])

# Perform TMM normalization
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of 
# uncertainty.
# Recommend including variables with a small number of categories 
# that explain a substantial amount of variation
design <- model.matrix( ~ Batch, info)

# Estimate precision weights for each gene and sample
# This models uncertainty in expression measurements
vobjGenes <- voom(gExpr, design )

# Define formula
form <- ~ (1|Individual) + (1|Tissue) + (1|Batch) + Age

# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE
varPart <- fitExtractVarPartModel( vobjGenes, form, info )

## ----echo=TRUE, cache=TRUE, results='hide'------------------------------------
library('DESeq2')

# create DESeq2 object from gene-level counts and metadata
dds <- DESeqDataSetFromMatrix(countData = geneCounts,
                              colData = info,
                              design = ~ 1)

# Estimate library size correction scaling factors
dds <- estimateSizeFactors(dds)

# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds)>1) >= 0.5 * ncol(dds)

# compute log2 Fragments Per Million
# Alternatively, fpkm(), vst() or rlog() could be used
quantLog <- log2( fpm( dds )[isexpr,] + 1)

# Define formula
form <- ~ (1|Individual) + (1|Tissue) + (1|Batch) + Age

# Run variancePartition analysis
varPart <- fitExtractVarPartModel( quantLog, form, info)

## ----echo=TRUE, cache=TRUE, results='hide', eval=FALSE------------------------
#  library('tximportData')
#  library('tximport')
#  library('readr')
#  
#  # Get data from folder where tximportData is installed
#  dir <- system.file("extdata", package = "tximportData")
#  samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
#  files <- file.path(dir, "kallisto", samples$run, "abundance.tsv")
#  names(files) <- paste0("sample", 1:6)
#  
#  tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
#  
#  # reads results from kallisto
#  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene,
#  	countsFromAbundance = "lengthScaledTPM")
#  
#  # define metadata (usually read from external source)
#  info_tximport <- data.frame( Sample = sprintf("sample%d", 1:6),
#  	Disease=c("case", "control")[c(rep(1, 3), rep(2, 3) )] )
#  
#  # Extract counts from kallisto
#  y <- DGEList( txi$counts )
#  
#  # compute library size normalization
#  y <- calcNormFactors(y)
#  
#  # apply voom to estimate precision weights
#  design <- model.matrix( ~ Disease, data = info_tximport)
#  vobj <- voom(y, design)
#  
#  # define formula
#  form <- ~ (1|Disease)
#  
#  # Run variancePartition analysis (on only 10 genes)
#  varPart_tx <- fitExtractVarPartModel( vobj[1:10,], form,
#  	info_tximport)
#  

## ----echo=TRUE, cache=TRUE, results='hide'------------------------------------
library('ballgown')

# Get data from folder where ballgown is installed
data_directory <- system.file('extdata', package='ballgown') 

# Load results of Cufflinks/Tablemaker
bg <- ballgown(dataDir=data_directory, samplePattern='sample', 
               meas='all')

# extract gene-level FPKM quantification
# Expression can be convert to log2-scale if desired
gene_expression <- gexpr(bg)

# extract transcript-level FPKM quantification
# Expression can be convert to log2-scale if desired
transcript_fpkm <- texpr(bg, 'FPKM')

# define metadata (usually read from external source)
info_ballgown <- data.frame( Sample = sprintf("sample%02d", 1:20), 
                             Batch = rep(letters[1:4], 5), 
                             Disease=c("case", "control")[c(rep(1, 10), rep(2, 10) )] )

# define formula
form <- ~ (1|Batch) + (1|Disease)

# Run variancePartition analysis 
# Gene-level analysis
varPart_gene <- fitExtractVarPartModel( gene_expression, form, 
                                        info_ballgown)

# Transcript-level analysis
varPart_transcript <- fitExtractVarPartModel( transcript_fpkm, form, 
                                              info_ballgown)

## ----echo=FALSE, cache=TRUE, results='hide'-----------------------------------
library('variancePartition')

sim_data = function( n, p, var_indiv, var_site){
  
  info = data.frame(ID = paste("s", 1:n, sep=''), 
                    Individual = factor(rep( round(seq(1, n/6, length.out=n/6)),6)), 
                    Tissue = factor(sort(rep(toupper(letters[1:3]), n/3))),
                    Age = rpois(n, 50))
  
  geneExpr = matrix(0, nrow=p, ncol=n)
  for( i in 1:p){
    beta_indiv = rnorm(nlevels(info$Individual), 0, sqrt(var_indiv))
    beta_site = rnorm(nlevels(info$Tissue), 0, sqrt(var_site))
    eta = model.matrix( ~ Individual, info) %*% beta_indiv + model.matrix( ~ Tissue+0, info) %*% beta_site
    
    noise = rbeta(1, 10, 100)
    errVar = var(eta) * (noise) / (1-noise)
    
    geneExpr[i,] = eta + rnorm(n, 0, sqrt(errVar))
  }
  colnames(geneExpr) = paste("s", 1:ncol(geneExpr), sep='')
  rownames(geneExpr) = paste("gene", 1:nrow(geneExpr), sep='')
  
  return( list(info=info, geneExpr=geneExpr))
}

plotVar = function( geneExpr, info ){
  
  form =  ~ (1|Individual) + (1|Tissue)
  varPart = fitExtractVarPartModel( geneExpr, form, info )
  plotVarPart( varPart, col=rainbow(8)[1:3] )
}

plotVarCross = function( geneExpr, info, label.angle=30 ){
  
  form =  ~ (Tissue+0|Individual) + (1|Tissue)
  
  #res = fitVarPartModel( geneExpr, form, info )
  #varPart = extractVarPart( res )
  
  varPart = fitExtractVarPartModel( geneExpr, form, info, showWarnings=FALSE )
  
  plotVarPart( varPart, label.angle=label.angle )
}

plotPCA = function(geneExpr, col){
  
  dcmp = prcomp(t(geneExpr))
  
  par(mar =  c(4, 4, 1, 1) + 0.1)
  plot(dcmp$x[,1:2], col=col)
  
}

plotTree = function(geneExpr, col){
  
  hc = hclust(dist(t(geneExpr)))
  library('dendextend')
  hcd = as.dendrogram(hc)
  labels_colors(hcd) <- col[match(labels(hcd), names(col))]
  
  par(mar =  c(4, 4, 1, 1) + 0.1, cex=.6)
  plot( hcd, yaxt='n', horiz=TRUE )
}

## ----siteDominant, echo=FALSE, cache=TRUE, results='hide', fig.height=5, fig.width=5----
set.seed(1)

n = 60
p = 200

data = sim_data( n, p, 1, 4)

colTissue = rainbow(nlevels(data$info$Tissue))[data$info$Tissue]
names(colTissue) = data$info$ID

colIndiv = rainbow(nlevels(data$info$Individual))[data$info$Individual]
names(colIndiv) = data$info$ID

plotPCA( data$geneExpr, colTissue)

plotPCA( data$geneExpr, colIndiv)

plotTree( data$geneExpr, colTissue)
legend("topleft", legend=levels(data$info$Tissue), fill=rainbow(nlevels(data$info$Tissue))[1:nlevels(data$info$Tissue)], title="Tissue")

plotTree( data$geneExpr, colIndiv)
legend("topleft", legend=levels(data$info$Individual), fill=rainbow(nlevels(data$info$Individual))[1:nlevels(data$info$Individual)], title="Individual")

plotVar( data$geneExpr, data$info )

plotVarCross( data$geneExpr, data$info, label.angle=60 )

## ----IndivDominant, echo=FALSE, cache=TRUE, results='hide', fig.height=5, fig.width=5----
set.seed(1)

n = 60
p = 200

data = sim_data( n, p, 3, 1)

colTissue = rainbow(nlevels(data$info$Tissue))[data$info$Tissue]
names(colTissue) = data$info$ID

colIndiv = rainbow(nlevels(data$info$Individual))[data$info$Individual]
names(colIndiv) = data$info$ID

plotPCA( data$geneExpr, colTissue)

plotPCA( data$geneExpr, colIndiv)

plotTree( data$geneExpr, colTissue)
legend("topleft", legend=levels(data$info$Tissue), fill=rainbow(nlevels(data$info$Tissue))[1:nlevels(data$info$Tissue)], title="Tissue")

plotTree( data$geneExpr, colIndiv)
legend("topleft", legend=levels(data$info$Individual), fill=rainbow(nlevels(data$info$Individual))[1:nlevels(data$info$Individual)], title="Individual")

plotVar( data$geneExpr, data$info )

plotVarCross( data$geneExpr, data$info, label.angle=60 )

## ----echo=FALSE, cache=FALSE--------------------------------------------------
# if( "cluster" %in% class(cl) ){
if( exists("cl") ){
  library('doParallel')
  
  # stop cluster and catch warning if invalid
  res = tryCatch( {stopCluster(cl)}, warning = function(x) {
  }, error = function(x) {
  }, finally={
  })
  # warning("STOP CLUSTER!!!!!!!!!!!!!!!!!\n")
}

## ----DE, echo=FALSE, cache=TRUE, fig.height=5, fig.width=5--------------------
set.seed(1) 

library('variancePartition')

n = 500
data = data.frame(Sex = c(rep('F', n), rep('M', n)))
data$expression = rnorm(2*n, (as.integer(data$Sex)-1)*2, .3)
data$Sex = factor(data$Sex)

fit = lm(expression ~ Sex, data)
# calcVarPart(fit)
# coef(fit)[2]

plotStratify( expression ~ Sex, data, ylim=c(-6, 9)) +
  annotate("text", x = 0.5, y = 9, hjust=0, label = paste("fold change:", format(coef(fit)[2], digits=3)), size=5.5) +
  annotate("text", x = 0.5, y = 8.2, hjust=0, label = paste("% variance of expression:", format(calcVarPart(fit)[1]*100, digits=3), "%"), size=5.5)


n = 500
data = data.frame(Sex = c(rep('F', n), rep('M', n)))
data$expression = rnorm(2*n, (as.integer(data$Sex)-1)*2, 2.01)
data$Sex = factor(data$Sex)

fit = lm(expression ~ Sex, data)
# calcVarPart(fit)
# coef(fit)[2] 

plotStratify( expression ~ Sex, data, ylim=c(-6, 9)) +
  annotate("text", x = 0.5, y = 9, hjust=0, label = paste("fold change:", format(coef(fit)[2], digits=3)), size=5.5) +
  annotate("text", x = 0.5, y = 8.2, hjust=0, label = paste("% variance of expression:", format(calcVarPart(fit)[1]*100, digits=3), "%"), size=5.5)

## ----sessInfo, results="asis", echo=FALSE-------------------------------------
toLatex(sessionInfo())

## ----resetOptions, results="hide", echo=FALSE---------------------------------
options(prompt="> ", continue="+ ")
