#load libraries
library(ape)
library(readxl)
library(tidyverse)
library(phylolm)

#load data (see FIGSHARE repository)
tree <- read.tree('SpeciesTree.txt')
traits <- read.delim('TraitMatrix.txt')
GCN <- read.delim('GeneFamilySize.txt')

#get eigenvectors from the variance-covariance matrix of the tree
projections = eigen(vcv(tree))$vectors
colnames(projections)=paste0("PC",1:length(tree$tip.label))
projections=bind_cols(tibble(Species=tree$tip.label), as_tibble(projections))

#merge everybody together
df <- merge(merge(traits, GCN, by=0), projections, by.x='Row.names', by.y="Species")
rownames(df) <- df$Row.names
df <- df %>% select(-Row.names)

#modified logistic regression with 2 leading eigenvectors included
run_models <- function(trait, OG) {
    formula <- as.formula(paste(trait, "~ sqrt(", OG, ")+", paste(names(df)[14767:14768], collapse= "+")))
    out <- try(phyloglm(formula, df, tree, method='logistic_MPLE'))
    try(cat(paste(trait, OG, summary(out)$coef[2,3], summary(out)$coef[2,4], "\r", sep = '\t'), file = 'lm_output.txt', append = TRUE))
}

#argument list of every trait+gene family combination
args <- expand.grid(trait=names(df)[1:56], OG=names(df)[57:14766])

#run
mapply(run_models, trait=paste(args$trait), OG=paste(args$OG), SIMPLIFY = FALSE)

#read in output to calculate false discovery rate
results <- read.delim('lm_output.txt', header=FALSE, col.names = c('trait', 'OG', 'z', 'p'))
results$fdr <- p.adjust(results$p, method = 'fdr')
write.delim('results.txt', quote=FALSE, row.names=FALSE, delim='\t')




