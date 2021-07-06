##### PCA with nucleotide base #####
source("mitogenome_visualization.R")
source("functions.R")

# Create normalised columns of nucleotide percentages of each accession number by length
weighted_PCG <- get_nucleotide_weighted_mean(full_PCG) 

# Creates unique values of weighted PCG nucleotide content
nucleotide_vector <- get_nucleotide_frequencies("Accession", weighted_PCG, unique, 'Taxonomy_Rank') 

# PCA with prcomp
nucleotide_vector_complete <- subset(nucleotide_vector, nucleotide_vector$Taxonomy_Rank != 'Batoidea incertae sedis')
pc <- prcomp(nucleotide_vector_complete[,2:13], center = TRUE, scale = TRUE)
summary(pc)


library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)

ggbiplot(pc, ellipse = TRUE, obs.scale = 1,  var.axes = TRUE, var.scale = TRUE, 
         groups = factor(nucleotide_vector_complete$Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes',
                                                                              'Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes',
                                                                              'Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes',
                                                                              'Chimaeriformes'))) +
    scale_color_manual(name = 'Taxonomic Rank', values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32',
                                                         '#3cb371','#006400','#7b68ee'))


#Create a subset only for members of Chimaeriformes
chimera_subset <- subset(PCG, PCG$Taxonomy_Rank == "Chimaeriformes")

# Create normalised columns of nucleotide percentages of each accession number by length
weighted_chimaeriformes <- get_nucleotide_weighted_mean(chimera_subset) 

# Creates unique values of weighted PCG nucleotide content
nucleotide_vector_chimaera <- get_nucleotide_frequencies("Accession", weighted_chimaeriformes, unique, 'Family')

pc <- prcomp(nucleotide_vector_chimaera[,2:13], center = TRUE, scale = TRUE)
ggbiplot(pc, ellipse = TRUE, obs.scale = 1,  var.axes = FALSE, var.scale = TRUE, 
         groups = factor(nucleotide_vector_chimaera$Family, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae'))) +
    scale_color_manual(name = 'Family', values=c("#00AFBB", "#E7B800", "#FC4E07"))

# With the FactoMineR

#install.packages("FactoMineR")
library(FactoMineR)

nucleotide.pca <- PCA(nucleotide_vector_chimaera[,2:13], scale.unit = TRUE)
ggbiplot(nucleotide.pca, ellipse = TRUE,obs.scale = 1, var.scale = 1, var.axes = TRUE, groups = nucleotide_vector_chimaera$Family) +
    theme(legend.position = "bottom")

###### Clustering ######

# Isolating GC_content on whole mitogenomes
new <- full_mitogenomes %>%
    group_by(Taxonomy_Rank) %>%
    mutate(GC_norm = weighted.mean(GC_content, Length))

full_mito_gc <- aggregate(GC_norm ~ Taxonomy_Rank, new, unique)
tax <- aggregate(Taxonomy_Rank ~ Name, full_mitogenomes, unique)
rownames(full_mito_gc) <- full_mito_gc[,1]
full_mito_gc$Taxonomy_Rank <- NULL
colnames(full_mito_gc) <- paste('Whole Mitogenome')

dist_wholemito <- dist(full_mito_gc, method = 'euclidean')
hclust_avg_wholemito <- hclust(dist_wholemito, method = 'average')
plot(hclust_avg_wholemito)

# Isolating GC content on PCGs
new <- full_PCG %>%
    group_by(Taxonomy_Rank) %>%
    mutate(GC_norm = weighted.mean(GC_content, Length))

full_gene_gc <- aggregate(GC_norm ~ Taxonomy_Rank, new, unique)
rownames(full_gene_gc) <- full_gene_gc[,1]
full_gene_gc <- subset(full_gene_gc, full_gene_gc$Taxonomy_Rank != 'Batoidea incertae sedis')
full_gene_gc$Taxonomy_Rank <- NULL
colnames(full_gene_gc) <- paste('PCG')

dist_pcg <- dist(full_gene_gc, method = 'euclidean')
hclust_avg_pcg <- hclust(dist_pcg, method = 'average')
plot(hclust_avg_pcg)


# Isolating GC content on rRNA genes
ribosomal_rnas <- c('rrnS','rrnL')

# Isolating all rRNA genes
rRNAs <- subset(Final_Results, Final_Results$State != 'partial')
rRNAs <- subset(Final_Results, Final_Results$Content %in% ribosomal_rnas)

# rRNA cluster analyses
new <- rRNAs %>%
    group_by(Taxonomy_Rank) %>%
    mutate(GC_norm = weighted.mean(GC_content, Length))

full_rrna_gc <- aggregate(GC_norm ~ Taxonomy_Rank, new, unique)
rownames(full_rrna_gc) <- full_rrna_gc[,1]
#full_gene_gc <- subset(full_gene_gc, full_gene_gc$Taxonomy_Rank != 'Batoidea incertae sedis')
full_rrna_gc$Taxonomy_Rank <- NULL
colnames(full_rrna_gc) <- paste('rRNAs')

dist_rrna <- dist(full_rrna_gc, method = 'euclidean')
hclust_avg_rrna <- hclust(dist_rrna, method = 'average')
plot(hclust_avg_rrna)


# Isolating GC content on tRNA genes
ribosomal_tnas <- c('tRNA-Phe','tRNA-Val','tRNA-Leu','tRNA-Ile','tRNA-Gln','tRNA-Met','tRNA-Trp','tRNA-Ala','tRNA-Asn','tRNA-Cys','tRNA-Tyr','tRNA-Ser','tRNA-Asp','tRNA-Lys','tRNA-Gly','tRNA-Arg','tRNA-His','tRNA-Ser2','tRNA-Leu2','tRNA-Glu','tRNA-Thr','tRNA-Pro')

# Isolating all tRNA genes
tRNAs <- subset(Final_Results, Final_Results$State != 'partial')
tRNAs <- subset(Final_Results, Final_Results$Content %in% ribosomal_tnas)

# tRNA cluster analyses
new <- tRNAs %>%
    group_by(Taxonomy_Rank) %>%
    mutate(GC_norm = weighted.mean(GC_content, Length))

full_trna_gc <- aggregate(GC_norm ~ Taxonomy_Rank, new, unique)
rownames(full_trna_gc) <- full_trna_gc[,1]
#full_gene_gc <- subset(full_gene_gc, full_gene_gc$Taxonomy_Rank != 'Batoidea incertae sedis')
full_trna_gc$Taxonomy_Rank <- NULL
colnames(full_trna_gc) <- paste('tRNAs')

dist_trna <- dist(full_trna_gc, method = 'euclidean')
hclust_avg_trna <- hclust(dist_trna, method = 'average')
plot(hclust_avg_trna, cex = 0.5)

full_mito_gc$PCG <- full_gene_gc$PCG
full_mito_gc$rRNAs <- full_rrna_gc$rRNAs
full_mito_gc$tRNAs <- full_trna_gc$tRNAs

#install.packages("gplots")

library("gplots")


some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(90, 30), l = c(30, 90), power = c(1/5, 1.5)))

par(cex.main=0.7)
heatmap.2(x=as.matrix(full_mito_gc),
          main = 'Heatmap for GC content',
          distfun = function(x) dist(x, method="euclidean"),
          hclustfun = function(x) hclust(x, method="average") ,
          srtCol = 20,
          dendrogram ='row',
          Colv = 'NA',
          Rowv = TRUE,
          trace = 'none',
          margins = c(10,10),
          key.xlab = '%',
          key.title = NULL,
          denscol = 'grey',
          density.info = 'density',
          col = some_col_func,
          cexCol = 1.2,
          keysize = 1.2)

# Neutrality plots
# GC content in the first two positions of codons vs GC content in the last position

# Getting new columns in full_PCG for GC3 and GC12
full_PCG$GC3 <- rowSums(cbind(full_PCG$`%G3`,full_PCG$`%C3`))

full_PCG$GC12 <- rowSums(cbind(full_PCG$`%G1`, full_PCG$`%C1`, full_PCG$`%G2`, full_PCG$`%C2`)/2)

# creating linear models

# For all data
lm_fit <- lm(GC12 ~ GC3, data=full_PCG)
summary(lm_fit)


ggplotRegression(lm_fit)

# Para Chimaeriformes

lm_carcharhiniformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Carcharhiniformes'))
ggplotRegression(lm_carcharhiniformes)

lm_chimaeriformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Chimaeriformes'))
ggplotRegression(lm_chimaeriformes)

ggplot(data = subset(full_PCG, full_PCG$Taxonomy_Rank == 'Chimaeriformes'), aes(x = GC3, y = GC12, color = Family)) +
    geom_point() +
    geom_smooth(method = "lm", col="red")


lm_heterodontiformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Heterodontiformes'))
ggplotRegression(lm_heterodontiformes)

lm_hexanchiformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Hexanchiformes'))
ggplotRegression(lm_hexanchiformes)

lm_lamniformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Lamniformes'))
ggplotRegression(lm_lamniformes)

lm_myliobatiformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Myliobatiformes'))
ggplotRegression(lm_myliobatiformes)

lm_orectolobiformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Orectolobiformes'))
ggplotRegression(lm_orectolobiformes)

lm_pristiophoriformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Pristiophoriformes'))
ggplotRegression(lm_pristiophoriformes)

lm_rajiformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Rajiformes'))
ggplotRegression(lm_rajiformes)

lm_rhinopristiformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Rhinopristiformes'))
ggplotRegression(lm_rhinopristiformes)

lm_squaliformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Squaliformes'))
ggplotRegression(lm_squaliformes)

lm_squatiniformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Squatiniformes'))
ggplotRegression(lm_squatiniformes)

lm_torpediniformes <- lm(GC12 ~ GC3, data=subset(full_PCG, full_PCG$Taxonomy_Rank == 'Torpediniformes'))
ggplotRegression(lm_torpediniformes)



