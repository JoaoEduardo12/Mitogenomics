##### PCA with nucleotide base #####
source("mitogenome_visualization.R")
source("functions.R")


##### Pairwise kolmogorov smirnov tests on RSCU values from different Orders #####

ks_class <- pairwise_statistical_test(class_cu, "mean_values", "ks.test")
ks_gene <- pairwise_statistical_test(gene_cu, "mean_values", "ks.test")

##### PCA on GC content #####
pc_nuc <- prcomp(nucleotide_vector_complete[,2:13], scale. = TRUE)
summary(pc_nuc)

# Scree plot
fviz_eig(pc_nuc)

p_pc_nuc <- fviz_pca_ind(pc_nuc,
                         axes = c(1,2),
             col.ind =  factor(nucleotide_vector_complete$Order, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes',
                                                                                    'Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes',
                                                                                    'Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes',
                                                                                    'Chimaeriformes')),
             label = "none",
             pointsize = 2.5,
             geom = "point",
             mean.point = FALSE,
             col.ind.sup = c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32',
                             '#3cb371','#006400','#7b68ee'),
             repel = TRUE 
)

p_pc_nuc + scale_color_manual(name = "Order", values = c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32',
                                  '#3cb371','#006400','#7b68ee')) + scale_shape_manual(values = c(rep(16, 13))) +
    guides(fill = "none" , shape = "none") 

## PCA loadings
fviz_pca_var(pc_nuc, col.var="steelblue")


###### Clustering on GC content ######

# Isolating GC_content on whole mitogenomes
weighted_mitogenomes_gc <- mitogenomes %>%
    group_by(Order) %>%
    mutate(GC_norm = weighted.mean(GC_content, Length))

mito_gc <- aggregate(GC_norm ~ Order, weighted_mitogenomes_gc, unique)
tax <- aggregate(Order ~ Name, mitogenomes, unique)
rownames(mito_gc) <- mito_gc[,1]
mito_gc$Order <- NULL
colnames(mito_gc) <- paste('Whole Mitogenome')

dist_wholemito <- dist(mito_gc, method = 'euclidean')
hclust_avg_wholemito <- hclust(dist_wholemito, method = 'average')
plot(hclust_avg_wholemito)

# Isolating GC content on PCGs
weighted_pcg_gc <- PCG %>%
    group_by(Order) %>%
    mutate(GC_norm = weighted.mean(GC_content, Length))

gene_gc <- aggregate(GC_norm ~ Order, weighted_pcg_gc, unique)
rownames(gene_gc) <- gene_gc[,1]
gene_gc <- subset(gene_gc, gene_gc$Order != 'Batoidea incertae sedis')
gene_gc$Order <- NULL
colnames(gene_gc) <- paste('PCG')

dist_pcg <- dist(gene_gc, method = 'euclidean')
hclust_avg_pcg <- hclust(dist_pcg, method = 'average')
plot(hclust_avg_pcg)


# Isolating GC content on rRNA genes
ribosomal_rnas <- c('rrnS','rrnL')

# Isolating all rRNA genes
rRNAs <- subset(dataset, dataset$Content %in% ribosomal_rnas)

# rRNA cluster analyses
weighted_rrna_gc <- rRNAs %>%
    group_by(Order) %>%
    mutate(GC_norm = weighted.mean(GC_content, Length))

rrna_gc <- aggregate(GC_norm ~ Order, weighted_rrna_gc, unique)
rownames(rrna_gc) <- rrna_gc[,1]
#full_gene_gc <- subset(full_gene_gc, full_gene_gc$Taxonomy_Rank != 'Batoidea incertae sedis')
rrna_gc$Order <- NULL
colnames(rrna_gc) <- paste('rRNAs')

dist_rrna <- dist(rrna_gc, method = 'euclidean')
hclust_avg_rrna <- hclust(dist_rrna, method = 'average')
plot(hclust_avg_rrna)


# Isolating GC content on tRNA genes
ribosomal_trnas <- c('tRNA-Phe','tRNA-Val','tRNA-Leu','tRNA-Ile','tRNA-Gln','tRNA-Met','tRNA-Trp','tRNA-Ala','tRNA-Asn','tRNA-Cys','tRNA-Tyr','tRNA-Ser','tRNA-Asp','tRNA-Lys','tRNA-Gly','tRNA-Arg','tRNA-His','tRNA-Ser2','tRNA-Leu2','tRNA-Glu','tRNA-Thr','tRNA-Pro')

# Isolating all tRNA genes
tRNAs <- subset(dataset, dataset$Content %in% ribosomal_trnas)

# tRNA cluster analyses
weighted_trna_gc <- tRNAs %>%
    group_by(Order) %>%
    mutate(GC_norm = weighted.mean(GC_content, Length))

trna_gc <- aggregate(GC_norm ~ Order, weighted_trna_gc, unique)
rownames(trna_gc) <- trna_gc[,1]
#full_gene_gc <- subset(full_gene_gc, full_gene_gc$Taxonomy_Rank != 'Batoidea incertae sedis')
trna_gc$Order <- NULL
colnames(trna_gc) <- paste('tRNAs')

dist_trna <- dist(trna_gc, method = 'euclidean')
hclust_avg_trna <- hclust(dist_trna, method = 'average')
plot(hclust_avg_trna, cex = 0.5)

mito_gc$PCG <- gene_gc$PCG
mito_gc$rRNAs <- rrna_gc$rRNAs
mito_gc$tRNAs <- trna_gc$tRNAs

#install.packages("gplots")

mito_gc <- mito_gc[-8,] # subsetting out Pristiophoriformes

library("gplots")

some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(90, 30), l = c(30, 90), power = c(1/5, 1.5)))

par(cex.main=0.7)
heatmap.2(x=as.matrix(mito_gc),
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

####### Neutrality plots #######
# GC content in the first two positions of codons vs GC content in the last position

# Getting new columns in full_PCG for GC3 and GC12
PCG$GC3 <- rowSums(cbind(PCG$G3,PCG$C3))

PCG$GC12 <- rowSums(cbind(PCG$G1, PCG$C1, PCG$G2, PCG$C2)/2)

# For all data

full_PCG_wo_atp8 = subset(full_PCG, full_PCG$Content != "atp8")


# Data extraction from GC12 and GC3 and standardization
GC_df <- cbind(PCG$GC12, PCG$GC3)
colnames(GC_df) <- c("GC12","GC3")
GC_p <- preProcess(GC_df, method = c("range"))
predict_values <- predict(GC_p,GC_df)
predict_values <- data.frame(predict_values)
PCG$GC12_norm <- predict_values$GC12
PCG$GC3_norm <- predict_values$GC3

# For all genes from every Order
lm_fit <- lm(GC12_norm ~ GC3_norm, data=PCG)
summary(lm_fit)
ggplotRegression(lm_fit, "All Genes")

# Carcharhiniformes
lm_carcharhiniformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Carcharhiniformes'))
ggplotRegression2(lm_carcharhiniformes)

lm_chimaeriformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Chimaeriformes'))
ggplotRegression2(lm_chimaeriformes)

# This was for separating into the different Chimaeriformes families
#ggplot(data = subset(PCG, PCG$Order == 'Chimaeriformes'), aes(x = GC3, y = GC12, color = Family)) +
#    geom_point() +
#    geom_smooth(method = "lm", col="red")

# Heterodontiformes
lm_heterodontiformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Heterodontiformes'))
ggplotRegression2(lm_heterodontiformes)

# Hexanchiformes
lm_hexanchiformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Hexanchiformes'))
ggplotRegression2(lm_hexanchiformes)

# Lamniformes
lm_lamniformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Lamniformes'))
ggplotRegression2(lm_lamniformes)

# Myliobatiformes
lm_myliobatiformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Myliobatiformes'))
ggplotRegression2(lm_myliobatiformes)

# Orectolobiformes
lm_orectolobiformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Orectolobiformes'))
ggplotRegression2(lm_orectolobiformes)

# Pristiophoriformes
lm_pristiophoriformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Pristiophoriformes'))
ggplotRegression2(lm_pristiophoriformes)

# Rajiformes
lm_rajiformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Rajiformes'))
ggplotRegression2(lm_rajiformes)

# Rhinopristiformes
lm_rhinopristiformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Rhinopristiformes'))
ggplotRegression2(lm_rhinopristiformes)

# Squaliformes
lm_squaliformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Squaliformes'))
ggplotRegression2(lm_squaliformes)

# Squatiniformes
lm_squatiniformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Squatiniformes'))
ggplotRegression2(lm_squatiniformes)

# Torpediniformes
lm_torpediniformes <- lm(GC12_norm ~ GC3_norm, data=subset(PCG, PCG$Order == 'Torpediniformes'))
ggplotRegression2(lm_torpediniformes)


####### PCA  on RSCU values #######
nucleotide_vector_complete <- subset(group_codon_usage, group_codon_usage$Order != 'Batoidea incertae sedis')
group_codon_usage <- subset(group_codon_usage, group_codon_usage$Order != 'NA')
pc_rscu <- prcomp(group_codon_usage[,3:61], scale. = TRUE)
summary(pc_rscu)

# Scree plot
fviz_eig(pc_rscu)

p_pc_rscu <- fviz_pca_ind(pc_rscu,
                         axes = c(1,2),
                         col.ind =  factor(group_codon_usage$Order, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes',
                                                                                                'Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes',
                                                                                                'Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes',
                                                                                                'Chimaeriformes')),
                         label = "none",
                         pointsize = 2.5,
                         geom = "point",
                         mean.point = FALSE,
                         col.ind.sup = c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32',
                                         '#3cb371','#006400','#7b68ee'),
                         repel = TRUE 
)


p_pc_rscu + scale_color_manual(name = "Order", values = c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32',
                                                         '#3cb371','#006400','#7b68ee')) + scale_shape_manual(values = c(rep(16, 13))) +
    guides(fill = "none" , shape = "none") 

codon_df <- cbind(rownames(pc_rscu$rotation))
colnames(codon_df) <- "codons"

codon_df <- data.frame(codon_df)
codon_df <- get_gc_info_codons(codon_df,"codons","colors")

fviz_pca_var(pc_rscu, axes = c(1,2), col.var = factor(codon_df$colors, levels = c("GC","AT"))) + scale_color_manual(name = "AT/GC encoding codons", values = c("#ff6961","#77dd77"))





###### Heatmap - clustering ########

# Preparing the dataset

# plot_all_codon <- plot_codon_heatmap(all_codon_usage,my_palette,color_codes_all, codon_color_codes, "row")
# 
# atp6_usage <- make_gene_codondataset(group_gene_usage, "atp6")
# acession_per_atp6 <- aggregate(Taxonomy_Rank ~ Accession, atp6_usage, unique)
# atp6_usage[["Accession"]] <- NULL
# atp6_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_atp6 <- c()
# for (element in acession_per_atp6[,2]) {
#     color_codes_atp6 <- c(color_codes_atp6, values(h, key = element))
# }
# plot_atp6 <- plot_codon_heatmap(atp6_usage,my_palette,color_codes_atp6, codon_color_codes[2:59],"column")
# 
# atp8_usage <- make_gene_codondataset(group_gene_usage, "atp8")
# acession_per_atp8 <- aggregate(Taxonomy_Rank ~ Accession, atp8_usage, unique)
# atp8_usage[["Accession"]] <- NULL
# atp8_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_atp8 <- c()
# for (element in acession_per_atp8[,2]) {
#     color_codes_atp8 <- c(color_codes_atp8, values(h, key = element))
# }
# plot_atp8 <- plot_codon_heatmap(atp8_usage,my_palette,color_codes_atp8, codon_color_codes[2:59],"column")
# 
# cob_usage <- make_gene_codondataset(group_gene_usage, "cob")
# acession_per_cob <- aggregate(Taxonomy_Rank ~ Accession, cob_usage, unique)
# cob_usage[["Accession"]] <- NULL
# cob_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_cob <- c()
# for (element in acession_per_cob[,2]) {
#     color_codes_cob <- c(color_codes_cob, values(h, key = element))
# }
# plot_cob <- plot_codon_heatmap(cob_usage,my_palette,color_codes_cob, codon_color_codes[2:59],"column")
# 
# 
# cox1_usage <- make_gene_codondataset(group_gene_usage, "cox1")
# acession_per_cox1 <- aggregate(Taxonomy_Rank ~ Accession, cox1_usage, unique)
# cox1_usage[["Accession"]] <- NULL
# cox1_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_cox1 <- c()
# for (element in acession_per_cox1[,2]) {
#     color_codes_cox1 <- c(color_codes_cox1, values(h, key = element))
# }
# plot_cox1 <- plot_codon_heatmap(cox1_usage,my_palette,color_codes_cox1, codon_color_codes[2:59],"column")
# 
# 
# cox2_usage <- make_gene_codondataset(group_gene_usage, "cox2")
# acession_per_cox2 <- aggregate(Taxonomy_Rank ~ Accession, cox2_usage, unique)
# cox2_usage[["Accession"]] <- NULL
# cox2_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_cox2 <- c()
# for (element in acession_per_cox2[,2]) {
#     color_codes_cox2 <- c(color_codes_cox2, values(h, key = element))
# }
# plot_cox2 <- plot_codon_heatmap(cox2_usage,my_palette,color_codes_cox2, codon_color_codes[2:59],"column")
# 
# 
# cox3_usage <- make_gene_codondataset(group_gene_usage, "cox3")
# acession_per_cox3 <- aggregate(Taxonomy_Rank ~ Accession, cox3_usage, unique)
# cox3_usage[["Accession"]] <- NULL
# cox3_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_cox3 <- c()
# for (element in acession_per_cox3[,2]) {
#     color_codes_cox3 <- c(color_codes_cox3, values(h, key = element))
# }
# plot_cox3 <- plot_codon_heatmap(cox3_usage,my_palette,color_codes_cox3, codon_color_codes[2:59],"column")
# 
# 
# nad1_usage <- make_gene_codondataset(group_gene_usage, "nad1")
# acession_per_nad1 <- aggregate(Taxonomy_Rank ~ Accession, nad1_usage, unique)
# nad1_usage[["Accession"]] <- NULL
# nad1_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_nad1 <- c()
# for (element in acession_per_nad1[,2]) {
#     color_codes_nad1 <- c(color_codes_nad1, values(h, key = element))
# }
# plot_nad1 <- plot_codon_heatmap(nad1_usage,my_palette,color_codes_nad1, codon_color_codes[2:59],"column")
# 
# 
# nad2_usage <- make_gene_codondataset(group_gene_usage, "nad2")
# acession_per_nad2 <- aggregate(Taxonomy_Rank ~ Accession, nad2_usage, unique)
# nad2_usage[["Accession"]] <- NULL
# nad2_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_nad2 <- c()
# for (element in acession_per_nad2[,2]) {
#     color_codes_nad2 <- c(color_codes_nad2, values(h, key = element))
# }
# plot_nad2 <- plot_codon_heatmap(nad2_usage,my_palette,color_codes_nad2, codon_color_codes[2:59],"column")
# 
# 
# nad3_usage <- make_gene_codondataset(group_gene_usage, "nad3")
# acession_per_nad3 <- aggregate(Taxonomy_Rank ~ Accession, nad3_usage, unique)
# nad3_usage[["Accession"]] <- NULL
# nad3_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_nad3 <- c()
# for (element in acession_per_nad3[,2]) {
#     color_codes_nad3 <- c(color_codes_nad3, values(h, key = element))
# }
# plot_nad3 <- plot_codon_heatmap(nad3_usage,my_palette,color_codes_nad3, codon_color_codes[2:59],"column")
# 
# 
# nad4_usage <- make_gene_codondataset(group_gene_usage, "nad4")
# acession_per_nad4 <- aggregate(Taxonomy_Rank ~ Accession, nad4_usage, unique)
# nad4_usage[["Accession"]] <- NULL
# nad4_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_nad4 <- c()
# for (element in acession_per_nad4[,2]) {
#     color_codes_nad4 <- c(color_codes_nad4, values(h, key = element))
# }
# plot_nad4 <- plot_codon_heatmap(nad4_usage,my_palette,color_codes_nad4, codon_color_codes[2:59],"column")
# 
# 
# nad4l_usage <- make_gene_codondataset(group_gene_usage, "nad4l")
# acession_per_nad4l <- aggregate(Taxonomy_Rank ~ Accession, nad4l_usage, unique)
# nad4l_usage[["Accession"]] <- NULL
# nad4l_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_nad4l <- c()
# for (element in acession_per_nad4l[,2]) {
#     color_codes_nad4l <- c(color_codes_nad4l, values(h, key = element))
# }
# plot_nad4l <- plot_codon_heatmap(nad4l_usage,my_palette,color_codes_nad4l, codon_color_codes[2:59],"column")
# 
# 
# nad5_usage <- make_gene_codondataset(group_gene_usage, "nad5")
# acession_per_nad5 <- aggregate(Taxonomy_Rank ~ Accession, nad5_usage, unique)
# nad5_usage[["Accession"]] <- NULL
# nad5_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_nad5 <- c()
# for (element in acession_per_nad5[,2]) {
#     color_codes_nad5 <- c(color_codes_nad5, values(h, key = element))
# }
# plot_nad5 <- plot_codon_heatmap(nad5_usage,my_palette,color_codes_nad5, codon_color_codes[2:59],"column")
# 
# 
# nad6_usage <- make_gene_codondataset(group_gene_usage, "nad6")
# acession_per_nad6 <- aggregate(Taxonomy_Rank ~ Accession, nad6_usage, unique)
# nad6_usage[["Accession"]] <- NULL
# nad6_usage[["Taxonomy_Rank"]] <- NULL
# color_codes_nad6 <- c()
# for (element in acession_per_nad6[,2]) {
#     color_codes_nad6 <- c(color_codes_nad6, values(h, key = element))
# }
# plot_nad6 <- plot_codon_heatmap(nad6_usage,my_palette,color_codes_nad6, codon_color_codes[2:59],"column")
# 
# 
# # plots
# 
# plot_all_codon <- plot_codon_heatmap(all_codon_usage,my_palette,color_codes_all, codon_color_codes,"column")
# plot_atp6 <- plot_codon_heatmap(atp6_usage,my_palette,color_codes_atp6, codon_color_codes[2:59],"column")
# plot_atp8 <- plot_codon_heatmap(atp8_usage,my_palette,color_codes_atp8, codon_color_codes[2:59],"column") ####
# plot_cob <- plot_codon_heatmap(cob_usage,my_palette,color_codes_cob, codon_color_codes[2:59],"column")
# plot_cox1 <- plot_codon_heatmap(cox1_usage,my_palette,color_codes_cox1, codon_color_codes[2:59],"column")
# plot_cox2 <- plot_codon_heatmap(cox2_usage,my_palette,color_codes_cox2, codon_color_codes[2:59],"column")
# plot_cox3 <- plot_codon_heatmap(cox3_usage,my_palette,color_codes_cox3, codon_color_codes[2:59],"column")
# plot_nad1 <- plot_codon_heatmap(nad1_usage,my_palette,color_codes_nad1, codon_color_codes[2:59],"column")
# plot_nad2 <- plot_codon_heatmap(nad2_usage,my_palette,color_codes_nad2, codon_color_codes[2:59],"column")
# plot_nad3 <- plot_codon_heatmap(nad3_usage,my_palette,color_codes_nad3, codon_color_codes[2:59],"column") ####
# plot_nad4 <- plot_codon_heatmap(nad4_usage,my_palette,color_codes_nad4, codon_color_codes[2:59],"column")
# plot_nad4l <- plot_codon_heatmap(nad4l_usage,my_palette,color_codes_nad4l, codon_color_codes[2:59],"column") ####
# plot_nad5 <- plot_codon_heatmap(nad5_usage,my_palette,color_codes_nad5, codon_color_codes[2:59],"column")
# plot_nad6 <- plot_codon_heatmap(nad6_usage,my_palette,color_codes_nad6, codon_color_codes[2:59],"column")  #####



