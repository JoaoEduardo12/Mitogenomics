### Script to analyse mitogenomes ###

# Separar pelo mesmo grupo taxonomico;
# Fazer testes de seleção

##### Functions ######

get_taxonomy_column <- function(tax_vector, dataset, column, new_col) {
    # "tax_vector is a vector that contains the taxonomic ranks to isolate
    # "dataset" is the dataset used to find taxonomic ranks
    # Given a vector of taxonomic ranks, it creates a new column in "dataset" with just
    # one specified taxonomic rank found in "tax_vector"
    dataset[, column] <- ''
    for(row in 1:nrow(dataset)) {
        for (element in tax_vector) {  #element aqui fica numa sequencia e nao num ponto
            if (element %in% str_split(dataset[row,5], ';')[[1]]) {
                dataset[row,new_col] <- element
            }
        }
    }
    dataset
}

create_dataframe <- function(my_dataset,vector,column, my_column) {
    taxonomy_names <- c()
    mean_values <- c()
    sd_values <- c()
    for (element in vector) {
        new_dataset <- subset(my_dataset, my_dataset[, column] == element)
        taxonomy_names <- c(taxonomy_names, element)
        mean_values <- c(mean_values, mean(new_dataset[[my_column]]))
        sd_values <- c(sd_values, sd(new_dataset[[my_column]]))
        df <- data.frame(taxonomy_names,mean_values,sd_values)
    }
    df
}

create_gene_dataframe <- function(my_dataset,vector, column, gene_vector, my_column) {
    taxonomy_names <- c()
    gene_names <- c()
    mean_values <- c()
    sd_values <- c()
    for (element in vector) {
        new_dataset <- subset(my_dataset, my_dataset[, column] == element)
        for (gene in gene_vector) {
            gene_dataset <- subset(new_dataset, new_dataset$Content == gene)
            taxonomy_names <- c(taxonomy_names, element)
            gene_names <- c(gene_names, gene)
            mean_values <- c(mean_values, mean(gene_dataset[[my_column]]))
            sd_values <- c(sd_values, sd(gene_dataset[[my_column]]))
        }
    }
    df <- data.frame(taxonomy_names,gene_names,mean_values,sd_values)
    df
}


##### Data Processing ####
# Set the right working directory
setwd("~/Desktop/Bioinformatica/Chondrichthyes/chondrichthyes")

# Import necessary library to open csv file
library(readr)

# Read the csv file
Final_Results <- read_csv("Final_Results.csv")

# Load tidyverse packages
library(plyr)
library(stringr)

# Vector to use in function to isolate specific taxonomic ranks
taxonomy_vector_1 <- c('Elasmobranchii', 'Holocephali')
taxonomy_vector_2 <- c('Galeomorphii','Squalomorphii','Batoidea','Chimaeridae','Rhinochimaeridae','Callorhinchidae')
taxonomy_vector_3 <- c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes','Batoidea incertae sedis')
taxonomy_vector_4 <- c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')
taxonomy_vector_5 <- c('Batoidea','Galeomorphii','Squalomorphii','Holocephali')

# Get an extra column in the dataset with values of specified taxonomic ranks
Final_Results <- get_taxonomy_column(taxonomy_vector_3,Final_Results, 'Taxonomy_Rank', 28)
Final_Results <- get_taxonomy_column(taxonomy_vector_5,Final_Results, 'Upper_Rank', 29)
Final_Results <- get_taxonomy_column(taxonomy_vector_4,Final_Results, 'Family', 30)


# Isolating complete mitogenomes
complete_mitogenomes <- subset(Final_Results, Final_Results$Content == 'Mitogenome')

# Creating a vector to isolate protein coding genes 
protein_coding <- c('atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6')

# Isolating all protein coding genes
PCG <- subset(Final_Results, Final_Results$Content %in% protein_coding)


##### GC skew vs AT skew #####

# load ggplot
library(ggplot2)

# Scatter plot com GC content a cores
# ggplot(data=complete_mitogenomes, aes(x = GC_skew, y = AT_skew)) +
#     xlab("GC skew") +  ylab("AT skew") +
#     geom_point(aes(color = GC_content, shape = factor(Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Echinorhiniformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes'))), size = 2.5, alpha = I(2/2)) +
#     geom_vline(aes(xintercept = mean(GC_skew)), color = "red", linetype = "dashed") +
#     geom_hline(aes(yintercept = mean(AT_skew)), color = "red", linetype = "dashed") +
#     scale_color_gradient(low = "yellow", high = "red") +
#     labs(color = 'GC content', shape = 'Taxonomy') +
#     xlim(c(-0.5,-0.1)) +
#     ggtitle('Compositional complete mitogenome differences') #+
#     #scale_shape_manual(values=c(19,17,15,23,3,10))

full_mitogenomes <- subset(complete_mitogenomes, complete_mitogenomes$State != 'partial')

# Scatter plot simples
ggplot(data=subset(complete_mitogenomes, complete_mitogenomes$Taxonomy_Rank != 'Batoidea incertae sedis'), aes(x = GC_skew, y = AT_skew)) +
    xlab("GC skew") +  ylab("AT skew") +
    geom_point(aes(color = factor(Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes',
                                                            'Lamniformes','Orectolobiformes','Heterodontiformes','Echinorhiniformes','Hexanchiformes',
                                                            'Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes'))
                   , shape = factor(Upper_Rank, levels = c('Batoidea','Galeomorphii','Squalomorphii','Holocephali'))), size = 2.5, alpha = I(2/2)) +
    scale_color_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                               '#006400','#7b68ee')) +
    geom_vline(aes(xintercept = mean(GC_skew)), color = "grey", linetype = "dashed") +
    geom_hline(aes(yintercept = mean(AT_skew)), color = "grey", linetype = "dashed") +
    xlim(c(-0.5,-0.1)) +
    labs(color = 'Taxonomy', shape = 'Group') +
    scale_shape_manual(values=c(17,3,10,15,15)) +
    ggtitle('Compositional complete mitogenome differences')

# Scatter plot sem parciais


ggplot(data=full_mitogenomes, aes(x = GC_skew, y = AT_skew)) +
    xlab("GC skew") +  ylab("AT skew") +
    geom_point(aes(color = factor(Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes',
                                                            'Lamniformes','Orectolobiformes','Heterodontiformes','Echinorhiniformes','Hexanchiformes',
                                                            'Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes'))
                   , shape = factor(Upper_Rank, levels = c('Batoidea','Galeomorphii','Squalomorphii','Holocephali'))), size = 2.5, alpha = I(2/2)) +
    scale_color_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                                '#006400','#7b68ee')) +
    geom_vline(aes(xintercept = mean(GC_skew)), color = "grey", linetype = "dashed") +
    geom_hline(aes(yintercept = mean(AT_skew)), color = "grey", linetype = "dashed") +
    xlim(c(-0.5,-0.1)) +
    labs(color = 'Taxonomy', shape = 'Group') +
    scale_shape_manual(values=c(17,3,10,15,15)) +
    ggtitle('Compositional complete mitogenome differences')


###### Mitogenome Length and Gene length comparison between taxonomic ranks #######

# Isolate the actual complete mitogenomes (those above 13000, thus excluding unverified sequences that would cause noise)
#mitogenome_len <- create_dataframe(complete_mitogenomes, taxonomy_vector_3, "Length")

# Barplot for mitogenome length
# ggplot(mitogenome_len[1:13,], aes(x=factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), y=mean_values, fill = factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')))) +
#     geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
#     ggtitle('Mitogenome Length') +
#     scale_fill_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
#                                 '#006400','#7b68ee')) +
#     labs(x = 'Taxonomic ranks', y = 'Mitogenome Length') +
#     theme(text = element_text(size=10),
#           axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
#     geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.1,
#                   position=position_dodge(.9))

mitogenome_len <- create_dataframe(full_mitogenomes, taxonomy_vector_3, 28, "Length")
mitogenome_len2 <- create_dataframe(full_mitogenomes, taxonomy_vector_4, 30, "Length")


# Sem parciais
ggplot(mitogenome_len[1:13,], aes(x=factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), y=mean_values, fill = factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')))) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    ggtitle('Mitogenome Length') +
    scale_fill_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                               '#006400','#7b68ee')) +
    labs(x = 'Taxonomic ranks', y = 'Mitogenome Length') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.1,
                  position=position_dodge(.9))

ggplot(mitogenome_len2, aes(x=factor(taxonomy_names, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')), y=mean_values, fill = factor(taxonomy_names, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')))) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .3) +
    ggtitle('Chimaeriformes') +
    #theme_classic() +
    scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +
    labs(x = 'Family', y = 'Mitogenome Length') +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.1,
                  position=position_dodge(.9))

# Isolating gene length
full_PCG <- subset(PCG, PCG$State != 'partial')

gene_len <- create_gene_dataframe(PCG, taxonomy_vector_3, protein_coding, 'Length')

# Plotting gene length
# ggplot(gene_len[1:169,], aes(x=factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), y=mean_values, fill = factor(gene_names, levels = c('atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6')))) +
#     geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
#     ggtitle('Gene Length') +
#     labs(x = '', y = 'Gene Length', fill = 'Gene') +
#     theme(text = element_text(size=12),
#           axis.text.x = element_text(angle=0, hjust=0.5)) +
#     geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.2,
#                   position=position_dodge(.7))

# Plotting gene length - more intuitive graph
# ggplot(gene_len[1:169,], aes(x = factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), y=mean_values, group = factor(gene_names, levels = c('nad5','cox1','nad4','cob','nad2','nad1','cox3','cox2','atp6','nad6','nad3','nad4l','atp8')), color=factor(gene_names, levels = c('nad5','cox1','nad4','cob','nad2','nad1','cox3','cox2','atp6','nad6','nad3','nad4l','atp8')))) +
#     geom_line(size = 1.1) +
#     geom_point(size = 3.2) +
#     #scale_fill_manual(values=c("#fcfc4b", "#3ace3a", "#80CEE1","#ff392e")) +
#     ggtitle('Gene Length') +
#     scale_color_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
#                                 '#006400','#7b68ee')) +
#     labs(x = '', y = 'Size', color= 'Gene') +
#     theme(text = element_text(size=15),
#           axis.text.x = element_text(angle=50, hjust=1)) +
#     geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.7,
#                   position=position_dodge(0))

gene_len <- create_gene_dataframe(full_PCG, taxonomy_vector_3, protein_coding, 'Length')

library(grid)
p <- ggplot(gene_len[1:169,], aes(x = factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes',
                                                                        'Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes',
                                                                        'Chimaeriformes')), 
                                  y=mean_values, 
                                  group = factor(gene_names, levels = c('nad5','cox1','nad4','cob','nad2','nad1','cox3','cox2','atp6','nad6','nad3','nad4l','atp8')), 
                                  color=factor(gene_names, levels = c('nad5','cox1','nad4','cob','nad2','nad1','cox3','cox2','atp6','nad6','nad3','nad4l','atp8')))) +
    scale_colour_discrete(guide = 'none')  +
    geom_line(size = 1.1) +
    geom_point(size = 3.2) +
    #geom_text(data = subset(gene_len[1:169,], taxonomy_names == "Chimaeriformes"), aes(label = gene_names, x = Inf, y = mean_values), hjust = -.1) +
    #scale_fill_manual(values=c("#fcfc4b", "#3ace3a", "#80CEE1","#ff392e")) +
    ggtitle('Gene Length') +
    theme(plot.margin = unit(c(1,3,1,1), "lines")) +
    #theme_classic() +
    scale_color_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                                '#006400','#7b68ee')) +
    labs(x = '', y = 'Size', color= 'Gene') +
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle=50, hjust=1)) +
    geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.7,
                  position=position_dodge(0))

p

gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

##### Mitogenome and Gene GC content ######

# Plotting mitogenome gc content
ggplot(full_mitogenomes, aes(x = factor(Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), y = `GC content`, color = factor(Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')))) +
    geom_jitter(aes(x = factor(Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), color = factor(Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes'))), 
                position = position_jitter(width = .05), alpha = 0.5) +
    geom_boxplot(outlier.colour = NA, position = "dodge") +
    ggtitle('Mitogenome GC content') +
    scale_color_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                                '#006400','#7b68ee')) +
    labs(x = '', y = 'GC content(%)') +
    theme(text = element_text(size=12),
          axis.text.x = element_text(angle=0, hjust=0.5),
          legend.position="none") #+
    #stat_summary(fun=mean, geom="point", shape=23, size=4)

ggplot(subset(full_mitogenomes, full_mitogenomes$Taxonomy_Rank == 'Chimaeriformes'), aes(x = factor(Family, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')), y = `GC content`, color = factor(Family, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')))) +
    geom_jitter(aes(x = factor(Family, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')), color = factor(Family, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae'))), 
                position = position_jitter(width = .05), alpha = 0.5) +
    geom_boxplot(position = "dodge", width = .3) +
    ggtitle('Chimaeriformes') +
    scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +
    labs(x = '', y = 'GC content(%)') +
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle=0, hjust=0.5),
          legend.position="none") #+
    #stat_summary(fun=mean, geom="point", shape=23, size=4)

gene_gc <- create_gene_dataframe(full_PCG, taxonomy_vector_3, protein_coding, 'GC content')
gene_gc2 <- create_gene_dataframe(PCG, taxonomy_vector_4, protein_coding, 'GC content')

ggplot(gene_gc[1:169,], aes(x = factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes',
                                                                   'Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes',
                                                                   'Chimaeriformes')), 
                             y=mean_values, 
                             group = factor(gene_names, levels = c('nad5','cox1','nad4','cob','nad2','nad1','cox3','cox2','atp6','nad6','nad3','nad4l','atp8')), 
                             color=factor(gene_names, levels = c('nad5','cox1','nad4','cob','nad2','nad1','cox3','cox2','atp6','nad6','nad3','nad4l','atp8')))) +
    geom_line(size = 0.5) +
    geom_point(size = 1.2) +
    #geom_text(data = subset(gene_len[1:169,], taxonomy_names == "Chimaeriformes"), aes(label = gene_names, x = Inf, y = mean_values), hjust = -.1) +
    #scale_fill_manual(values=c("#fcfc4b", "#3ace3a", "#80CEE1","#ff392e")) +
    ggtitle('Gene Length') +
    theme(plot.margin = unit(c(1,3,1,1), "lines")) +
    #theme_classic() +
    scale_color_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                                '#006400','#7b68ee')) +
    labs(x = '', y = 'Size', color= 'Gene') +
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle=50, hjust=1)) +
    geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.7,
                  position=position_dodge(0))

ggplot(gene_gc2, aes(x = factor(gene_names), 
                            y=mean_values, 
                            group = factor(taxonomy_names, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')), 
                            color=factor(taxonomy_names, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')))) +
    geom_line(size = 1.2) +
    geom_point(size = 3.2) +
    #geom_text(data = subset(gene_len[1:169,], taxonomy_names == "Chimaeriformes"), aes(label = gene_names, x = Inf, y = mean_values), hjust = -.1) +
    #scale_fill_manual(values=c("#fcfc4b", "#3ace3a", "#80CEE1","#ff392e")) +
    ggtitle('Gene GC content') +
    theme(plot.margin = unit(c(1,3,1,1), "lines")) +
    #theme_classic() +
    scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +
    labs(x = 'protein-coding genes', y = '%GC', color= 'Family') +
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle=50, hjust=1)) +
    geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.7,
                  position=position_dodge(0))

mitogenome_gc <- create_dataframe(full_mitogenomes, taxonomy_vector_3, "GC content")


# Barplots for each taxonomic rank
# ggplot(gene_gc[1:169,], aes(x=factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), y=mean_values, fill = factor(gene_names, levels = c('atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6')))) +
#     geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .9) +
#     ggtitle('') +
#     labs(x = '', y = 'Gene GC content', fill = 'Gene') +
#     theme(text = element_text(size=8),
#           axis.text.x = element_text(angle=0, hjust=0.5)) +
#     geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.2,
#                   position=position_dodge(.9))


# Overall PCG GC content
# ggplot(full_PCG, aes(x=factor(Content, levels = c('atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6')), y=`GC content`, color = factor(Content, levels = c('atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6')), width = 2)) +
#     #geom_jitter(aes(x = factor(Content, levels = c('atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6')), color = factor(Content, levels = c('atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6'))), 
#      #           position = position_jitter(width = .05), alpha = 0.5) +
#     geom_boxplot(outlier.colour = NA, position = "dodge", size = .9) +
#     ggtitle('Gene GC content') +
#     labs(x = 'Gene', y = 'GC content (%)') +
#     theme(text = element_text(size=8),
#           axis.text.x = element_text(angle=0, hjust=0.5),
#           legend.position="none") +
#     stat_summary(fun=mean, geom="point", shape=23, size=4)


##### Start and Stop Codons #####

# Load tidyverse packages
library(dplyr)
library(tidyr)

# Preprocess to obtain counts and percentages of two categorical variables
start_perc <- PCG %>% 
    group_by(Content, Start) %>% 
    tally() %>% 
    complete(Start, fill = list(n = 0)) %>% 
    mutate(percentage = n / sum(n) * 100)


# Start Codon Frequency
ggplot(start_perc, aes(x=Content, y = percentage , fill = factor(Start))) +
    geom_bar(stat="identity",position='stack') +
    ggtitle('Start Codon Frequency per Gene') +
    labs(x = 'Genes', y = 'Frequency (%)',fill = 'Start Codons') +
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle=0, hjust=0.5))

# ditto
stop_perc <- PCG %>% 
    group_by(Content, Stop) %>% 
    tally() %>% 
    complete(Stop, fill = list(n = 0)) %>% 
    mutate(percentage = n / sum(n) * 100)

# Stop Codon Frequency
ggplot(stop_perc, aes(x=Content, y = percentage , fill = factor(Stop))) +
    geom_bar(stat="identity",position='stack', color = '#505050') +
    ggtitle('Stop Codon Frequency per Gene') +
    labs(x = 'Genes', y = 'Frequency (%)',fill = 'Stop Codons') +
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle=0, hjust=0.5))

# Get taxonomy rank also grouped
start_perc <- PCG %>% 
    group_by(Taxonomy_Rank, Content, Start) %>% 
    tally() %>% 
    complete(Start, fill = list(n = 0)) %>% 
    mutate(percentage = n / sum(n) * 100)

# Start Codon Frequency seperated by taxonomic rank
ggplot(start_perc, aes(x=interaction(Taxonomy_Rank,Content, lex.order = TRUE), y = percentage , group = 1, fill = factor(Start))) +
    geom_bar(stat="identity",position='stack', color = '#505050') +
    #annotate(geom = "text", x = seq_len(nrow(start_perc)), y = 34, label = start_perc$Content, size = 4) +
    #annotate(geom = "text", x = 0.5 + 4 * (0:5), y = 32, label = unique(start_perc$Taxonomy_Rank), size = 6) +
    #coord_cartesian(ylim = c(35, 65), expand = FALSE, clip = "off") +
    ggtitle('Start Codon Frequency per Gene') +
    labs(x = 'Genes', y = 'Frequency (%)',fill = 'Start Codons') +
    theme(text = element_text(size=10),
          #axis.text.x = element_text(angle=0, hjust=0.5),
          #plot.margin = unit(c(1, 1, 4, 1), "lines"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90, hjust=0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())

# Ditto

stop_perc <- PCG %>% 
    group_by(Taxonomy_Rank, Content, Stop) %>% 
    tally() %>% 
    complete(Stop, fill = list(n = 0)) %>% 
    mutate(percentage = n / sum(n) * 100)

# Stop Codon Frequency separated by taxonomic rank

ggplot(stop_perc, aes(x=interaction(Taxonomy_Rank,Content, lex.order = TRUE), y = percentage , group = 1, fill = factor(Stop))) +
    geom_bar(stat="identity",position='stack', color = '#505050') +
    #annotate(geom = "text", x = seq_len(nrow(start_perc)), y = 34, label = start_perc$Content, size = 4) +
    #annotate(geom = "text", x = 0.5 + 4 * (0:5), y = 32, label = unique(start_perc$Taxonomy_Rank), size = 6) +
    #coord_cartesian(ylim = c(35, 65), expand = FALSE, clip = "off") +
    ggtitle('Start Codon Frequency per Gene') +
    labs(x = 'Genes', y = 'Frequency (%)',fill = 'Start Codons') +
    theme(text = element_text(size=10),
          #axis.text.x = element_text(angle=0, hjust=0.5),
          #plot.margin = unit(c(1, 1, 4, 1), "lines"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90, hjust=0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())


##### PCA with nucleotide base #####

make_comp_def <- function(main_df, col1, col2, col3 = 'Length') {
    for (id in unique(main_df[,col1])) {
        for (lines in subset(main_df, main_df[,col1] == id)) {
            
        }
    }
}

new <- full_PCG %>%
    group_by(Accession) %>%
    mutate(G1 = weighted.mean(`%G1`, Length)) %>%
    mutate(C1 = weighted.mean(`%C1`, Length)) %>%
    mutate(A1 = weighted.mean(`%A1`, Length)) %>%
    mutate(T1 = weighted.mean(`%T1`, Length)) %>%
    mutate(G2 = weighted.mean(`%G2`, Length)) %>%
    mutate(C2 = weighted.mean(`%C2`, Length)) %>%
    mutate(A2 = weighted.mean(`%A2`, Length)) %>%
    mutate(T2 = weighted.mean(`%T2`, Length)) %>%
    mutate(G3 = weighted.mean(`%G3`, Length)) %>%
    mutate(C3 = weighted.mean(`%C3`, Length)) %>%
    mutate(A3 = weighted.mean(`%A3`, Length)) %>%
    mutate(T3 = weighted.mean(`%T3`, Length))
    
nucleotide_vector <- aggregate(G1 ~ Accession, new, unique)
c1 <- aggregate(C1 ~ Accession, new, unique)
nucleotide_vector$C1 <- cbind(c1$C1)
a1 <- aggregate(A1 ~ Accession, new, unique)
nucleotide_vector$A1 <- cbind(a1$A1)
t1 <- aggregate(T1 ~ Accession, new, unique)
nucleotide_vector$T1 <- cbind(t1$T1)
g2 <- aggregate(G2 ~ Accession, new, unique)
nucleotide_vector$G2 <- cbind(g2$G2)
c2 <- aggregate(C2 ~ Accession, new, unique)
nucleotide_vector$C2 <- cbind(c2$C2)
a2 <- aggregate(A2 ~ Accession, new, unique)
nucleotide_vector$A2 <- cbind(a2$A2)
t2 <- aggregate(T2 ~ Accession, new, unique)
nucleotide_vector$T2 <- cbind(t2$T2)
g3 <- aggregate(G3 ~ Accession, new, unique)
nucleotide_vector$G3 <- cbind(g3$G3)
c3 <- aggregate(C3 ~ Accession, new, unique)
nucleotide_vector$C3 <- cbind(c3$C3)
a3 <- aggregate(A3 ~ Accession, new, unique)
nucleotide_vector$A3 <- cbind(a3$A3)
t3 <- aggregate(T3 ~ Accession, new, unique)
nucleotide_vector$T3 <- cbind(t3$T3)

# aggregate taxonomic info
tax <- aggregate(Taxonomy_Rank ~ Accession, new, unique)
nucleotide_vector$Taxonomy_Rank <- cbind(tax$Taxonomy_Rank)

# concatenate PCGs
# nucleotide_vector <- aggregate(`%G1` ~ Accession, PCG, mean)
# c1 <- aggregate(`%C1` ~ Accession, PCG, mean)
# nucleotide_vector$`%C1` <- cbind(c1$`%C1`)
# a1 <- aggregate(`%A1` ~ Accession, PCG, mean)
# nucleotide_vector$`%A1` <- cbind(a1$`%A1`)
# t1 <- aggregate(`%T1` ~ Accession, PCG, mean)
# nucleotide_vector$`%T1` <- cbind(t1$`%T1`)
# g2 <- aggregate(`%G2` ~ Accession, PCG, mean)
# nucleotide_vector$`%G2` <- cbind(g2$`%G2`)
# c2 <- aggregate(`%C2` ~ Accession, PCG, mean)
# nucleotide_vector$`%C2` <- cbind(c2$`%C2`)
# a2 <- aggregate(`%A2` ~ Accession, PCG, mean)
# nucleotide_vector$`%A2` <- cbind(a2$`%A2`)
# t2 <- aggregate(`%T2` ~ Accession, PCG, mean)
# nucleotide_vector$`%T2` <- cbind(t2$`%T2`)
# g3 <- aggregate(`%G3` ~ Accession, PCG, mean)
# nucleotide_vector$`%G3` <- cbind(g3$`%G3`)
# c3 <- aggregate(`%C3` ~ Accession, PCG, mean)
# nucleotide_vector$`%C3` <- cbind(c3$`%C3`)
# a3 <- aggregate(`%A3` ~ Accession, PCG, mean)
# nucleotide_vector$`%A3` <- cbind(a3$`%A3`)
# t3 <- aggregate(`%T3` ~ Accession, PCG, mean)
# nucleotide_vector$`%T3` <- cbind(t3$`%T3`)

# aggregate taxonomic info
#tax <- aggregate(Taxonomy_Rank ~ Accession, PCG, unique)
#nucleotide_vector$Taxonomy_Rank <- cbind(tax$Taxonomy_Rank)

# PCA with prcomp
nucleotide_vector2 <- subset(nucleotide_vector, nucleotide_vector$Taxonomy_Rank != 'Batoidea incertae sedis')
pc <- prcomp(nucleotide_vector2[,2:13], center = TRUE, scale = TRUE)
summary(pc)


library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)

ggbiplot(pc, ellipse = TRUE, obs.scale = 1,  var.axes = TRUE, var.scale = TRUE, 
         groups = factor(nucleotide_vector2$Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes',
                                                                      'Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes',
                                                                      'Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes',
                                                                      'Chimaeriformes'))) +
         scale_color_manual(name = 'Taxonomic Rank', values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32',
                                           '#3cb371','#006400','#7b68ee'))

new <- subset(PCG, PCG$Taxonomy_Rank == 'Chimaeriformes') %>%
    group_by(Accession) %>%
    mutate(G1 = weighted.mean(`%G1`, Length)) %>%
    mutate(C1 = weighted.mean(`%C1`, Length)) %>%
    mutate(A1 = weighted.mean(`%A1`, Length)) %>%
    mutate(T1 = weighted.mean(`%T1`, Length)) %>%
    mutate(G2 = weighted.mean(`%G2`, Length)) %>%
    mutate(C2 = weighted.mean(`%C2`, Length)) %>%
    mutate(A2 = weighted.mean(`%A2`, Length)) %>%
    mutate(T2 = weighted.mean(`%T2`, Length)) %>%
    mutate(G3 = weighted.mean(`%G3`, Length)) %>%
    mutate(C3 = weighted.mean(`%C3`, Length)) %>%
    mutate(A3 = weighted.mean(`%A3`, Length)) %>%
    mutate(T3 = weighted.mean(`%T3`, Length))

nucleotide_vector <- aggregate(G1 ~ Accession, new, unique)
c1 <- aggregate(C1 ~ Accession, new, unique)
nucleotide_vector$C1 <- cbind(c1$C1)
a1 <- aggregate(A1 ~ Accession, new, unique)
nucleotide_vector$A1 <- cbind(a1$A1)
t1 <- aggregate(T1 ~ Accession, new, unique)
nucleotide_vector$T1 <- cbind(t1$T1)
g2 <- aggregate(G2 ~ Accession, new, unique)
nucleotide_vector$G2 <- cbind(g2$G2)
c2 <- aggregate(C2 ~ Accession, new, unique)
nucleotide_vector$C2 <- cbind(c2$C2)
a2 <- aggregate(A2 ~ Accession, new, unique)
nucleotide_vector$A2 <- cbind(a2$A2)
t2 <- aggregate(T2 ~ Accession, new, unique)
nucleotide_vector$T2 <- cbind(t2$T2)
g3 <- aggregate(G3 ~ Accession, new, unique)
nucleotide_vector$G3 <- cbind(g3$G3)
c3 <- aggregate(C3 ~ Accession, new, unique)
nucleotide_vector$C3 <- cbind(c3$C3)
a3 <- aggregate(A3 ~ Accession, new, unique)
nucleotide_vector$A3 <- cbind(a3$A3)
t3 <- aggregate(T3 ~ Accession, new, unique)
nucleotide_vector$T3 <- cbind(t3$T3)

# aggregate taxonomic info
tax <- aggregate(Family ~ Accession, new, unique)
nucleotide_vector$Family <- cbind(tax$Family)

pc <- prcomp(nucleotide_vector[,2:13], center = TRUE, scale = TRUE)

ggbiplot(pc, ellipse = TRUE, obs.scale = 1,  var.axes = FALSE, var.scale = TRUE, 
         groups = factor(nucleotide_vector$Family, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae'))) +
    scale_color_manual(name = 'Family', values=c("#00AFBB", "#E7B800", "#FC4E07"))

# With the FactoMineR

install.packages("FactoMineR")
library(FactoMineR)

nucleotide.pca <- PCA(nucleotide_vector[,2:13])
ggbiplot(nucleotide.pca, ellipse = TRUE,obs.scale = 1, var.scale = 1, var.axes = TRUE, groups = nucleotide_vector$Taxonomy_Rank) +
    theme(legend.position = "bottom")

###### Clustering ######

# Isolating GC_content on whole mitogenomes
new <- full_mitogenomes %>%
    group_by(Taxonomy_Rank) %>%
    mutate(GC_norm = weighted.mean(`GC content`, Length))

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
    mutate(GC_norm = weighted.mean(`GC content`, Length))

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
    mutate(GC_norm = weighted.mean(`GC content`, Length))

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
    mutate(GC_norm = weighted.mean(`GC content`, Length))

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


