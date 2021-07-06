### Script to analyse mitogenomes ###

# import function scripts
source("functions.R")

##### Data Processing ####
# Set the right working directory
setwd("~/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/chondrichthyes")

# Import necessary library to open csv file
library(readr)

# Read the csv file
Final_Results <- read_csv("Final_Results_2.csv")
unique_Final_Results <- read_csv("~/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Uniques/unique_Final_Results.csv")
# 
# Load tidyverse packages
library(plyr)
library(stringr)

# Vector to use in function to isolate specific taxonomic ranks
taxonomy_vector_1 <- c('Elasmobranchii', 'Holocephali') # Subclasses
taxonomy_vector_2 <- c('Galeomorphii','Squalomorphii','Batoidea','Chimaeridae','Rhinochimaeridae','Callorhinchidae') # Superorders
taxonomy_vector_3 <- c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes','Batoidea incertae sedis') # Orders
taxonomy_vector_4 <- c('Callorhinchidae','Chimaeridae','Rhinochimaeridae') # Families
taxonomy_vector_5 <- c('Batoidea','Galeomorphii','Squalomorphii','Holocephali')

# Get extra columns in the dataset with values of specified taxonomic ranks
# In this case I want info on order, superclass, and family of chimaeriformes
Final_Results <- get_taxonomy_column(taxonomy_vector_3, unique_Final_Results, 'Taxonomy_Rank')
Final_Results <- get_taxonomy_column(taxonomy_vector_5, Final_Results, 'Upper_Rank')
Final_Results <- get_taxonomy_column(taxonomy_vector_4, Final_Results, 'Family')


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

# isolate complete mitogenomes excluding partial
full_mitogenomes <- subset(complete_mitogenomes, complete_mitogenomes$State != 'partial')

# Simple scatter plot for the total mitogenomes
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

# Simple scatter plot without partial mitogenomes

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

# Create a particular dataset with length mean and sd that we can use to plot in ggplot
mitogenome_len <- create_meansd_dataframe(full_mitogenomes, taxonomy_vector_3, "Taxonomy_Rank", "Length")
# Also dividing on the families of Chimaeriformes so we can also see differences between distinct families
mitogenome_len_families <- create_meansd_dataframe(full_mitogenomes, taxonomy_vector_4, "Family", "Length")

# Plot mitogenome length for all taxonomy ranks
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

# Plot mitogenome length for each of the families
ggplot(mitogenome_len_families, aes(x=factor(taxonomy_names, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')), y=mean_values, fill = factor(taxonomy_names, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')))) +
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
# Create unique gene data frame with mean and sd values that can later be used in ggplot
gene_len <- create_gene_dataframe(full_PCG, taxonomy_vector_3, 'Taxonomy_Rank', protein_coding, 'Length')

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

# Plotting full mitogenome gc content
ggplot(full_mitogenomes, aes(x = factor(Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), y = GC_content, color = factor(Taxonomy_Rank, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')))) +
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

# Plotting the same gc graph but only for the families in Chimaeriformes
ggplot(subset(full_mitogenomes, full_mitogenomes$Taxonomy_Rank == 'Chimaeriformes'), aes(x = factor(Family, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')), y = GC_content, color = factor(Family, levels = c('Callorhinchidae','Chimaeridae','Rhinochimaeridae')))) +
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

# Creating a dataframe to seperate the genes according to taxonomic rank
#The following line separates by superclass
gene_gc <- create_gene_dataframe(full_PCG, taxonomy_vector_3, 'Taxonomy_Rank', protein_coding, 'GC_content')
# This next one separates by families of Chimaeriformes
gene_gc_families <- create_gene_dataframe(full_PCG, taxonomy_vector_4, 'Family', protein_coding, 'GC_content')

# Gene GC content
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
    ggtitle('Gene GC content') +
    theme(plot.margin = unit(c(1,3,1,1), "lines")) +
    #theme_classic() +
    scale_color_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                                '#006400','#7b68ee')) +
    labs(x = '', y = 'GC content', color= 'Gene') +
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle=50, hjust=1)) +
    geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.7,
                  position=position_dodge(0))

ggplot(gene_gc_families, aes(x = factor(gene_names), 
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

# ditto for stop codons
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

# Ditto for stop codons

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


## Codon Analyses

codon_usage <- read.csv("/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Codon_Analyses.csv")

codon_usage <- get_taxonomy_column(taxonomy_vector_3, codon_usage, 'Taxonomy_Rank')
codon_usage <- get_taxonomy_column(taxonomy_vector_5, codon_usage, 'Upper_Rank')
codon_usage <- get_taxonomy_column(taxonomy_vector_4, codon_usage, 'Family')

group_codon_usage <- codon_usage %>% 
    group_by(Taxonomy_Rank, colnames(codon_usage[7:66]))

myliobatiformes_codon_usage <- create_codon_dataframe(codon_usage, 'Myliobatiformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
rajiformes_codon_usage <- create_codon_dataframe(codon_usage, 'Rajiformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
rhinopristiformes_codon_usage <- create_codon_dataframe(codon_usage, 'Rhinopristiformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
torpediniformes_codon_usage <- create_codon_dataframe(codon_usage, 'Torpediniformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
carcharhiniformes_codon_usage <- create_codon_dataframe(codon_usage, 'Carcharhiniformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
lamniformes_codon_usage <- create_codon_dataframe(codon_usage, 'Lamniformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
orectolobiformes_codon_usage <- create_codon_dataframe(codon_usage, 'Orectolobiformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
heterodontiformes_codon_usage <- create_codon_dataframe(codon_usage, 'Heterodontiformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
hexanchiformes_codon_usage <- create_codon_dataframe(codon_usage, 'Hexanchiformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
pristiophoriformes_codon_usage <- create_codon_dataframe(codon_usage, 'Pristiophoriformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
squaliformes_codon_usage <- create_codon_dataframe(codon_usage, 'Squaliformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
squatiniformes_codon_usage <- create_codon_dataframe(codon_usage, 'Squatiniformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
chimaeriformes_codon_usage <- create_codon_dataframe(codon_usage, 'Chimaeriformes', 'Taxonomy_Rank',7:66,'aminoacids',"2")
atp6_codon_usage <- create_codon_dataframe(codon_usage, 'atp6', 'Gene',7:66,'aminoacids',"2")
atp8_codon_usage <- create_codon_dataframe(codon_usage, 'atp8', 'Gene',7:66,'aminoacids',"2")
cox1_codon_usage <- create_codon_dataframe(codon_usage, 'cox1', 'Gene',7:66,'aminoacids',"2")
cox2_codon_usage <- create_codon_dataframe(codon_usage, 'cox2', 'Gene',7:66,'aminoacids',"2")
cox3_codon_usage <- create_codon_dataframe(codon_usage, 'cox3', 'Gene',7:66,'aminoacids',"2")
nad1_codon_usage <- create_codon_dataframe(codon_usage, 'nad1', 'Gene',7:66,'aminoacids',"2")
nad2_codon_usage <- create_codon_dataframe(codon_usage, 'nad2', 'Gene',7:66,'aminoacids',"2")
nad3_codon_usage <- create_codon_dataframe(codon_usage, 'nad3', 'Gene',7:66,'aminoacids',"2")
nad4_codon_usage <- create_codon_dataframe(codon_usage, 'nad4', 'Gene',7:66,'aminoacids',"2")
nad5_codon_usage <- create_codon_dataframe(codon_usage, 'nad5', 'Gene',7:66,'aminoacids',"2")
nad6_codon_usage <- create_codon_dataframe(codon_usage, 'nad6', 'Gene',7:66,'aminoacids',"2")
nad4l_codon_usage <- create_codon_dataframe(codon_usage, 'nad4l', 'Gene',7:66,'aminoacids',"2")
cob_codon_usage <- create_codon_dataframe(codon_usage, 'cob', 'Gene',7:66,'aminoacids',"2")

# myliobatiformes

myliobatiformes_rscu <- plot_rscu(myliobatiformes_codon_usage, 'Myliobatiformes')
myliobatiformes_rscu
ggsave("myliobatiformes_rscu.png",myliobatiformes_rscu)

# rhinopristiformes

rhinopristiformes_rscu <- plot_rscu(rhinopristiformes_codon_usage, 'Rhinopristiformes')
rhinopristiformes_rscu
ggsave("rhinopristiformes_rscu.png",rhinopristiformes_rscu)

#torpediniformes

torpediniformes_rscu <- plot_rscu(torpediniformes_codon_usage, 'Torpediniformes')
torpediniformes_rscu
ggsave("torpediniformes_rscu.png",torpediniformes_rscu)

# carcharhiniformes 

carcharhiniformes_rscu <- plot_rscu(carcharhiniformes_codon_usage, 'Carcharhiniformes')
carcharhiniformes_rscu
ggsave("carcharhiniformes_rscu.png",carcharhiniformes_rscu)

# lamniformes

lamniformes_rscu <- plot_rscu(lamniformes_codon_usage, 'Lamniformes')
lamniformes_rscu
ggsave("lamniformes_rscu.png",lamniformes_rscu)

# orectolobiformes

orectolobiformes_rscu <- plot_rscu(orectolobiformes_codon_usage, 'Orectolobiformes')
orectolobiformes_rscu
ggsave("orectolobiformes_rscu.png",orectolobiformes_rscu)

# heterodontiformes

heterodontiformes_rscu <- plot_rscu(heterodontiformes_codon_usage, 'Heterodontiformes')
heterodontiformes_rscu
ggsave("heterodontiformes_rscu.png",heterodontiformes_rscu)

# hexanchiformes

hexanchiformes_rscu <- plot_rscu(hexanchiformes_codon_usage, 'Hexanchiformes')
hexanchiformes_rscu
ggsave("hexanchiformes_rscu.png",hexanchiformes_rscu)

# pristiophoriformes

pristiophoriformes_rscu <- plot_rscu(pristiophoriformes_codon_usage, 'Pristiophoriformes')
pristiophoriformes_rscu
ggsave("pristiophoriformes_rscu.png",pristiophoriformes_rscu)

# squaliformes

squaliformes_rscu <- plot_rscu(squaliformes_codon_usage, 'Squaliformes')
squaliformes_rscu
ggsave("squaliformes_rscu.png",squaliformes_rscu)

# squatiniformes

squatiniformes_rscu <- plot_rscu(squatiniformes_codon_usage, 'Squatiniformes')
squatiniformes_rscu
ggsave("squatiniformes_rscu.png",squatiniformes_rscu)

# chimaeriformes

chimaeriformes_rscu <- plot_rscu(chimaeriformes_codon_usage, 'Chimaeriformes')
chimaeriformes_rscu
ggsave("chimaeriformes_rscu.png",chimaeriformes_rscu)

#atp6

atp6_rscu <- plot_rscu(atp6_codon_usage, 'ATP6')
atp6_rscu
ggsave("atp6_rscu.png",atp6_rscu)

atp8_rscu <- plot_rscu(atp8_codon_usage, 'ATP8')
atp8_rscu
ggsave("atp8_rscu.png",atp6_rscu)

cox1_rscu <- plot_rscu(cox1_codon_usage, 'COX1')
cox1_rscu
ggsave("atp8_rscu.png",atp6_rscu)

cox2_rscu <- plot_rscu(cox2_codon_usage, 'COX2')
cox2_rscu
ggsave("atp8_rscu.png",atp6_rscu)

cox3_rscu <- plot_rscu(cox3_codon_usage, 'COX3')
cox3_rscu

cob_rscu <- plot_rscu(cob_codon_usage, 'COB')
cob_rscu

nad1_rscu <- plot_rscu(nad1_codon_usage, 'NAD1')
nad1_rscu

nad2_rscu <- plot_rscu(nad2_codon_usage, 'NAD2')
nad2_rscu

nad3_rscu <- plot_rscu(nad3_codon_usage, 'NAD3')
nad3_rscu

nad4_rscu <- plot_rscu(nad4_codon_usage, 'NAD4')
nad4_rscu

nad4l_rscu <- plot_rscu(nad4l_codon_usage, 'NAD4L')
nad4l_rscu

nad5_rscu <- plot_rscu(nad5_codon_usage, 'NAD5')
nad5_rscu

nad6_rscu <- plot_rscu(nad6_codon_usage, 'NAD6')
nad6_rscu

# Pairwise kolmogorov smirnov tests

class_cu <- c("myliobatiformes_codon_usage", "rajiformes_codon_usage", "rhinopristiformes_codon_usage", "torpediniformes_codon_usage", "carcharhiniformes_codon_usage",
              "lamniformes_codon_usage", "orectolobiformes_codon_usage", "heterodontiformes_codon_usage", "hexanchiformes_codon_usage", "pristiophoriformes_codon_usage",
              "squaliformes_codon_usage", "squatiniformes_codon_usage", "chimaeriformes_codon_usage")

gene_cu <- c("atp6_codon_usage", "atp8_codon_usage", "cob_codon_usage", "cox1_codon_usage", "cox2_codon_usage",
              "cox3_codon_usage", "nad1_codon_usage", "nad2_codon_usage", "nad3_codon_usage", "nad4_codon_usage",
              "nad4l_codon_usage", "nad5_codon_usage", "nad6_codon_usage")

ks_class <- pairwise_statistical_test(class_cu, "mean_values", "ks.test")
ks_gene <- pairwise_statistical_test(gene_cu, "mean_values", "ks.test")


# heatmap - clustering

library(dplyr)

group_codon_usage <- codon_usage %>%
    group_by(Accession, Taxonomy_Rank) %>%
    summarise_at(c(7:66), mean, na.rm = TRUE)

group_gene_usage <- codon_usage %>%
    group_by(Accession, Gene, Taxonomy_Rank) %>%
    summarise_at(c(7:66), mean, na.rm = TRUE)

group_codon_usage <- data.frame(group_codon_usage)
rownames(group_codon_usage) <- group_codon_usage[,1]

all_codon_usage <- group_codon_usage
all_codon_usage$Accession <- NULL
all_codon_usage$Taxonomy_Rank <- NULL
all_codon_usage$Upper_Rank <- NULL

acession_per_taxonomy <- codon_usage %>%
    group_by(Accession, Taxonomy_Rank) %>%
    summarise(Taxonomy_Rank)

acession_per_taxonomy <- aggregate(Taxonomy_Rank ~ Accession, codon_usage, unique)



library("gplots")

some_col_func <- function(n) rev(colorspace::heat_hcl(n, h = c(90,180), c = c(100, 30), l = c(30, 80), power = c(1/5, 1.5)))

par(cex.main=0.05, mar = rep(30,4))

my_palette <- colorRampPalette(c("green","red"))(n = 1000)

library(hash)
h <- hash()
h[c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes',
    'Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes',
    'Chimaeriformes')] <- c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                            '#006400','#7b68ee')

color_codes_all <- c()
for (element in acession_per_taxonomy[,2]) {
    color_codes_all <- c(color_codes_all, values(h, key = element))
}


plot_all_codon <- plot_codon_heatmap(all_codon_usage,my_palette,color_codes_all)

atp6_usage <- make_gene_codondataset(group_gene_usage, "atp6")
acession_per_atp6 <- aggregate(Taxonomy_Rank ~ Accession, atp6_usage, unique)
atp6_usage[["Accession"]] <- NULL
atp6_usage[["Taxonomy_Rank"]] <- NULL
color_codes_atp6 <- c()
for (element in acession_per_atp6[,2]) {
    color_codes_atp6 <- c(color_codes_atp6, values(h, key = element))
}
plot_atp6 <- plot_codon_heatmap(atp6_usage,my_palette,color_codes_atp6)

atp8_usage <- make_gene_codondataset(group_gene_usage, "atp8")
acession_per_atp8 <- aggregate(Taxonomy_Rank ~ Accession, atp8_usage, unique)
atp8_usage[["Accession"]] <- NULL
atp8_usage[["Taxonomy_Rank"]] <- NULL
color_codes_atp8 <- c()
for (element in acession_per_atp8[,2]) {
    color_codes_atp8 <- c(color_codes_atp8, values(h, key = element))
}
plot_atp8 <- plot_codon_heatmap(atp8_usage,my_palette,color_codes_atp8)

cob_usage <- make_gene_codondataset(group_gene_usage, "cob")
acession_per_cob <- aggregate(Taxonomy_Rank ~ Accession, cob_usage, unique)
cob_usage[["Accession"]] <- NULL
cob_usage[["Taxonomy_Rank"]] <- NULL
color_codes_cob <- c()
for (element in acession_per_cob[,2]) {
    color_codes_cob <- c(color_codes_cob, values(h, key = element))
}
plot_cob <- plot_codon_heatmap(cob_usage,my_palette,color_codes_cob)


cox1_usage <- make_gene_codondataset(group_gene_usage, "cox1")
acession_per_cox1 <- aggregate(Taxonomy_Rank ~ Accession, cox1_usage, unique)
cox1_usage[["Accession"]] <- NULL
cox1_usage[["Taxonomy_Rank"]] <- NULL
color_codes_cox1 <- c()
for (element in acession_per_cox1[,2]) {
    color_codes_cox1 <- c(color_codes_cox1, values(h, key = element))
}
plot_cox1 <- plot_codon_heatmap(cox1_usage,my_palette,color_codes_cox1)


cox2_usage <- make_gene_codondataset(group_gene_usage, "cox2")
acession_per_cox2 <- aggregate(Taxonomy_Rank ~ Accession, cox2_usage, unique)
cox2_usage[["Accession"]] <- NULL
cox2_usage[["Taxonomy_Rank"]] <- NULL
color_codes_cox2 <- c()
for (element in acession_per_cox2[,2]) {
    color_codes_cox2 <- c(color_codes_cox2, values(h, key = element))
}
plot_cox2 <- plot_codon_heatmap(cox2_usage,my_palette,color_codes_cox2)


cox3_usage <- make_gene_codondataset(group_gene_usage, "cox3")
acession_per_cox3 <- aggregate(Taxonomy_Rank ~ Accession, cox3_usage, unique)
cox3_usage[["Accession"]] <- NULL
cox3_usage[["Taxonomy_Rank"]] <- NULL
color_codes_cox3 <- c()
for (element in acession_per_cox3[,2]) {
    color_codes_cox3 <- c(color_codes_cox3, values(h, key = element))
}
plot_cox3 <- plot_codon_heatmap(cox3_usage,my_palette,color_codes_cox3)


nad1_usage <- make_gene_codondataset(group_gene_usage, "nad1")
acession_per_nad1 <- aggregate(Taxonomy_Rank ~ Accession, nad1_usage, unique)
nad1_usage[["Accession"]] <- NULL
nad1_usage[["Taxonomy_Rank"]] <- NULL
color_codes_nad1 <- c()
for (element in acession_per_nad1[,2]) {
    color_codes_nad1 <- c(color_codes_nad1, values(h, key = element))
}
plot_nad1 <- plot_codon_heatmap(nad1_usage,my_palette,color_codes_nad1)


nad2_usage <- make_gene_codondataset(group_gene_usage, "nad2")
acession_per_nad2 <- aggregate(Taxonomy_Rank ~ Accession, nad2_usage, unique)
nad2_usage[["Accession"]] <- NULL
nad2_usage[["Taxonomy_Rank"]] <- NULL
color_codes_nad2 <- c()
for (element in acession_per_nad1[,2]) {
    color_codes_nad2 <- c(color_codes_nad2, values(h, key = element))
}
plot_nad2 <- plot_codon_heatmap(nad2_usage,my_palette,color_codes_nad2)


nad3_usage <- make_gene_codondataset(group_gene_usage, "nad3")
acession_per_nad3 <- aggregate(Taxonomy_Rank ~ Accession, nad3_usage, unique)
nad3_usage[["Accession"]] <- NULL
nad3_usage[["Taxonomy_Rank"]] <- NULL
color_codes_nad3 <- c()
for (element in acession_per_nad3[,2]) {
    color_codes_nad3 <- c(color_codes_nad3, values(h, key = element))
}
plot_nad3 <- plot_codon_heatmap(nad3_usage,my_palette,color_codes_nad3)


nad4_usage <- make_gene_codondataset(group_gene_usage, "nad4")
acession_per_nad4 <- aggregate(Taxonomy_Rank ~ Accession, nad4_usage, unique)
nad4_usage[["Accession"]] <- NULL
nad4_usage[["Taxonomy_Rank"]] <- NULL
color_codes_nad4 <- c()
for (element in acession_per_nad4[,2]) {
    color_codes_nad4 <- c(color_codes_nad4, values(h, key = element))
}
plot_nad4 <- plot_codon_heatmap(nad4_usage,my_palette,color_codes_nad4)


nad4l_usage <- make_gene_codondataset(group_gene_usage, "nad4l")
acession_per_nad4l <- aggregate(Taxonomy_Rank ~ Accession, nad4l_usage, unique)
nad4l_usage[["Accession"]] <- NULL
nad4l_usage[["Taxonomy_Rank"]] <- NULL
color_codes_nad4l <- c()
for (element in acession_per_nad4l[,2]) {
    color_codes_nad4l <- c(color_codes_nad4l, values(h, key = element))
}
plot_nad4l <- plot_codon_heatmap(nad4l_usage,my_palette,color_codes_nad4l)


nad5_usage <- make_gene_codondataset(group_gene_usage, "nad5")
acession_per_nad5 <- aggregate(Taxonomy_Rank ~ Accession, nad5_usage, unique)
nad5_usage[["Accession"]] <- NULL
nad5_usage[["Taxonomy_Rank"]] <- NULL
color_codes_nad5 <- c()
for (element in acession_per_nad5[,2]) {
    color_codes_nad5 <- c(color_codes_nad5, values(h, key = element))
}
plot_nad5 <- plot_codon_heatmap(nad5_usage,my_palette,color_codes_nad5)


nad6_usage <- make_gene_codondataset(group_gene_usage, "nad6")
acession_per_nad6 <- aggregate(Taxonomy_Rank ~ Accession, nad6_usage, unique)
nad6_usage[["Accession"]] <- NULL
nad6_usage[["Taxonomy_Rank"]] <- NULL
color_codes_nad6 <- c()
for (element in acession_per_nad6[,2]) {
    color_codes_nad6 <- c(color_codes_nad6, values(h, key = element))
}
plot_nad6 <- plot_codon_heatmap(nad6_usage,my_palette,color_codes_nad6)


# plots

plot_all_codon <- plot_codon_heatmap(all_codon_usage,my_palette,color_codes_all)
plot_atp6 <- plot_codon_heatmap(atp6_usage,my_palette,color_codes_atp6)
plot_atp8 <- plot_codon_heatmap(atp8_usage,my_palette,color_codes_atp8) ####
plot_cob <- plot_codon_heatmap(cob_usage,my_palette,color_codes_cob)
plot_cox1 <- plot_codon_heatmap(cox1_usage,my_palette,color_codes_cox1)
plot_cox2 <- plot_codon_heatmap(cox2_usage,my_palette,color_codes_cox2)
plot_cox3 <- plot_codon_heatmap(cox3_usage,my_palette,color_codes_cox3)
plot_nad1 <- plot_codon_heatmap(nad1_usage,my_palette,color_codes_nad1)
plot_nad2 <- plot_codon_heatmap(nad2_usage,my_palette,color_codes_nad2)
plot_nad3 <- plot_codon_heatmap(nad3_usage,my_palette,color_codes_nad3) ####
plot_nad4 <- plot_codon_heatmap(nad4_usage,my_palette,color_codes_nad4)
plot_nad4l <- plot_codon_heatmap(nad4l_usage,my_palette,color_codes_nad4l) ####
plot_nad5 <- plot_codon_heatmap(nad5_usage,my_palette,color_codes_nad5)
plot_nad6 <- plot_codon_heatmap(nad6_usage,my_palette,color_codes_nad6)  #####

