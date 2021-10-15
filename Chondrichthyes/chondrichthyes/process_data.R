### Script to analyse mitogenomes ###

####
# This are 4 main scripts to analyse our data. The first depicts all pure visualisation modules
# while the rest is to perform statistical analyses, analyse codeml output and apply machine learning modules

# import function scripts
source("/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/chondrichthyes/functions.R")

##### Data Processing ####
# Set the right working directory
setwd("~/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Deep_Mitogenomics")

# Import necessary library to open csv file
library(readr)
# Load tidyverse packages
library(plyr)
library(dplyr)
library(stringr)
# load ggplot
library(ggplot2)
library(EnvStats)
library(grid)
library("gplots")
library(devtools)
# creating linear models
library(caret)
library("FactoMineR")
library("factoextra")
# library(ggbiplot)

# Read the csv file
# This file has compositional analyses
dataset <- read_csv("Chondrichthyes.csv")    # EDIT HERE TO UPLOAD YOUR FILE
# This file has codon properties
codon_dataset <- read.csv("Chondrichthyes_codons_.csv") # EDIT HERE TO UPLOAD YOUR FILE

# We need to seperate this dataset into complete mitogenomes dataset and protein-coding dataset, and gather the useful information from both datasets
# Creating a vector to isolate protein coding genes 
PCG_name_vector <- c('atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6')

# Isolating all protein coding genes
PCG <- subset(dataset, dataset$Content %in% PCG_name_vector)

# merge the 2 datasets; PCG and codon_dataset
PCG[,34:93] <- codon_dataset[,7:66]

# Creating new datasets according to each individual protein coding gene, for later purposes :)
atp6 = subset(PCG, PCG$Content == "atp6")
atp8 = subset(PCG, PCG$Content == "atp8")
cox1 = subset(PCG, PCG$Content == "cox1")
cox2 = subset(PCG, PCG$Content == "cox2")
cox3 = subset(PCG, PCG$Content == "cox3")
cob = subset(PCG, PCG$Content == "cob")
nad1 = subset(PCG, PCG$Content == "nad1")
nad2 = subset(PCG, PCG$Content == "nad2")
nad3 = subset(PCG, PCG$Content == "nad3")
nad4 = subset(PCG, PCG$Content == "nad4")
nad5 = subset(PCG, PCG$Content == "nad5")
nad6 = subset(PCG, PCG$Content == "nad6")
nad4l = subset(PCG, PCG$Content == "nad4l")

# Vector to use in function to isolate specific taxonomic ranks
# Here are some examples of what you can isolate based on the column "Taxonomy_Rank":
# taxonomy_vector <- c('Elasmobranchii', 'Holocephali') # Subclasses
# taxonomy_vector <- c('Galeomorphii','Squalomorphii','Batoidea','Chimaeridae','Rhinochimaeridae','Callorhinchidae') # Superorders
# taxonomy_vector <- c('Callorhinchidae','Chimaeridae','Rhinochimaeridae') # Families
# taxonomy_vector <- c('Batoidea','Galeomorphii','Squalomorphii','Holocephali') # optimal subgroups

# We are using this vector, which isolates every Order in Chondrichthyes
taxonomy_vector <- c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes','Batoidea incertae sedis') # Orders
# this is to figure out if sequencing tech was producing different results from the measured variables
tech_vector <- c("Sanger dideoxy sequencing", "454", "Illumina","IonTorrent","Nanopore","PacBio") # the sequencing technology used

# Get extra columns in the dataset with values of specified taxonomic ranks
# In this case I want info on Order
dataset <- get_taxonomy_column(taxonomy_vector, dataset, 'Order')
# Same for genes
PCG <- get_taxonomy_column(taxonomy_vector, PCG, 'Order')


# Isolating mitogenomes in a single dataset
mitogenomes <- subset(dataset, dataset$Content == 'Mitogenome') # there are 3 mitogenomes which do not information regarding their genes in this particular dataset

# Create a particular dataset with length mean and sd that we can use to plot in ggplot
mitogenome_len <- create_meansd_dataframe(mitogenomes, taxonomy_vector, "Order", "Length")

# Integrated NGS data so we may see differences between different NGS techs
mitogenome_len_seq <- create_meansd_dataframe(mitogenomes, tech_vector, "Technology", "Length")

# Create unique gene data frame with mean and sd values that can later be used in ggplot
gene_len <- create_gene_dataframe(PCG, taxonomy_vector, 'Order', PCG_name_vector, 'Length')

# Creating a dataframe to seperate the genes according to taxonomic rank
#The following line separates by superclass
gene_gc <- create_gene_dataframe(PCG, taxonomy_vector, 'Order', PCG_name_vector, 'GC_content')

# ??? I dunno what I did there
group_codon_usage <- PCG %>% 
    group_by(Order, colnames(PCG[7:66]))

# this is for creating specific datasets with RSCU values per Chondrichthyes Order
myliobatiformes_codon_usage <- create_codon_dataframe(PCG, 'Myliobatiformes', 'Order',35:93,'aminoacids',"2")
rajiformes_codon_usage <- create_codon_dataframe(PCG, 'Rajiformes', 'Order',35:93,'aminoacids',"2")
rhinopristiformes_codon_usage <- create_codon_dataframe(PCG, 'Rhinopristiformes', 'Order',35:93,'aminoacids',"2")
torpediniformes_codon_usage <- create_codon_dataframe(PCG, 'Torpediniformes', 'Order',35:93,'aminoacids',"2")
carcharhiniformes_codon_usage <- create_codon_dataframe(PCG, 'Carcharhiniformes', 'Order',35:93,'aminoacids',"2")
lamniformes_codon_usage <- create_codon_dataframe(PCG, 'Lamniformes', 'Order',35:93,'aminoacids',"2")
orectolobiformes_codon_usage <- create_codon_dataframe(PCG, 'Orectolobiformes', 'Order',35:93,'aminoacids',"2")
heterodontiformes_codon_usage <- create_codon_dataframe(PCG, 'Heterodontiformes', 'Order',35:93,'aminoacids',"2")
hexanchiformes_codon_usage <- create_codon_dataframe(PCG, 'Hexanchiformes', 'Order',35:93,'aminoacids',"2")
pristiophoriformes_codon_usage <- create_codon_dataframe(PCG, 'Pristiophoriformes', 'Order',35:93,'aminoacids',"2")
squaliformes_codon_usage <- create_codon_dataframe(PCG, 'Squaliformes', 'Order',35:93,'aminoacids',"2")
squatiniformes_codon_usage <- create_codon_dataframe(PCG, 'Squatiniformes', 'Order',35:93,'aminoacids',"2")
chimaeriformes_codon_usage <- create_codon_dataframe(PCG, 'Chimaeriformes', 'Order',35:93,'aminoacids',"2")
atp6_codon_usage <- create_codon_dataframe(PCG, 'atp6', 'Content',35:93,'aminoacids',"2")
atp8_codon_usage <- create_codon_dataframe(PCG, 'atp8', 'Content',35:93,'aminoacids',"2")
cox1_codon_usage <- create_codon_dataframe(PCG, 'cox1', 'Content',35:93,'aminoacids',"2")
cox2_codon_usage <- create_codon_dataframe(PCG, 'cox2', 'Content',35:93,'aminoacids',"2")
cox3_codon_usage <- create_codon_dataframe(PCG, 'cox3', 'Content',35:93,'aminoacids',"2")
nad1_codon_usage <- create_codon_dataframe(PCG, 'nad1', 'Content',35:93,'aminoacids',"2")
nad2_codon_usage <- create_codon_dataframe(PCG, 'nad2', 'Content',35:93,'aminoacids',"2")
nad3_codon_usage <- create_codon_dataframe(PCG, 'nad3', 'Content',35:93,'aminoacids',"2")
nad4_codon_usage <- create_codon_dataframe(PCG, 'nad4', 'Content',35:93,'aminoacids',"2")
nad5_codon_usage <- create_codon_dataframe(PCG, 'nad5', 'Content',35:93,'aminoacids',"2")
nad6_codon_usage <- create_codon_dataframe(PCG, 'nad6', 'Content',35:93,'aminoacids',"2")
nad4l_codon_usage <- create_codon_dataframe(PCG, 'nad4l', 'Content',35:93,'aminoacids',"2")
cob_codon_usage <- create_codon_dataframe(PCG, 'cob', 'Content',35:93,'aminoacids',"2")


## Just vectors to later perform statistical tests

class_cu <- c("myliobatiformes_codon_usage", "rajiformes_codon_usage", "rhinopristiformes_codon_usage", "torpediniformes_codon_usage", "carcharhiniformes_codon_usage",
              "lamniformes_codon_usage", "orectolobiformes_codon_usage", "heterodontiformes_codon_usage", "hexanchiformes_codon_usage", "pristiophoriformes_codon_usage",
              "squaliformes_codon_usage", "squatiniformes_codon_usage", "chimaeriformes_codon_usage")

gene_cu <- c("atp6_codon_usage", "atp8_codon_usage", "cob_codon_usage", "cox1_codon_usage", "cox2_codon_usage",
             "cox3_codon_usage", "nad1_codon_usage", "nad2_codon_usage", "nad3_codon_usage", "nad4_codon_usage",
             "nad4l_codon_usage", "nad5_codon_usage", "nad6_codon_usage")

# Create normalised columns of nucleotide percentages of each accession number by length
weighted_PCG <- get_nucleotide_weighted_mean(PCG) 

# Creates unique values of weighted PCG nucleotide content
nucleotide_vector <- get_nucleotide_frequencies("Accession", weighted_PCG, unique, 'Order') 

# PCA with prcomp
nucleotide_vector_complete <- subset(nucleotide_vector, nucleotide_vector$Order != 'Batoidea incertae sedis')

group_codon_usage <- PCG %>%
    group_by(Accession, Order) %>%
    summarise_at(c(35:93), mean, na.rm = TRUE)

group_gene_usage <- PCG %>%
    group_by(Content, Order) %>%
    summarise_at(c(35:93), mean, na.rm = TRUE)



# Preparing the dataset for heatmaps for RSCU values
group_codon_usage <- data.frame(group_codon_usage)
rownames(group_codon_usage) <- group_codon_usage[,1]

all_codon_usage <- group_codon_usage
all_codon_usage$Accession <- NULL
all_codon_usage$Order <- NULL
all_codon_usage$Upper_Rank <- NULL

acession_per_taxonomy <- PCG %>%
    group_by(Accession, Order) %>%
    summarise(Order)

acession_per_taxonomy <- aggregate(Order ~ Accession, PCG, unique)



some_col_func <- function(n) rev(colorspace::heat_hcl(n, h = c(90,180), c = c(100, 30), l = c(30, 80), power = c(1/5, 1.5)))

par(cex.main=0.05, mar = rep(30,4))

my_palette <- colorRampPalette(c("#fdfd96","#FF6700"))(n = 1000)

library(hash)
h <- hash()
h[c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes',
    'Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes',
    'Chimaeriformes','Batoidea incertae sedis')] <- c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                                                      '#006400','#7b68ee',"black")

color_codes_all <- c()
for (element in acession_per_taxonomy[,2]) {
    color_codes_all <- c(color_codes_all, values(h, key = element))
}

codon_df$col <- ""
for(i in 1:nrow(codon_df)) {
    if (codon_df[i,2] == "AT") {
        codon_df[i,3] = "#77dd77"
    } else {
        codon_df[i,3] = "#ff6961"
    }
}

codon_color_codes <- as.character(codon_df$col)
names(codon_color_codes) <- codon_df$codons

## all for now :)
