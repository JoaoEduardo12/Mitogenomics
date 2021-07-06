
get_taxonomy_column <- function(tax_vector, dataset, column) {
    # "tax_vector is a vector that contains the taxonomic ranks to isolate
    # "dataset" is the dataset used to find taxonomic ranks
    # Given a vector of taxonomic ranks, it creates a new column in "dataset" with just
    # one specified taxonomic rank found in "tax_vector"
    dataset[, column] <- ''
    for(row in 1:nrow(dataset)) {
        for (element in tax_vector) {  
            if (element %in% str_split(dataset[row,which(colnames(dataset) == "Taxonomy")], ';')[[1]]) {
                dataset[row,length(dataset)] <- element
            }
        }
    }
    dataset
}

create_meansd_dataframe <- function(my_dataset,vector,column,my_column) {
    # Creates a new dataframe with only the mean and sd values of a specified column and
    # choosing a particular taxonomy to display the results
    # my_dataset = the dataset in use for the analysis
    # vector = taxonomy names
    # column = column name of the taxonomy levels to separate
    # my_column = target column name to which extract the mean and sd
    taxonomy_names <- c()
    mean_values <- c()
    sd_values <- c()
    for (element in vector) { 
        new_dataset <- subset(my_dataset, my_dataset[,which(colnames(my_dataset) == column)] == element)
        taxonomy_names <- c(taxonomy_names, element)
        mean_values <- c(mean_values, mean(new_dataset[[my_column]]))
        sd_values <- c(sd_values, sd(new_dataset[[my_column]]))
        df <- data.frame(taxonomy_names,mean_values,sd_values)
    }
    df
}

create_gene_dataframe <- function(my_dataset,vector, column, gene_vector, my_column) {
    # Creates a new dataframe with only the mean and sd values of a specified column and
    # choosing a particular taxonomy to display the results, but integrates genes instead of
    # complete mitogenomes
    # my_dataset = the dataset in use for the analysis
    # vector = taxonomy names
    # column = column name of the taxonomy levels to separate
    # gene_vector = vector with the gene names
    # my_column = target column name to which extract the mean and sd
    taxonomy_names <- c()
    gene_names <- c()
    mean_values <- c()
    sd_values <- c()
    for (element in vector) {
        new_dataset <- subset(my_dataset, my_dataset[, which(colnames(my_dataset) == column)] == element)
        for (gene in gene_vector) {
            gene_dataset <- subset(new_dataset, new_dataset[, which(colnames(new_dataset) == "Content")] == gene)
            taxonomy_names <- c(taxonomy_names, element)
            gene_names <- c(gene_names, gene)
            mean_values <- c(mean_values, mean(gene_dataset[[my_column]]))
            sd_values <- c(sd_values, sd(gene_dataset[[my_column]]))
        }
    }
    df <- data.frame(taxonomy_names,gene_names,mean_values,sd_values)
    df
}

make_comp_def <- function(main_df, col1, col2, col3 = 'Length') {
    for (id in unique(main_df[,col1])) {
        for (lines in subset(main_df, main_df[,col1] == id)) {
            
        }
    }
}

get_nucleotide_weighted_mean <- function(dataset) {
    # Creates additonal columns with the mean nucleotide percentage weighted over its
    # gene length, so the results can be comprared fairly
    # uses tidyverse
    PCG_weighted_mean <- dataset %>%
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
    PCG_weighted_mean
}

get_nucleotide_frequencies <- function(column, dataset, my_function, tax) {
    # gets unique values from nucleotide frequencies per a specified column
    # column = is the column to which get the unique values
    # my_function = is the function to be applied to all members
    base_vector <- c("G1","C1","A1","T1","G2","C2","A2","T2","G3","C3","A3","T3",tax)
    for (i in base_vector) {
        if (i == "G1") {
            nucleotide_vector <- aggregate(dataset[[i]] ~ dataset[[column]], dataset, my_function)
            colnames(nucleotide_vector) <- c("Accession", i)
        } else {
            intermidiate <- aggregate(dataset[[i]] ~ dataset[[column]], dataset, my_function)
            colnames(intermidiate) <- c("Accession", i)
            nucleotide_vector[[i]] <- cbind(intermidiate[[i]])
        }
    }
    nucleotide_vector
}

ggplotRegression <- function (fit) {
    
    require(ggplot2)
    
    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
        geom_point(col = "blue") +
        stat_smooth(method = "lm", col = "red") +
        labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                           "Intercept =",signif(fit$coef[[1]],5 ),
                           " Slope =",signif(fit$coef[[2]], 5),
                           " P =",signif(summary(fit)$coef[2,4], 5)))
}

get_mean_dnds <- function(gene_list, method) {
    final_vector <- c()
    for (gene in gene_list) {
        temp <- read.csv(paste0("/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/prev/KaKs/",gene,"/",gene, "_dnds.csv"))
        temp_2 <- temp[,-1]
        rownames(temp_2) <- temp[,1]
        temp_2 <- as.matrix(temp_2)
        temp_2[temp_2 > 3] <- NA
        if (method == "mean") {
            final_vector <- c(final_vector, mean(temp_2, na.rm = TRUE))
        } else {
            final_vector <- c(final_vector, sd(temp_2, na.rm = TRUE))
        }
    }
    final_vector
}

get_dnds_pergroup <- function(gene_list, method, look_up_table, group1, group2) {
    final_vector <- c()
    invalid_accessions <- c()
    if (group1 == group2) {
        for (i in 1:nrow(look_up_table)) {
            if (look_up_table[i,2] != group1) {
                invalid_accessions <- c(invalid_accessions, look_up_table[i,1])
            }
        }
        for (gene in gene_list) {
            temp <- read.csv(paste0("/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/prev/KaKs/",gene,"/",gene, "_dnds.csv"))
            temp_2 <- temp[,-1]
            rownames(temp_2) <- temp[,1]
            temp_2[temp_2 > 3] <- NA
            cont <- 1
            for (column in colnames(temp_2)) {
                if (column %in% invalid_accessions) {
                    temp_2[[column]] <- NA
                    temp_2[cont] <- NA
                }
                cont <- cont + 1
            }
            if (method == "mean") {
                final_vector <- c(final_vector, mean(as.matrix(temp_2), na.rm = TRUE))
            } else {
                final_vector <- c(final_vector, sd(as.matrix(temp_2), na.rm = TRUE))
            }
        }
    }
    final_vector
}


create_codon_dataframe <- function(my_dataset,taxonomy_name,column,columns,new_name, genetic_code) {
    # Creates a new dataframe with only the mean and sd values of a specified column and
    # choosing a particular taxonomy to display the results
    # my_dataset = the dataset in use for the analysis
    # vector = taxonomy names
    # column = column name of the taxonomy levels to separate
    # my_column = target column name to which extract the mean and sd
    require(Biostrings)
    new_dataset <- subset(my_dataset, my_dataset[,which(colnames(my_dataset) == column)] == taxonomy_name)
    mean_values <- colMeans(new_dataset[,columns])
    codons <- names(mean_values)
    mean_values <- unname(mean_values)
    df <- data.frame(codons, mean_values)
    aminoacids <- c()
    genetic_code_values <- getGeneticCode(genetic_code)
    for (row in df[,1]) {
        aminoacids <- c(aminoacids, toString(translate(DNAString(row), genetic.code = genetic_code_values)))
    }
    df[[new_name]] <- aminoacids
    df
}

plot_rscu <- function(codon_usage, name) {
    library(dplyr)
    library(patchwork)
    library(forcats)
    
    #codon_usage[["Col"]] <- as.integer(factor(codon_usage[["aminoacids"]], levels = unique(codon_usage[["aminoacids"]])))
    #codon_usage <- codon_usage %>%
        #mutate(Col = as.factor(Col))
    #codon_usage <- codon_usage %>%
    #    group_by(aminoacids) %>%
     #   mutate(Col = lag(cumsum(Col), default = 0)) %>%
      #  mutate(Col = as.factor(Col))
    rep <- c()
    colors <- c()
    for (element in codon_usage[["aminoacids"]]) {
        rep <- c(rep, element)
        colors <- c(colors, length(which(rep == element)))
    }
    codon_usage["Col"] <- colors
    codon_usage <- codon_usage %>%
        mutate(Col = as.factor(Col))
     
    p0 <- ggplot(codon_usage, aes(x=factor(aminoacids), y = mean_values, fill = Col)) +
        geom_bar(stat = "identity", color = '#505050', width = .7, show.legend = FALSE) +
        ggtitle(name) +
        labs(y = 'RSCU') +
        theme(legend.position = "none", legend.title = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank())
    
    p1 <- ggplot(data = codon_usage) +
        geom_tile(aes(x=aminoacids, y = fct_rev(Col), fill = Col), col = '#505050') +
        geom_text(aes(x = aminoacids, y = Col, label = codons), col = 'white', fontface = 'bold') +
        labs(x = "Aminoacids") +
        theme_void() +
        theme(legend.position = 'none', text = element_text(size = 13),
              axis.text.x = element_text(angle = 0, hjust = 0.5),
              axis.title.y = element_blank()) +
        scale_x_discrete(labels = c('Ala','Cys','Asp','Glu','Phe','Gly','His','Lys','Leu','Met','Asn','Pro','Gln','Arg','Ser','Thr','Val','Trp','Tyr'))
    
    p= p0/p1
    p = p + plot_layout(heights = c(4,1))
    p
}


pairwise_statistical_test <- function(group, column, test) {
    pvalues <- c()
    names <- c()
    for (i in seq(group)) {
        for (j in seq(group)) {
            if (j > i) {
                pvalues <- c(pvalues, eval(parse(text = test))(eval(parse(text = paste(group[[i]], "$", column, sep = ''))), eval(parse(text = paste(group[[j]], "$", column, sep = ''))))$p.value)
                names <- c(names, paste(strsplit(group[i],"_codon_usage"), " x ", paste(strsplit(group[j],"_codon_usage"))))
            } 
        }
    }
    df <- data.frame(names, pvalues)
    df
}


plot_codon_heatmap <- function(dataset, palette, col_codes) {
    heatmap_rscu <- heatmap.2(x=t(as.matrix(dataset)),
                              trace = 'none',
                              #scale = 'row',
                              ColSideColors = col_codes,
                              distfun = function(x) dist(x, method="euclidean"),
                              hclustfun = function(x) hclust(x, method="average"),
                              key.title = '',
                              key.ylab = '',
                              key.xlab = '',
                              colCol = NULL,
                              xlab = '',
                              denscol = 'grey',
                              density.info = 'density',
                              col = palette,
                              keysize = 1.5,
                              labCol = '',
                              key.par=list(),
                              cexRow = 0.7)
    heatmap_rscu
}


make_gene_codondataset <- function(dataset, gene) {
    new_dataset <- subset(dataset, dataset[["Gene"]] == gene)
    new_dataset <- data.frame(new_dataset)
    rownames(new_dataset) <- new_dataset[,1]
    new_dataset[["Upper_Rank"]] <- NULL
    new_dataset[["Family"]] <- NULL
    new_dataset[["Gene"]] <- NULL
    new_dataset
}
