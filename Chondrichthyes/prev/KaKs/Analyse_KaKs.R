source("functions.R")

# M0 model of PAML


# Create a lookup table for accession and taxonomy_rank
look_up_table <- aggregate(Taxonomy_Rank ~ Accession, PCG, unique)
look_up_table_chimaeriformes <- aggregate(Family ~ Accession, PCG, unique)
look_up_table_chimaeriformes[look_up_table_chimaeriformes == ''] <- NA
look_up_table_chimaeriformes <- look_up_table_chimaeriformes[complete.cases(look_up_table_chimaeriformes),]

gene_list <- c("atp6","atp8","cob","cox1","cox2","cox3","nad1","nad2","nad3","nad4","nad4l","nad5","nad6")
dnds_mean <- get_mean_dnds(gene_list, "mean")
dnds_sd <- get_mean_dnds(gene_list, "sd")

dnds <- data.frame(gene_list,dnds_mean,dnds_sd)

ggplot(dnds, aes(x=factor(gene_list), y=dnds_mean)) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    ggtitle('Chondrichthyes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
         axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_mean-dnds_sd/6, ymax=dnds_mean+dnds_sd/6), width=.1,
                  position=position_dodge(.9))



dnds_myliobatiformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Myliobatiformes","Myliobatiformes")
dnds_myliobatiformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Myliobatiformes","Myliobatiformes")

dnds_myliobatiformes <- data.frame(gene_list,dnds_myliobatiformes_mean,dnds_myliobatiformes_sd)

ggplot(dnds_myliobatiformes, aes(x=factor(gene_list), y=dnds_myliobatiformes_mean, fill = '#f50000')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#f50000') +
    ggtitle('Myliobatiformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_myliobatiformes_mean-dnds_myliobatiformes_sd/6, ymax=dnds_myliobatiformes_mean+dnds_myliobatiformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_rajiformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Rajiformes","Rajiformes")
dnds_rajiformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Rajiformes","Rajiformes")

dnds_rajiformes <- data.frame(gene_list,dnds_rajiformes_mean,dnds_rajiformes_sd)

ggplot(dnds_rajiformes, aes(x=factor(gene_list), y=dnds_rajiformes_mean, fill = '#cc0000')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#cc0000') +
    ggtitle('Rajiformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_rajiformes_mean-dnds_rajiformes_sd/6, ymax=dnds_rajiformes_mean+dnds_rajiformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_rhinopristiformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Rhinopristiformes","Rhinopristiformes")
dnds_rhinopristiformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Rhinopristiformes","Rhinopristiformes")

dnds_rhinopristiformes <- data.frame(gene_list,dnds_rhinopristiformes_mean,dnds_rhinopristiformes_sd)

ggplot(dnds_rhinopristiformes, aes(x=factor(gene_list), y=dnds_rhinopristiformes_mean, fill = '#a30000')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#a30000') +
    ggtitle('Rhinopristiformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_rhinopristiformes_mean-dnds_rhinopristiformes_sd/6, ymax=dnds_rhinopristiformes_mean+dnds_rhinopristiformes_sd/6), width=.1,
                  position=position_dodge(.9))


dnds_torpediniformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Torpediniformes","Torpediniformes")
dnds_torpediniformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Torpediniformes","Torpediniformes")

dnds_torpediniformes <- data.frame(gene_list,dnds_torpediniformes_mean,dnds_torpediniformes_sd)

ggplot(dnds_torpediniformes, aes(x=factor(gene_list), y=dnds_torpediniformes_mean, fill = '#4e1609')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#4e1609') +
    ggtitle('Torpediniformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_torpediniformes_mean-dnds_torpediniformes_sd/6, ymax=dnds_torpediniformes_mean+dnds_torpediniformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_carcharhiniformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Carcharhiniformes","Carcharhiniformes")
dnds_carcharhiniformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Carcharhiniformes","Carcharhiniformes")

dnds_carcharhiniformes <- data.frame(gene_list,dnds_carcharhiniformes_mean,dnds_rhinopristiformes_sd)

ggplot(dnds_carcharhiniformes, aes(x=factor(gene_list), y=dnds_carcharhiniformes_mean, fill = '#4f86f7')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#4f86f7') +
    ggtitle('Carcharhiniformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_carcharhiniformes_mean-dnds_carcharhiniformes_sd/6, ymax=dnds_carcharhiniformes_mean+dnds_rhinopristiformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_lamniformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Lamniformes","Lamniformes")
dnds_lamniformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Lamniformes","Lamniformes")

dnds_lamniformes <- data.frame(gene_list,dnds_lamniformes_mean,dnds_lamniformes_sd)

ggplot(dnds_lamniformes, aes(x=factor(gene_list), y=dnds_lamniformes_mean, fill = '#45b1e8')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#45b1e8') +
    ggtitle('Lamniformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_lamniformes_mean-dnds_lamniformes_sd/6, ymax=dnds_lamniformes_mean+dnds_lamniformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_orectolobiformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Orectolobiformes","Orectolobiformes")
dnds_orectolobiformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Orectolobiformes","Orectolobiformes")

dnds_orectolobiformes <- data.frame(gene_list,dnds_orectolobiformes_mean,dnds_orectolobiformes_sd)

ggplot(dnds_orectolobiformes, aes(x=factor(gene_list), y=dnds_orectolobiformes_mean, fill = '#87ceeb')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#87ceeb') +
    ggtitle('Orectolobiformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_orectolobiformes_mean-dnds_orectolobiformes_sd/6, ymax=dnds_orectolobiformes_mean+dnds_orectolobiformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_heterodontiformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Heterodontiformes","Heterodontiformes")
dnds_heterodontiformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Heterodontiformes","Heterodontiformes")

dnds_heterodontiformes <- data.frame(gene_list,dnds_heterodontiformes_mean,dnds_heterodontiformes_sd)

ggplot(dnds_heterodontiformes, aes(x=factor(gene_list), y=dnds_heterodontiformes_mean, fill = '#30d5c8')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#30d5c8') +
    ggtitle('Heterodontiformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_heterodontiformes_mean-dnds_heterodontiformes_sd/6, ymax=dnds_heterodontiformes_mean+dnds_heterodontiformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_hexachinformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Hexachinformes","Hexachinformes")
dnds_hexachinformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Hexachinformes","Hexachinformes")

dnds_hexachinformes <- data.frame(gene_list,dnds_hexachinformes_mean,dnds_hexachinformes_sd)

ggplot(dnds_hexachinformes, aes(x=factor(gene_list), y=dnds_hexachinformes_mean, fill = '#adff2f')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#adff2f') +
    ggtitle('Hexachinformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_hexachinformes_mean-dnds_hexachinformes_sd/6, ymax=dnds_hexachinformes_mean+dnds_hexachinformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_pristiophoriformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Pristiophoriformes","Pristiophoriformes")
dnds_pristiophoriformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Pristiophoriformes","Pristiophoriformes")

dnds_pristiophoriformes <- data.frame(gene_list,dnds_pristiophoriformes_mean,dnds_pristiophoriformes_sd)

ggplot(dnds_pristiophoriformes, aes(x=factor(gene_list), y=dnds_pristiophoriformes_mean, fill = '#32cd32')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#32cd32') +
    ggtitle('Pristiophoriformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_pristiophoriformes_mean-dnds_pristiophoriformes_sd/6, ymax=dnds_pristiophoriformes_mean+dnds_pristiophoriformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_squaliformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Squaliformes","Squaliformes")
dnds_squaliformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Squaliformes","Squaliformes")

dnds_squaliformes <- data.frame(gene_list,dnds_squaliformes_mean,dnds_squaliformes_sd)

ggplot(dnds_squaliformes, aes(x=factor(gene_list), y=dnds_squaliformes_mean, fill = '#3cb371')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#3cb371') +
    ggtitle('Squaliformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_squaliformes_mean-dnds_squaliformes_sd/6, ymax=dnds_squaliformes_mean+dnds_squaliformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_squatiniformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Squatiniformes","Squatiniformes")
dnds_squatiniformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Squatiniformes","Squatiniformes")

dnds_squatiniformes <- data.frame(gene_list,dnds_squatiniformes_mean,dnds_squatiniformes_sd)

ggplot(dnds_squatiniformes, aes(x=factor(gene_list), y=dnds_squatiniformes_mean, fill = '#006400')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#006400') +
    ggtitle('Squatiniformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_squatiniformes_mean-dnds_squatiniformes_sd/6, ymax=dnds_squatiniformes_mean+dnds_squatiniformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_chimaeriformes_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Chimaeriformes","Chimaeriformes")
dnds_chimaeriformes_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Chimaeriformes","Chimaeriformes")

dnds_chimaeriformes <- data.frame(gene_list,dnds_chimaeriformes_mean,dnds_chimaeriformes_sd)

ggplot(dnds_chimaeriformes, aes(x=factor(gene_list), y=dnds_chimaeriformes_mean, fill = '#7b68ee')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#7b68ee') +
    ggtitle('Chimaeriformes') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_chimaeriformes_mean-dnds_chimaeriformes_sd/6, ymax=dnds_chimaeriformes_mean+dnds_chimaeriformes_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_callorhinchidae_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table, "Callorhinchidae","Callorhinchidae")
dnds_callorhinchidae_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table, "Callorhinchidae","Callorhinchidae")

dnds_callorhinchidae <- data.frame(gene_list,dnds_callorhinchidae_mean,dnds_callorhinchidae_sd)

ggplot(dnds_callorhinchidae, aes(x=factor(gene_list), y=dnds_callorhinchidae_mean, fill = '#00AFBB')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#00AFBB') +
    ggtitle('Callorhinchidae') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_callorhinchidae_mean-dnds_callorhinchidae_sd/6, ymax=dnds_callorhinchidae_mean+dnds_callorhinchidae_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_chimaeridae_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table_chimaeriformes, "Chimaeridae","Chimaeridae")
dnds_chimaeridae_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table_chimaeriformes, "Chimaeridae","Chimaeridae")

dnds_chimaeridae <- data.frame(gene_list,dnds_chimaeridae_mean,dnds_chimaeridae_sd)

ggplot(dnds_chimaeridae, aes(x=factor(gene_list), y=dnds_chimaeridae_mean, fill = '#E7B800')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#E7B800') +
    ggtitle('Chimaeridae') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_chimaeridae_mean-dnds_chimaeridae_sd/6, ymax=dnds_chimaeridae_mean+dnds_chimaeridae_sd/6), width=.1,
                  position=position_dodge(.9))

dnds_rhinochimaeridae_mean <- get_dnds_pergroup(gene_list, "mean", look_up_table_chimaeriformes, "Rhinochimaeridae","Rhinochimaeridae")
dnds_rhinochimaeridae_sd <- get_dnds_pergroup(gene_list, "sd", look_up_table_chimaeriformes, "Rhinochimaeridae","Rhinochimaeridae")

dnds_rhinochimaeridae <- data.frame(gene_list,dnds_rhinochimaeridae_mean,dnds_rhinochimaeridae_sd)

ggplot(dnds_rhinochimaeridae, aes(x=factor(gene_list), y=dnds_rhinochimaeridae_mean, fill = '#FC4E07')) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    scale_fill_manual(values = '#FC4E07') +
    ggtitle('Rhinochimaeridae') +
    labs(x = 'Protein-Coding Genes', y = 'dN/dS') +
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=dnds_rhinochimaeridae_mean-dnds_rhinochimaeridae_sd/6, ymax=dnds_rhinochimaeridae_mean+dnds_rhinochimaeridae_sd/6), width=.1,
                  position=position_dodge(.9))
