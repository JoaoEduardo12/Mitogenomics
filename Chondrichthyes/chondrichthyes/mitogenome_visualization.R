
#### Load dataset #####
source("get_dataset.R")
# This scripts imports the dataset but also performs some processing so hold tight!

##### GC skew vs AT skew #####

# Scatter plot com GC content a cores
ggplot(data=mitogenomes, aes(x = GC_skew, y = AT_skew)) +
    xlab("GC skew") +  ylab("AT skew") +
    geom_point(aes(color = GC_content, shape = factor(Order, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Echinorhiniformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes'))), size = 2.5, alpha = I(2/2)) +
    geom_vline(aes(xintercept = mean(GC_skew)), color = "red", linetype = "dashed") +
    geom_hline(aes(yintercept = mean(AT_skew)), color = "red", linetype = "dashed") +
    scale_color_gradient(low = "yellow", high = "red") +
    labs(color = 'GC content', shape = 'Taxonomy') +
    xlim(c(-0.5,-0.1)) +
    ggtitle('Compositional complete mitogenome differences') #+
    #scale_shape_manual(values=c(19,17,15,23,3,10))

# Simple scatter plot for the total mitogenomes
# You can't run this portion of the code if you didnt name Upper_Rank as one of the columns by specifying a specific taxonomic_group corresponding to subclasses
ggplot(data=subset(mitogenomes, mitogenomes$Order != 'Batoidea incertae sedis'), aes(x = GC_skew, y = AT_skew)) +
   xlab("GC skew") +  ylab("AT skew") +
   geom_point(aes(color = factor(Order, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes',
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

# ditto
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

# Plot mitogenome length for all taxonomy ranks
ggplot(mitogenome_len[1:13,], aes(x=factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), y=mean_values, fill = factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')))) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    ggtitle('Mitogenome Length') +
    scale_fill_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                               '#006400','#7b68ee')) +
    labs(x = 'Taxonomic ranks', y = 'Mitogenome Length') +
    annotate("text", x = 1, y = 19500, label = "n=41") + annotate("text", x = 2, y = 19500, label = "n=64") + annotate("text", x = 3, y = 19500, label = "n=12") + annotate("text", x = 4, y = 19500, label = "n=4") + annotate("text", x = 5, y = 19500, label = "n=260") + annotate("text", x = 6, y = 19500, label = "n=172") + annotate("text", x = 7, y = 19500, label = "n=13") + annotate("text", x =8, y = 19500, label = "n=2") + annotate("text", x = 9, y = 19500, label = "n=7") + annotate("text", x = 10, y = 19500, label = "n=1") + annotate("text", x = 11, y = 19500, label = "n=12") + annotate("text", x = 12, y = 19500, label = "n=6") + annotate("text", x = 13, y= 23000, label = "n=11") + 
    #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.1,
                  position=position_dodge(.9))

# Mitogenome length according to each sequencing technology
ggplot(mitogenome_len_seq, aes(x=factor(taxonomy_names, levels = c('Sanger dideoxy sequencing','454','Illumina','IonTorrent','Nanopore','PacBio')), y=mean_values, fill = factor(taxonomy_names, levels = c('Sanger dideoxy sequencing','454','Illumina','IonTorrent','Nanopore','PacBio')))) +
    geom_bar(stat="identity",position=position_dodge(), color = '#505050', width = .7) +
    ggtitle('Mitogenome Length') +
    labs(x = 'Sequencing Technology', y = 'Mitogenome Length') +
    annotate("text", x = 1, y = 18500, label = "n=29") + annotate("text", x = 2, y = 18500, label = "n=3") + annotate("text", x = 3, y = 18500, label = "n=257") + annotate("text", x = 4, y = 18500, label = "n=16") + annotate("text", x = 5, y = 18500, label = "n=6") + annotate("text", x = 6, y = 18500, label = "n=1") +  #theme_classic() +
    theme(text = element_text(size=13),
          axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
    geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.1,
                  position=position_dodge(.9))


# Gene length by Chondrichthyes Order
ggplot(gene_len[1:169,], aes(x = factor(taxonomy_names, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes',
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
    annotate("text", x = 1, y = 2000, label = "n=41") + annotate("text", x = 2, y = 2000, label = "n=64") + annotate("text", x = 3, y = 2000, label = "n=12") + annotate("text", x = 4, y = 2000, label = "n=4") + annotate("text", x = 5, y = 2000, label = "n=260") + annotate("text", x = 6, y = 2000, label = "n=172") + annotate("text", x = 7, y = 2000, label = "n=13") + annotate("text", x =8, y = 2000, label = "n=2") + annotate("text", x = 9, y = 2000, label = "n=7") + annotate("text", x = 10, y = 2000, label = "n=1") + annotate("text", x = 11, y = 2000, label = "n=12") + annotate("text", x = 12, y = 2000, label = "n=6") + annotate("text", x = 13, y= 2000, label = "n=11") + 
    theme(plot.margin = unit(c(1,3,1,1), "lines")) +
    #theme_classic() +
    scale_color_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                                '#006400','#7b68ee')) +
    labs(x = '', y = 'Size', color= 'Gene') +
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle=50, hjust=1)) +
    geom_errorbar(aes(ymin=mean_values-sd_values, ymax=mean_values+sd_values), width=.7,
                  position=position_dodge(0))


##### Mitogenome and Gene GC content ######

# Plotting full mitogenome gc content
ggplot(mitogenomes, aes(x = factor(Order, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), y = GC_content, color = factor(Order, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')))) +
    geom_jitter(aes(x = factor(Order, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes')), color = factor(Order, levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes'))), 
                position = position_jitter(width = .05), alpha = 0.5) +
    geom_boxplot(outlier.colour = NA, position = "dodge") +
    ggtitle('Mitogenome GC content') +
    scale_color_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                                '#006400','#7b68ee')) +
    stat_n_text() +
    labs(x = '', y = 'GC content(%)') +
    theme(text = element_text(size=12),
          axis.text.x = element_text(angle=0, hjust=0.5),
          legend.position="none") #+
    #stat_summary(fun=mean, geom="point", shape=23, size=4)


# Integrated by Sequencing tech
ggplot(mitogenomes, aes(x = factor(Technology, levels = c('Sanger dideoxy sequencing','454','Illumina','IonTorrent','Nanopore','PacBio')), y = GC_content, color = factor(Technology, levels = c('Sanger dideoxy sequencing','454','Illumina','IonTorrent','Nanopore','PacBio')))) +
    geom_jitter(aes(x = factor(Technology, levels = c('Sanger dideoxy sequencing','454','Illumina','IonTorrent','Nanopore','PacBio'))), 
                position = position_jitter(width = .05), alpha = 0.5) +
    geom_boxplot(outlier.colour = NA, position = "dodge") +
    ggtitle('Mitogenome GC content') +
    stat_n_text() +
    labs(x = '', y = 'GC content(%)') +
    theme(text = element_text(size=12),
          axis.text.x = element_text(angle=0, hjust=0.5),
          legend.position="none") #+
#stat_summary(fun=mean, geom="poin

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
    geom_line(size = 0.7) +
    geom_point(size = 1.3) +
    #geom_text(data = subset(gene_len[1:169,], taxonomy_names == "Chimaeriformes"), aes(label = gene_names, x = Inf, y = mean_values), hjust = -.1) +
    #scale_fill_manual(values=c("#fcfc4b", "#3ace3a", "#80CEE1","#ff392e")) +
    ggtitle('Gene GC content') +
    annotate("text", x = 1, y = 51, label = "n=41") + annotate("text", x = 2, y = 51, label = "n=64") + annotate("text", x = 3, y = 51, label = "n=12") + annotate("text", x = 4, y = 51, label = "n=4") + annotate("text", x = 5, y = 51, label = "n=260") + annotate("text", x = 6, y = 51, label = "n=172") + annotate("text", x = 7, y = 51, label = "n=13") + annotate("text", x =8, y = 51, label = "n=2") + annotate("text", x = 9, y = 51, label = "n=7") + annotate("text", x = 10, y = 51, label = "n=1") + annotate("text", x = 11, y = 51, label = "n=12") + annotate("text", x = 12, y = 51, label = "n=6") + annotate("text", x = 13, y= 51, label = "n=11") + 
    theme(plot.margin = unit(c(1,3,1,1), "lines")) +
    #theme_classic() +
    scale_color_manual(values=c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32','#3cb371',
                                '#006400','#7b68ee')) +
    labs(x = '', y = 'GC content', color= 'Gene') +
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle=50, hjust=1))


####### Codon Analyses ##########


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

