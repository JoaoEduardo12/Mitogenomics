#install.packages('reshape2')
library(reshape2)
library(plyr)
#install.packages("devtools")
library(devtools)
#install_github("easyGgplot2", "kassambara")
library(ggplot2)

setwd("~/Desktop/Mitogenomas/All_Mitogenomas/Mitogenome_Statistics")
base_statistics <- read.csv('mitogenome_statistics.xlsx', row.names = 1)

library('ggplot2')

##### Gene Orders ######

base_statistics[base_statistics == ""] = NA

base_statistics <- base_statistics[complete.cases(base_statistics), ]

gene_orders <- tibble::rownames_to_column(base_statistics, "sample")
gene_orders <- gene_orders[,1:2]


gene_orders <- data.frame(table(gene_orders$Gene.Order))
gene_orders$values <- gene_orders$Freq / sum(gene_orders$Freq) * 100

library("dplyr")
gene_orders <- gene_orders %>%
  arrange(desc(Var1)) %>%
  mutate(lab.ypos = cumsum(values) - 0.5*values)

gene_orders[3,4] <- 1.4
gene_orders[1,4] <- 14
gene_orders$Val2 <- c('Unionidae: Gonideinae','Unionidae: Ambleminae; Anodontinae; Discomyiinae; Parreysiinae; Unioninae \nMycetopodidae: Leillinae; Monocondylaeinae \nIridinidae: Aspathariinae; Iridininae','Margaritiferidae: Gibbosulinae; Margaritiferinae')


#mycols <- c("#CD534CFF","#EFC000FF","#0073C2FF")
mycols <- c("#CD534CFF","#EFC000FF","#4191c4")
ggplot(gene_orders, aes(x=2, y=values, fill=factor(Val2, levels = c('Unionidae: Ambleminae; Anodontinae; Discomyiinae; Parreysiinae; Unioninae \nMycetopodidae: Leillinae; Monocondylaeinae \nIridinidae: Aspathariinae; Iridininae','Unionidae: Gonideinae','Margaritiferidae: Gibbosulinae; Margaritiferinae')))) +
  geom_bar(stat="identity",width = 1, color = '#505050') +
  coord_polar(theta = "y", start=0) +
  geom_text(aes(y = lab.ypos, label = factor(Var1, levels = c('UF1', 'UF2','MF1'))), color = "white") +
  theme_void() +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = mycols) +
  labs(x = '', y = 'Frequency of Gene Orders',fill = 'Families: Subfamilies') #+
#geom_text(data=base_statistics$Gene.Order, aes(x=factor(1), y=position, label=quantity) , size=5) + labs(x="", y="")
# theme_void() +
# scale_fill_brewer(palette="Set1")


##### GC content ao longo da sequencia   #########

base_statistics <- read.csv('mitogenome_statistics.xlsx', row.names = 1)



# Boxplot com as frequencias nucleotidicas

nucleotide_freq = base_statistics[,3:6]
colnames(nucleotide_freq) <- c('C','G','T','A')

ggplot(stack(nucleotide_freq), aes(x = ind, y = values, fill = ind)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = 'Nucleotides', y = 'Percentage (%)', fill = 'Nucleotides')

# Scattter plot

base_statistics$AT.content <- 100 - base_statistics$GC.content

#AT skew vs AT content

ggplot(data=base_statistics, aes(x = AT_skew, y = AT.content)) +
  geom_point(aes(color=Mitogenome.Length, shape=factor(Family, levels = c('Mulleriidae','Margaritiferidae','Unionidae','Unionidae:Gonideinae','Iridinidae')))) +
  xlab("AT skew") +  ylab("AT content") +
  geom_point(aes(color = Mitogenome.Length, shape = Family), size = 2.5, alpha = I(2/2)) +
  geom_vline(aes(xintercept = mean(AT_skew)), color = "red", linetype = "dashed") +
  geom_hline(aes(yintercept = mean(AT.content)), color = "red", linetype = "dashed") +
  scale_color_gradient(low = "yellow", high = "red") +
  labs(shape = 'Family',color = 'Mitogenome Length') +
  scale_shape_manual(values=c(19,17,15,23,3)) +
  xlab("AT skew") +  ylab("AT content")

#GC skew vs GC content

ggplot(data=base_statistics, aes(x = GC_skew, y = GC.content)) +
  geom_point(aes(color=Mitogenome.Length, shape=factor(Family, levels = c('Mulleriidae','Margaritiferidae','Unionidae','Unionidae:Gonideinae','Iridinidae')))) +
  xlab("GC skew") +  ylab("GC content") +
  geom_point(aes(color = Mitogenome.Length, shape = Family), size = 2.5, alpha = I(2/2)) +
  geom_vline(aes(xintercept = mean(GC_skew)), color = "red", linetype = "dashed") +
  geom_hline(aes(yintercept = mean(GC.content)), color = "red", linetype = "dashed") +
  scale_color_gradient(low = "yellow", high = "red") +
  labs(shape = 'Family',color = 'Mitogenome Length') +
  scale_shape_manual(values=c(19,17,15,23,3)) +
  xlab("GC skew") +  ylab("GC content")

# Sem Length

ggplot(data=base_statistics, aes(x = GC_skew, y = GC.content)) +
  geom_point(aes(color=factor(Gene.Order, levels = c('UF1','UF2','MF1')), shape=factor(Gene.Order, levels = c('UF1','UF2','MF1'))),size = 2.5, alpha = I(2/2)) +
  xlab("GC skew") +  ylab("GC content") +
  #geom_point(aes(color = Family, shape = Family), size = 2.5, alpha = I(2/2)) +
  geom_vline(aes(xintercept = mean(GC_skew)), color = "gray", linetype = "dashed") +
  geom_hline(aes(yintercept = mean(GC.content)), color = "gray", linetype = "dashed") +
  labs(color = 'Gene Order', shape = 'Gene Order') +
  scale_shape_manual(values = c(15,19,17)) +
  scale_color_manual(values = c('#00BFC4','#7CAE00','#F8766D'))
xlab("GC skew") +  ylab("GC content")

ggplot(data=base_statistics, aes(x = AT_skew, y = AT.content)) +
  geom_point(aes(color=factor(Gene.Order, levels = c('UF1','UF2','MF1')), shape=factor(Gene.Order, levels = c('UF1','UF2','MF1'))),size = 2.5, alpha = I(2/2)) +
  xlab("GC skew") +  ylab("GC content") +
  #geom_point(aes(color = Family, shape = Family), size = 2.5, alpha = I(2/2)) +
  geom_vline(aes(xintercept = mean(AT_skew)), color = "gray", linetype = "dashed") +
  geom_hline(aes(yintercept = mean(AT.content)), color = "gray", linetype = "dashed") +
  labs(color = 'Gene Order', shape = 'Gene Order') +
  scale_shape_manual(values = c(15,19,17)) +
  scale_color_manual(values = c('#00BFC4','#7CAE00','#F8766D'))
xlab("GC skew") +  ylab("GC content")

# Com subfamilias

ggplot(data=base_statistics, aes(x = GC_skew, y = GC.content)) +
  geom_point(aes(color=factor(Subfamily, levels = c('Anodontinae','Unioninae','Ambleminae','Gonideinae','Parreysiinae','Gibbosulinae','Margaritiferinae','Iridininae','Aspathariinae','Leilinae'))),size = 2.5, alpha = I(2/2), shape=factor(Family, levels = c('Mulleriidae','Margaritiferidae','Unionidae','Unionidae:Gonideinae','Iridinidae'))) +
  xlab("GC skew") +  ylab("GC content") +
  #geom_point(aes(color = Family, shape = Family), size = 2.5, alpha = I(2/2)) +
  geom_vline(aes(xintercept = mean(GC_skew)), color = "gray", linetype = "dashed") +
  geom_hline(aes(yintercept = mean(GC.content)), color = "gray", linetype = "dashed") +
  labs(color = 'Subfamily') +
  scale_shape_manual(values = c(19,17,15,23,3)) +
  #scale_color_manual(values = c('#00BFC4','#7CAE00','#F8766D'))
  xlab("GC skew") +  ylab("GC content")

ggplot(data=base_statistics, aes(x = AT_skew, y = AT.content)) +
  geom_point(aes(color=factor(Gene.Order, levels = c('UF1','UF2','MF1')), shape=factor(Gene.Order, levels = c('UF1','UF2','MF1'))),size = 2.5, alpha = I(2/2)) +
  xlab("GC skew") +  ylab("GC content") +
  #geom_point(aes(color = Family, shape = Family), size = 2.5, alpha = I(2/2)) +
  geom_vline(aes(xintercept = mean(AT_skew)), color = "gray", linetype = "dashed") +
  geom_hline(aes(yintercept = mean(AT.content)), color = "gray", linetype = "dashed") +
  labs(color = 'Gene Order', shape = 'Gene Order') +
  scale_shape_manual(values = c(15,19,17)) +
  scale_color_manual(values = c('#00BFC4','#7CAE00','#F8766D'))
xlab("GC skew") +  ylab("GC content")

###### Frequency of start and stop codons #####
#Start codons
start_freq <- read.csv('gene_start.xlsx',row.names=1)
cds_freq <- start_freq[1:13]
cds_freq[cds_freq == ""] = NA
cds_freq[4,] = NA
cds_freq[149,] <- NA
cds_freq <- cds_freq[complete.cases(cds_freq), ]

cds_freq_t <- t(cds_freq)
cds_freq_t <- data.frame(cds_freq_t)
cds <- apply(cds_freq_t,MARGIN=1,table)

library (plyr)
cds <- ldply(cds, data.frame) #CONVERT THE LIST INTO DATA FRAME

ggplot(cds, aes(x=.id, y =Freq , fill = factor(Var1))) +
  geom_bar(stat="identity",position='stack', color = '#505050') +
  ggtitle('Start Codon Frequency per Gene') +
  labs(x = 'Genes', y = 'Frequency (n)',fill = 'Start Codons') +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=0, hjust=0.5)) +
  scale_fill_manual(values=c("#d3a7cb", "#adc3f4", "#b5dfaf", "#f7eb89", "#d9d9d9", "#ffcf8d", "#e19b93", "#c3c3c3", "#aaaaaa", "#939393", "#797979", "#636363", "#515151", "#404040"))
#scale_fill_manual(values=c("#353e7c", "#007094", "#0eb6af", "#00be7d", "#5087C1", "#96d84b", "#fde333", "#61C074", "#a4d764", "#CCCCCC", "#Af5baa", "#CC78AF", "#CBA079", "#FFFFFF"))

#Stop codons
stop_freq <- read.csv('gene_stop.xlsx',row.names=1)
cds_freq <- stop_freq[1:13]
cds_freq[cds_freq == ""] = NA
cds_freq[4,] = NA
cds_freq[149,] <- NA
cds_freq[114,] <- NA

cds_freq <- cds_freq[complete.cases(cds_freq), ]

cds_freq_t <- t(cds_freq)
cds_freq_t <- data.frame(cds_freq_t)
cds <- apply(cds_freq_t,MARGIN=1,table)

library (plyr)
cds <- ldply(cds, data.frame) #CONVERT THE LIST INTO DATA FRAME

ggplot(cds, aes(x=.id, y =Freq , fill = factor(Var1))) +
  geom_bar(stat="identity",position='stack', color = '#505050') +
  ggtitle('Stop Codon Frequency per Gene') +
  labs(x = 'Genes', y = 'Frequency (n)',fill = 'Stop Codons') +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=0, hjust=0.5)) +
  scale_fill_manual(values = c('#ffab78','#ede863','#f0f0f0','#e5e5e5','#d9d9d9','#cccccc','#c0c0c0','#b4b4b4','#a8a8a8','#9c9c9c','#8f8f8f','#828282','#767676','#6a6a6a','#5e5e5e','#535353'))

####### Base frequencies per gene #####
library(dplyr)

#A
base_a <- read.csv('gene_A.xlsx', row.names = 1)

base_a[,38] <- NULL
base_a[base_a == ""] = NA

base_a <- base_a[complete.cases(base_a), ]

base_a <- t(base_a)
base_a <- data.frame(base_a)
base_a <- tibble::rownames_to_column(base_a, "gene")
base_a$base <- 'A'

base_a$values <- rowMeans(base_a[,2:119])
base_a$values
base_a[,2:119] <- NULL


#T
base_t <- read.csv('gene_T.xlsx', row.names = 1)

base_t[,38] <- NULL
base_t[base_t == ""] = NA

base_t <- base_t[complete.cases(base_t), ]
base_t <- t(base_t)
base_t <- data.frame(base_t)
base_t <- tibble::rownames_to_column(base_t, "gene")
base_t$base <- 'T'

base_t$values <- rowMeans(base_t[,2:119])
base_t[,2:119] <- NULL

#G
base_g <- read.csv('gene_G.xlsx', row.names = 1)

base_g[,38] <- NULL
base_g[base_g == ""] = NA

base_g <- base_g[complete.cases(base_g), ]
base_g <- t(base_g)
base_g <- data.frame(base_g)
base_g <- tibble::rownames_to_column(base_g, "gene")
base_g$base <- 'G'

base_g$values <- rowMeans(base_g[,2:119])
base_g[,2:119] <- NULL


#C
base_c <- read.csv('gene_C.xlsx', row.names = 1)

base_c[,38] <- NULL
base_c[base_c == ""] = NA

base_c <- base_c[complete.cases(base_c), ]
base_c <- t(base_c)
base_c <- data.frame(base_c)
base_c <- tibble::rownames_to_column(base_c, "gene")
base_c$base <- 'C'

base_c$values <- rowMeans(base_c[,2:119])
base_c[,2:119] <- NULL

#Concatenation

all_bases <- rbind(base_a,base_t,base_g,base_c)

#tirar rownames e converte-lo na primeira coluna
all_bases <- data.frame(t(all_bases))



ggplot(all_bases, aes(fill=base, y=values, x=gene)) + 
  geom_bar(position="stack", stat="identity")

#Grafico com todos os genes
ggplot(all_bases, aes(x=gene, y =values , fill = factor(base, levels = c('A','T','G','C')))) +
  geom_bar(stat="identity",position='stack', color = 'black') +
  scale_fill_manual(values=c("#fcfc4b", "#3ace3a", "#80CEE1","#ff392e")) +
  ggtitle('Base Frequency per Gene') +
  labs(x = 'Gene', y = 'Frequency (%)',fill = 'Base') +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=0, hjust=0.5))

ggplot(all_bases, aes(x=gene, y =values , fill = factor(base, levels = c('A','T','G','C')))) +
  geom_bar(stat="identity",position='stack', color = 'black') +
  scale_fill_manual(values=c("#fcfc4b", "#3ace3a", "#80CEE1","#ff392e")) +
  ggtitle('Base Frequency per Gene') +
  labs(x = 'Gene', y = 'Frequency (%)',fill = 'Base') +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=0, hjust=0.5))

####### Base frequencies per coding + rRNA #####

library(dplyr)
#A
base_a <- read.csv('gene_A.xlsx', row.names = 1)

base_a[,38] <- NULL
base_a[base_a == ""] = NA

base_a <- base_a[complete.cases(base_a), ]
base_a <- base_a[,1:13]

base_a <- t(base_a)
base_a <- data.frame(base_a)
base_a <- tibble::rownames_to_column(base_a, "gene")
base_a$base <- 'A'

base_a$values <- rowMeans(base_a[,2:119])
base_a$sd <- apply(base_a[,2:119],1,sd)


#T
base_t <- read.csv('gene_T.xlsx', row.names = 1)

base_t[,38] <- NULL
base_t[base_t == ""] = NA

base_t <- base_t[complete.cases(base_t), ]
base_t <- base_t[,1:13]

base_t <- t(base_t)
base_t <- data.frame(base_t)
base_t <- tibble::rownames_to_column(base_t, "gene")
base_t$base <- 'T'

base_t$values <- rowMeans(base_t[,2:119])
base_t$sd <- apply(base_t[,2:119],1,sd)

#G
base_g <- read.csv('gene_G.xlsx', row.names = 1)

base_g[,38] <- NULL
base_g[base_g == ""] = NA

base_g <- base_g[complete.cases(base_g), ]
base_g <- base_g[,1:13]

base_g <- t(base_g)
base_g <- data.frame(base_g)
base_g <- tibble::rownames_to_column(base_g, "gene")
base_g$base <- 'G'

base_g$values <- rowMeans(base_g[,2:119])
base_g$sd <- apply(base_g[,2:119],1,sd)


#C
base_c <- read.csv('gene_C.xlsx', row.names = 1)

base_c[,38] <- NULL
base_c[base_c == ""] = NA

base_c <- base_c[complete.cases(base_c), ]
base_c <- base_c[,1:13]

base_c <- t(base_c)
base_c <- data.frame(base_c)
base_c <- tibble::rownames_to_column(base_c, "gene")
base_c$base <- 'C'

base_c$values <- rowMeans(base_c[,2:119])
base_c$sd <- apply(base_c[,2:119],1,sd)

#Concatenation

all_bases <- rbind(base_a,base_t,base_g,base_c)

#Grafico com so os cds
ggplot(all_bases, aes(x=gene, y=values , fill = factor(base, levels = c('A','T','G','C')))) +
  geom_bar(stat="identity",position=position_dodge(), color = '#505050') +
  scale_fill_manual(values=c("#f5d96e", "#99c796", "#aac5f2","#f29164")) +
  ggtitle('Base Frequency per Gene') +
  labs(x = 'Gene', y = 'Frequency (%)',fill = 'Base') +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=0, hjust=0.5)) +
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.2,
                position=position_dodge(.9))


##### Gene Length ####

gene_len <- read.csv('gene_length.xlsx', row.names = 1)
View(gene_len)
gene_len[,38] <- NULL
gene_len <- gene_len[1:13]
gene_len[gene_len == ""] = NA

gene_len <- gene_len[complete.cases(gene_len), ]

gene_len = t(gene_len)
gene_len <- data.frame(gene_len)
gene_len <- tibble::rownames_to_column(gene_len, "gene")
gene_len$values <- rowMeans(gene_len[,2:137])
gene_len$sd <- apply(gene_len[,2:137],1,sd)

ggplot(gene_len, aes(x=gene, y=values, color = gene)) +
  geom_bar(stat="identity",position=position_dodge(), fill = rgb(1,1,1)) +
  ggtitle('Gene Length') +
  labs(x = 'Gene', y = 'Gene Length') +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.2,
                position=position_dodge(.9))

ggplot(gene_len, aes(x=gene, y=values)) +
  geom_bar(stat="identity",position=position_dodge()) +
  ggtitle('Gene Length') +
  labs(x = 'Gene', y = 'Gene Length') +
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=0.5), legend.position="none") +
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.2,
                position=position_dodge(.9))


##### PCG GC #####

gene_gc <- read.csv('gene_gc.xlsx', row.names = 1)
gene_gc[,38] <- NULL
gene_gc <- gene_gc[,1:13]
gene_gc[gene_gc == ""] = NA

gene_gc <- gene_gc[complete.cases(gene_gc), ]

ggplot(stack(gene_gc), aes(x=ind, y=values, color = ind)) +
  geom_jitter(aes(x = ind, color = ind), 
              position = position_jitter(width = .05), alpha = 0.5) +
  geom_boxplot(outlier.colour = NA, position = "dodge") +
  ggtitle('Gene GC content') +
  labs(x = 'Gene', y = 'GC content (%)') +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=0, hjust=0.5),
        legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=23, size=4)

##### AT and GC skews per Gene #####

at_skew <- read.csv('gene_ATskew.xlsx',row.names=1)
at_skew <- at_skew[1:13]
at_skew[at_skew == ""] = NA
at_skew[4,] = NA
at_skew[149,] <- NA
at_skew <- at_skew[complete.cases(at_skew), ]
at_skew <- t(at_skew)
at_skew <- data.frame(at_skew)
at_skew <- tibble::rownames_to_column(at_skew, 'genes')
at_skew$values <- rowMeans(at_skew[,2:135])
at_skew$sd <- apply(at_skew[,2:135],1,sd)
at_skew$supp <- 'AT'


gc_skew <- read.csv('gene_GCskew.xlsx',row.names=1)
gc_skew <- gc_skew[1:13]
gc_skew[gc_skew == ""] = NA
gc_skew[4,] = NA
gc_skew[149,] <- NA
gc_skew <- gc_skew[complete.cases(gc_skew), ]
gc_skew <- t(gc_skew)
gc_skew <- data.frame(gc_skew)
gc_skew <- tibble::rownames_to_column(gc_skew, 'genes')
gc_skew$values <- rowMeans(gc_skew[,2:135])
gc_skew$sd <- apply(gc_skew[,2:135],1,sd)
gc_skew$supp <- 'GC'

all_skews <- rbind(at_skew,gc_skew)

# para todos os genes
ggplot(all_skews, aes(x=genes, y=values, group = supp, shape = factor(supp, levels = c('GC','AT')), color = factor(supp, levels = c('GC','AT')))) +
  geom_line(size = 1.1) +
  geom_point(size = 3.2) +
  #scale_fill_manual(values=c("#fcfc4b", "#3ace3a", "#80CEE1","#ff392e")) +
  ggtitle('AT and GC skews') +
  labs(x = 'Genes', y = 'AT/GC skew',color = 'Skews', shape = 'Skews') +
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.2,)#+
#theme(text = element_text(size=20),
#axis.text.x = element_text(angle=0, hjust=0.5)) +
#geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.2,
#position=position_dodge(.9))

all_skews_H <- rbind(all_skews[1:2,],all_skews[4:6,],all_skews[9:10,],all_skews[11:12,],all_skews[14:15,],all_skews[17:19,],all_skews[22:23,],all_skews[24:25,])
all_skews_L <- rbind(all_skews[3,],all_skews[7:8,],all_skews[13,],all_skews[16,],all_skews[20:21,],all_skews[26,])

# plot com os da H-strand
ggplot(all_skews_H, aes(x=genes, y=values, group = supp, shape = factor(supp, levels = c('GC','AT')), color = factor(supp, levels = c('GC','AT')))) +
  geom_line(size = 1.1) +
  geom_point(size = 3.2) +
  #scale_fill_manual(values=c("#fcfc4b", "#3ace3a", "#80CEE1","#ff392e")) +
  ggtitle('AT and GC skews') +
  labs(x = 'H-Strand Genes', y = 'AT/GC skew',color = 'Skews', shape = 'Skews') +
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.1,)

ggplot(all_skews_L, aes(x=genes, y=values, group = supp, shape = factor(supp, levels = c('GC','AT')), color = factor(supp, levels = c('GC','AT')))) +
  geom_line(size = 1.1) +
  geom_point(size = 3.2) +
  #scale_fill_manual(values=c("#fcfc4b", "#3ace3a", "#80CEE1","#ff392e")) +
  ggtitle('AT and GC skews') +
  labs(x = 'L-strand Genes', y = 'AT/GC skew',color = 'Skews', shape = 'Skews') +
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.1,)


