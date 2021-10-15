# ML for Chondrichthyes evolutionary pattern identification

# import necessary functions
source("/home/edu/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/chondrichthyes/functions.R")

# Set the right working directory
setwd("~/Desktop/Bioinformatica/Mitogenomics/Chondrichthyes/Deep_Mitogenomics")

# Import necessary library to open csv file
library(readr)

# Read the csv file
dataset <- read_csv("Chondrichthyes.csv")

dataset = subset(dataset, dataset$State == "complete")

codon_dataset <- read.csv("Chondrichthyes_codons_.csv")

codon_dataset = subset(codon_dataset, codon_dataset$State == "complete")


complete_mitogenomes = subset(dataset, dataset$Content == "Mitogenome")
complete_mitogenomes = complete_mitogenomes[-c(4:7)]
complete_mitogenomes = complete_mitogenomes[-c(12:25)]
# Creating a vector to isolate protein coding genes 
protein_coding <- c('atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6')

# Isolating all protein coding genes
PCG <- subset(dataset, dataset$Content %in% protein_coding)

# merge the 2 datasets

PCG[,34:93] <- codon_dataset[,7:66]


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

check_missing_genes(data.frame(PCG)) # "KT698052.1" "KC992792.1"

taxonomy_vector_3 <- c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes','Batoidea incertae sedis') # Orders
PCG <- get_taxonomy_column(taxonomy_vector_3, PCG, 'Order')

# Feature selection

# clean up columns that don't matter

PCG <- PCG[-c(1:2,4:7)]

nrow(subset(PCG_clean, PCG_clean$Order == "Myliobatiformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Rajiformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Torpediniformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Pristiophoriformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Carcharhiniformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Lamniformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Chimaeriformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Orectolobiformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Rhinopristiformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Hexanchiformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Heterodontiformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Squaliformes"))
nrow(subset(PCG_clean, PCG_clean$Order == "Squatiniformes"))

PCG_test <- subset(PCG_clean, PCG_clean$Order == "Myliobatiformes" | PCG_clean$Order == "Rajiformes" | PCG_clean$Order == "Carcharhiniformes" | PCG_clean$Order == "Lamniformes")

library(caret)

library(Boruta)

# Perform Boruta search

boruta_output <- Boruta(factor(Order) ~ ., data = na.omit(PCG_clean), doTrace = 0)

PCG_dclean = PCG_clean[-c(9:10)]

PCG_scaled <- scale(PCG_dclean[1:80])

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

library(dplyr)




#PCG_grouped <- PCG_dclean %>%
#    group_by(Accession) %>%
#    for (i in seq_along(8:80)) {
#        mutate()
#    }



############### New

library(readr)

# Read the csv file
dataset <- read_csv("Chondrichthyes.csv")

dataset = subset(dataset, dataset$State == "complete")

codon_dataset <- read.csv("Chondrichthyes_codons.csv")

codon_dataset = subset(codon_dataset, codon_dataset$State == "complete")


complete_mitogenomes = subset(dataset, dataset$Content == "Mitogenome")
complete_mitogenomes = complete_mitogenomes[-c(4:7)]
complete_mitogenomes = complete_mitogenomes[-c(12:25)]
# Creating a vector to isolate protein coding genes 
protein_coding <- c('atp6','atp8','cox1','cox2','cox3','cob','nad1','nad2','nad3','nad4','nad4l','nad5','nad6')

# Isolating all protein coding genes
PCG <- subset(dataset, dataset$Content %in% protein_coding)

# merge the 2 datasets

PCG[,30:89] <- codon_dataset[,7:66]


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

check_missing_genes(data.frame(PCG)) # "KT698052.1" "KC992792.1"

taxonomy_vector_3 <- c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes','Carcharhiniformes','Lamniformes','Orectolobiformes','Heterodontiformes','Hexanchiformes','Pristiophoriformes','Squaliformes','Squatiniformes','Chimaeriformes','Batoidea incertae sedis') # Orders
PCG <- get_taxonomy_column(taxonomy_vector_3, PCG, 'Order')

# Feature selection

library(mltools)
library(data.table)
library(caret)

# clean up columns that don't matter

PCG_clean <- PCG[-c(1:2)]
PCG_clean <- PCG_clean[-c(2:5)]

# One hot encoding categorical data

dummy <- dummyVars(" ~ Content + Start + Stop", data=PCG_clean)
newdata <- data.frame(predict(dummy, newdata = PCG_clean))

PCG_clean[,30:52] <- newdata[,1:23]

PCG_clean <- PCG_clean[-c(1,10:11)]


library(caret)

library(Boruta)

# Perform Boruta search

boruta_output <- Boruta(factor(Order) ~ ., data = na.omit(PCG_clean), doTrace = 0)

print(boruta_output)

boruta_signif <- getSelectedAttributes(boruta_output, withTentative = FALSE)
print(boruta_signif)

plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")  

PCG_clean$GRAVY <- gsub(",", "", as.character(PCG_clean$GRAVY)) # erro no script de python, os gravy ficam com virgula! corrigir later
PCG_clean$GRAVY <- as.numeric(PCG_clean$GRAVY)

# Random Forest Variable Importance

library(e1071)
rrfMod <- train(Order ~ ., data = na.omit(PCG_clean), method = "RRF")
rrfImp <- varImp(rrfMod, scale = T)

rrfImp
plot(rrfImp, top = 15, main='Variable Importance')


PCG_dclean = PCG_clean[-c(9:10)]

PCG_scaled <- scale(PCG_dclean[1:80])

corr.test(PCG_clean$GC_skew, PCG_clean$GC_content, method="spearman", adjust = "bonferroni")


cox1$GRAVY <- gsub(",", "", as.character(cox1$GRAVY)) # erro no script de python, os gravy ficam com virgula! corrigir later
cox1$GRAVY <- as.numeric(cox1$GRAVY)

cox1[,1:7] <- NULL
cox1[,9:10] <- NULL

library(e1071)
rrfMod <- train(Order ~ ., data = na.omit(cox1), method = "RRF")
rrfImp <- varImp(rrfMod, scale = T)

rrfImp
plot(rrfImp, top = 15, main='Variable Importance')

cox2[,1:7] <- NULL
cox2[,9:10] <- NULL

library(e1071)
rrfMod <- train(Order ~ ., data = na.omit(cox2), method = "RRF")
rrfImp <- varImp(rrfMod, scale = T)

rrfImp
plot(rrfImp, top = 15, main='Variable Importance')


cox3[,1:7] <- NULL
cox3[,9:10] <- NULL

library(e1071)
rrfMod <- train(Order ~ ., data = na.omit(cox3), method = "RRF")
rrfImp <- varImp(rrfMod, scale = T)

rrfImp
plot(rrfImp, top = 15, main='Variable Importance')



atp8[,1:7] <- NULL
atp8[,9:10] <- NULL

library(e1071)
rrfMod <- train(Order ~ ., data = na.omit(atp8), method = "RRF")
rrfImp <- varImp(rrfMod, scale = T)

rrfImp
plot(rrfImp, top = 15, main='Variable Importance')



atp6[,1:7] <- NULL
atp6[,9:10] <- NULL

library(e1071)
rrfMod <- train(Order ~ ., data = na.omit(atp6), method = "RRF")
rrfImp <- varImp(rrfMod, scale = T)

rrfImp
plot(rrfImp, top = 15, main='Variable Importance')



nad1[,1:7] <- NULL
nad1[,9:10] <- NULL

library(e1071)
rrfMod <- train(Order ~ ., data = na.omit(nad1), method = "RRF")
rrfImp <- varImp(rrfMod, scale = T)

rrfImp
plot(rrfImp, top = 15, main='Variable Importance')



## PCA with all data

PCG_2 <- na.omit(PCG)
nad5 <- na.omit(nad5)

pca <- prcomp(nad5[-c(1,10:11,88)], scale. = TRUE)
summary(pca)

fviz_eig(pca)

library(devtools)
library(ggbiplot) 
library(factoextra)


p_pca <- fviz_pca_ind(pca,
                         axes = c(1,2),
                         col.ind =  factor(na.omit(nad5$Order), levels = c('Myliobatiformes','Rajiformes','Rhinopristiformes','Torpediniformes',
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

p_pca + scale_color_manual(name = "Order", values = c("#f50000", "#cc0000", "#a30000","#4e1609",'#4f86f7','#45b1e8','#87ceeb','#30d5c8','#adff2f','#32cd32',
                                                         '#3cb371','#006400','#7b68ee')) + scale_shape_manual(values = c(rep(16, 13))) +
    guides(fill = "none" , shape = "none") 

fviz_pca_var(pca, col.var="steelblue")


### Support Vector Machines

library(e1071)
library(caTools)

# Splitting the data

PCG_ML <- subset(PCG_ML, PCG_ML$Order != "Pristiophoriformes")
PCG_ML[,109] <- PCG_ML[,85]
PCG_ML$Order <- NULL
PCG_ML[,109] <- PCG[,88]
colnames(PCG_ML)[108] <- "Order"

PCG_ML$Order <- as.factor(PCG_ML$Order)

split = sample.split(PCG_ML$Order, SplitRatio = 0.8)
training_set = subset(PCG_ML, split == TRUE)
test_set = subset(PCG_ML, split == FALSE)

#only PCG with RSCU data

PCG_ML <- na.omit(PCG_ML)

PCG_RSCU <- PCG_ML[-c(1:25,85:107)]
PCG_woRSCU <- PCG_ML[-c(26:84)]

split = sample.split(PCG_woRSCU$Order, SplitRatio = 0.8)
training_set = subset(PCG_woRSCU, split == TRUE)
test_set = subset(PCG_woRSCU, split == FALSE)

# Feature scale!

training_set[-c(26:49)] <- scale(training_set[-c(26:49)])
test_set[-c(26:49)] <- scale(test_set[-c(26:49)])

gammalist <- c(0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05)
PCG_ML$StartATT <- NULL
svmfit <- tune.svm(as.factor(Order) ~., data = training_set, type = "C-classification", kernel = "radial", cost = 2^(-1:5), gamma = gammalist, scale = FALSE)


svmfit_poly <- tune.svm(as.factor(Order) ~., data = training_set, type = "C-classification", kernel = "polynomial", degree = 2, cost = 2^(-1:5), gamma = gammalist, scale = FALSE)






summary(svmfit)
summary(svmfit$best.model)

test_set <- na.omit(test_set)
svm <- predict(svmfit$best.model, na.omit(test_set[,1:48]))
confusionMatrix(svm, as.factor(na.omit(test_set[,49])))

svm3 <- predict(svmfit_poly$best.model, na.omit(test_set[,1:104]))
confusionMatrix(svm3, as.factor(na.omit(test_set[1:547,105])))

svm4 <- predict(svmfit$best.model, na.omit(test_set[,1:72]))
confusionMatrix(svm4, as.factor(na.omit(test_set[,73])))

svmfit_poly$train.ind

## training with only RSCU data






y_pred <- predict(svmfit, newdata = test_set[,1:107])

dummy <- dummyVars(" ~ Content + Start + Stop", data=PCG)
newdata <- data.frame(predict(dummy, newdata = PCG))

PCG_ML <- as.data.frame(cbind(PCG, newdata))

PCG_ML <- PCG_ML[-c(1,10:11)]
cm = table(test_set$Order, y_pred)
