library(readxl)
setwd('/home/edu/Desktop/Mitogenomas/All_Mitogenomas/Mitogenome_Statistics')

atp6_kaks <- read.excel('atp6_kaks.xlsx')

atp6_kaks[atp6_kaks == 0.0000000] = NA

atp6.mean <- atp6_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(atp6.mean[2:119]))

atp8_kaks <- read.excel('atp8_kaks.xlsx')

atp8_kaks[atp8_kaks == 0.0000000] = NA
atp8_kaks[2,114] <- NA
atp8_kaks[2,115] <- NA
atp8_kaks[,16] <- NULL

atp8.mean <- atp8_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(atp8.mean[2:118]))


cox1_kaks[cox1_kaks == 0.0000000] = NA

cox1.mean <- cox1_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(cox1.mean[2:119]))


cob_kaks[cob_kaks == 0.0000000] = NA
cob_kaks[37,83] <- NA

cob.mean <- cob_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(cob.mean[2:119]))

cox2_kaks[cox2_kaks == 0.0000000] = NA
cox2_kaks[65,66] <- NA

cox2.mean <- cox2_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(cox2.mean[2:119]))


cox3_kaks[cox3_kaks == 0.0000000] = NA

cox3.mean <- cox3_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(cox3.mean[2:119]))


nad1_kaks[nad1_kaks == 0.0000000] = NA

nad1.mean <- nad1_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(nad1.mean[2:119]))



nad2_kaks[nad2_kaks == 0.0000000] = NA

nad2.mean <- nad2_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(nad2.mean[2:119]))


nad3_kaks[nad3_kaks == 0.0000000] = NA
nad3_kaks[37,84] <- NA

nad3.mean <- nad3_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(nad3.mean[2:119]))


nad4_kaks[nad4_kaks == 0.0000000] = NA

nad4.mean <- nad4_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(nad4.mean[2:119]))


nad4l_kaks[nad4l_kaks == 0.0000000] = NA

nad4l.mean <- nad4l_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(nad4l.mean[2:119]))


nad5_kaks[nad5_kaks == 0.0000000] = NA

nad5.mean <- nad5_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(nad5.mean[2:119]))


nad6_kaks[nad6_kaks == 0.0000000] = NA

nad6.mean <- nad6_kaks %>%
  summarise_all(mean, na.rm = TRUE)

mean(as.numeric(nad6.mean[2:119]))


kaks_names <- c('atp6','atp8','cob','cox1','cox2','cox3','nad1','nad2','nad3','nad4','nad4l','nad5','nad6') 
kaks_values <- c(mean(as.numeric(atp6.mean[2:119])),mean(as.numeric(atp8.mean[2:118])),mean(as.numeric(cob.mean[2:119])),mean(as.numeric(cox1.mean[2:119])),mean(as.numeric(cox2.mean[2:119])),mean(as.numeric(cox3.mean[2:119])),mean(as.numeric(nad1.mean[2:119])),mean(as.numeric(nad2.mean[2:119])),mean(as.numeric(nad3.mean[2:119])),mean(as.numeric(nad4.mean[2:119])),mean(as.numeric(nad4l.mean[2:119])),mean(as.numeric(nad5.mean[2:119])),mean(as.numeric(nad6.mean[2:119])))

kaks <- data.frame(name = c('atp6','atp8','cob','cox1','cox2','cox3','nad1','nad2','nad3','nad4','nad4l','nad5','nad6') ,
                   values = c(mean(as.numeric(atp6.mean[2:119])),mean(as.numeric(atp8.mean[2:118])),mean(as.numeric(cob.mean[2:119])),mean(as.numeric(cox1.mean[2:119])),mean(as.numeric(cox2.mean[2:119])),mean(as.numeric(cox3.mean[2:119])),mean(as.numeric(nad1.mean[2:119])),mean(as.numeric(nad2.mean[2:119])),mean(as.numeric(nad3.mean[2:119])),mean(as.numeric(nad4.mean[2:119])),mean(as.numeric(nad4l.mean[2:119])),mean(as.numeric(nad5.mean[2:119])),mean(as.numeric(nad6.mean[2:119]))),
                   stand = c(mean(as.numeric(atp6.sd[3:119])),mean(as.numeric(atp8.sd[3:118])),mean(as.numeric(cob.sd[3:119])),mean(as.numeric(cox1.sd[3:119])),mean(as.numeric(cox2.sd[3:119])),mean(as.numeric(cox3.sd[3:119])),mean(as.numeric(nad1.sd[3:119])),mean(as.numeric(nad2.sd[3:119])),mean(as.numeric(nad3.sd[3:119]), na.rm = TRUE),mean(as.numeric(nad4.sd[3:119])),mean(as.numeric(nad4l.sd[3:119])),mean(as.numeric(nad5.sd[3:119])),mean(as.numeric(nad6.sd[3:119]))))


ggplot(kaks, aes(x=name, y =values)) +
  geom_bar(stat="identity", fill= 'steelblue', color = '#505050', width = 0.85) +
  ggtitle('Ka/Ks') +
  labs(x = 'Genes', y = 'Ka/Ks value') +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=0, hjust=0.5)) +
  geom_errorbar(aes(ymin=values-stand/2.5, ymax=values+stand/2.5), width=.2,
                position=position_dodge(.9), color = '#707070')

atp6.sd <- atp6_kaks %>%
  summarise_all(sd, na.rm = TRUE)

atp8.sd <- atp8_kaks %>%
  summarise_all(sd, na.rm = TRUE)

cox1.sd <- cox1_kaks %>%
  summarise_all(sd, na.rm = TRUE)

cox2.sd <- cox2_kaks %>%
  summarise_all(sd, na.rm = TRUE)

cox3.sd <- cox3_kaks %>%
  summarise_all(sd, na.rm = TRUE)

nad1.sd <- nad1_kaks %>%
  summarise_all(sd, na.rm = TRUE)

nad3.sd <- nad3_kaks %>%
  summarise_all(sd, na.rm = TRUE)

nad3.sd[16] <- NULL

nad2.sd <- nad2_kaks %>%
  summarise_all(sd, na.rm = TRUE)

nad4.sd <- nad4_kaks %>%
  summarise_all(sd, na.rm = TRUE)

nad4l.sd <- nad4l_kaks %>%
  summarise_all(sd, na.rm = TRUE)

nad5.sd <- nad5_kaks %>%
  summarise_all(sd, na.rm = TRUE)

nad6.sd <- nad6_kaks %>%
  summarise_all(sd, na.rm = TRUE)

all_skews_H <- rbind(all_skews[1:2,],all_skews[4:6,],all_skews[9:10,],all_skews[11:12,],all_skews[14:15,],all_skews[17:19,],all_skews[22:23,],all_skews[24:25,])

mycetopodidae <- atp6_kaks[18,115]
margaritiferidae <- mean(as.numeric(atp6_kaks[31,49],atp6_kaks[31,50]))
iridinidae <- mean(as.numeric(atp6_kaks[57,58],atp6_kaks[57,60],atp6_kaks[58,60]), na.rm = TRUE)


mycetopodidae <- data.frame(name = c('atp6','atp8','cob','cox1','cox2','cox3','nad1','nad2','nad3','nad4','nad4l','nad5','nad6') ,
                   values = c(atp6_kaks[18,115],atp8_kaks[18,115],cob_kaks[18,115],cox1_kaks[18,115],cox2_kaks[18,115],cox3_kaks[18,115],nad1_kaks[18,115],nad2_kaks[18,115],nad3_kaks[18,115],nad4_kaks[18,115],nad4l_kaks[18,115],nad5_kaks[18,115],nad6_kaks[18,115])
                   
mycetopodidae <- data.frame(name = c('atp6','atp8','cob','cox1','cox2','cox3','nad1','nad2','nad3','nad4','nad4l','nad5','nad6') ,
                   values = c(mean(as.numeric(atp6.mean[2:119])),mean(as.numeric(atp8.mean[2:118])),mean(as.numeric(cob.mean[2:119])),mean(as.numeric(cox1.mean[2:119])),mean(as.numeric(cox2.mean[2:119])),mean(as.numeric(cox3.mean[2:119])),mean(as.numeric(nad1.mean[2:119])),mean(as.numeric(nad2.mean[2:119])),mean(as.numeric(nad3.mean[2:119])),mean(as.numeric(nad4.mean[2:119])),mean(as.numeric(nad4l.mean[2:119])),mean(as.numeric(nad5.mean[2:119])),mean(as.numeric(nad6.mean[2:119])))