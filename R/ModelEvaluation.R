library(pROC)
library(ggplot2)
library(reshape)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(bio3d)
library(randomForest)

source("Blomap.R")

aucs <- data.frame(
  "Factor" = rep(NA, 100),
  "Hydro" = NA,
  "Blomap" = NA,
  "AAC" = NA,
  "QSO" = NA,
  "PAAC" = NA,
  "APAAC" = NA,
  "Blosum62" = NA,
  "Blosum35" = NA,
  "Wen" = NA
)


for(i in 1:100){
  
  table.name <- paste("results/iteration_", i, ".txt", sep="")
  d <- read.csv(table.name, sep=" ")
  
  ## AUCS ##
  aucs$Factor[i] <- auc(d$Immunogenicity, d$Factor)
  aucs$Hydro[i] <- auc(d$Immunogenicity, d$Hydro)
  aucs$Blomap[i] <- auc(d$Immunogenicity, d$Blomap)
  aucs$AAC[i] <- auc(d$Immunogenicity, d$AAC)
  aucs$QSO[i] <- auc(d$Immunogenicity, d$QSO)
  aucs$PAAC[i] <- auc(d$Immunogenicity, d$PAAC)
  aucs$APAAC[i] <- auc(d$Immunogenicity, d$APAAC)
  aucs$Blosum62[i] <- auc(d$Immunogenicity, d$Blosum62)
  aucs$Blosum35[i] <- auc(d$Immunogenicity, d$Blosum35)
  aucs$Wen[i] <- auc(d$Immunogenicity, d$Wen)
}

colnames(aucs) <- c("Factor", "Hydro", "BLOMAP", "AAC", "QSO", "PseAAC", "APseAAC", "Blosum62", "Blosum35", "Wen")

############
### PLOT ###
############

aucsLong <- melt(aucs)

#### DISTRIBUTION OF AUCS ####

ggplot(aes(x=variable, fill=alpha("#003361", 0.5)), data=aucsLong) + 
  geom_boxplot(aes(y = value)) +
  ggtitle("Bootstrap Results", subtitle = "AUC by Model") +
  xlab("Model") + 
  ylab("AUC") +
  scale_fill_manual(values=alpha(rep("#003361",10), 0.5)) +
  theme(legend.position = "none", axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        title = element_text(size = 20, face ="bold"),
        plot.subtitle =element_text(size=18, face="italic"))


#############
### TESTS ###
#############

aucs$id <- 1:100

aucsTest <- data.frame(
  "id" = 1:100,
  "Factor" = aucs$Factor,
  "Hydro" = aucs$Hydro,
  "BLOMAP" = aucs$BLOMAP
)


aucsFull <- aucs %>%
  gather(key = "model", value = "AUC", Factor, Hydro, BLOMAP, AAC, QSO, PseAAC, APseAAC, Blosum62, Blosum35, Wen) %>%
  convert_as_factor(id, model)

aucsFull %>% group_by(model) %>% summarize(n = n(),
                                           min = fivenum(AUC)[1], 
                                           Q1 = fivenum(AUC)[2],
                                           median = fivenum(AUC)[3],
                                           Q3 = fivenum(AUC)[4],
                                           max = fivenum(AUC)[5])

aucsTest <- aucsTest %>%
  gather(key = "model", value = "AUC", Factor, Hydro, BLOMAP) %>%
  convert_as_factor(id, model)

### NORMALITY ###
aucsTest %>% group_by(model) %>% shapiro_test(AUC)

ggqqplot(aucsTest, "AUC", facet.by = "model")

### Repeated Measures ANOVA ###

ggplot(aes(x=variable, fill=alpha("#003361", 0.5)), data=aucsLong) + 
  geom_boxplot(aes(y = value)) +
  ggtitle("Bootstrap Results", subtitle = "AUC by Model") +
  xlab("Model") + 
  ylab("AUC") +
  scale_fill_manual(values=alpha(rep("#003361",10), 0.5)) +
  theme(legend.position = "none",
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size = 20, face ="bold"),
        plot.subtitle =element_text(size=16, face="italic"),
        axis.text = element_text(size = 14))

aucs.aov <- anova_test(data = aucsTest, dv = AUC, wid = id, within = model)
get_anova_table(aucs.aov)

### Post Hoc Test ###
pwc <- aucsTest %>%
  pairwise_t_test(
    AUC ~ model, paired = TRUE,
    p.adjust.method = "bonferroni"
  )

pwc <- pwc %>% add_xy_position(x = "model")
plt <- ggboxplot(aucsTest, x = "model", y = "AUC", add = "point", fill=alpha("#003361", 0.5)) + 
  theme_grey() + 
  stat_pvalue_manual(pwc) +
  xlab("") +
  labs(
    subtitle = get_test_label(aucs.aov, detailed = TRUE, correction = "GG"),
    caption = get_pwc_label(pwc)
  ) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 16, face = "italic"),
        plot.caption = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "none") + 
  ggtitle("Repeated Measures ANOVA") +
  xlab("Model")

plt

##########
### RF ###
##########

data <- read.csv("data/Immunogenicity_Dataset.csv")

peptides <- data.frame(
  "MutatedPeptide" = data$Peptide,
  "Immunogenicity" = data$Immunogenicity
)
peptides$MutatedPeptide <- as.character(peptides$MutatedPeptide)

#### FACTOR ####
peptidesFactor <- do.call(rbind.data.frame, lapply(peptides$MutatedPeptide, function(x){
  strsplit(str_pad(x, 11, side = "right", pad = 'X'), split="")[[1]]
}))
colnames(peptidesFactor) <- paste("P", 1:11, sep="_")
peptidesFactor$Immunogenicity <- peptides$Immunogenicity

#### HYDROPHOBICITY / BULKINESS / POLARITY ####
pepPhys <- do.call(rbind.data.frame, lapply(peptides$MutatedPeptide, function(x){
  aas <- strsplit(str_pad(x, 11, side = "right", pad = 'X'), split="")[[1]]
  as.vector(rbind(
    "hydrophobicity" = aa2index(aas, index = "KYTJ820101", window=1),
    "bulkiness" = aa2index(aas, index = "ZIMJ680102", window = 1),
    "polarity" = aa2index(aas, index = "GRAR740102", window = 1)
  ))
}))
pepPhys[is.na(pepPhys)] <- 0

colnames(pepPhys) <- as.vector(t(outer(
  c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11"),
  c("Hydro", "Bulk", "Polarity"),
  paste, sep="_")))
pepPhys$Immunogenicity <- peptides$Immunogenicity

#### BLOMAP #####
pepBlomap <- do.call(rbind.data.frame, lapply(peptides$MutatedPeptide, function(x){
  aas <- strsplit(str_pad(x, 11, side = "right", pad = 'X'), split="")[[1]]
  as.vector(sapply(aas, blomap))
}))
colnames(pepBlomap) <- as.vector(t(outer(
  c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11"),
  c("X1", "X2", "X3", "X4", "X5"),
  paste, sep="_")))
pepBlomap$Immunogenicity <- peptides$Immunogenicity

set.seed(01015742)
rfF <- randomForest(Immunogenicity ~ ., data=peptidesFactor)
rfBm <- randomForest(Immunogenicity ~ ., data=pepBlomap)
rfH <- randomForest(Immunogenicity ~ ., data=pepPhys)

importance.f <- data.frame(importance(rfF))
importance.f$grp <- 1:11
importance.f$n <- rep(1,11)

importance.bm <- data.frame(importance(rfBm))
importance.bm$gr <- rep(1:11, each=5)
importance.bm$n <- rep(rep(1:5),11)

importance.H <- data.frame(importance(rfH))
importance.H$gr <- rep(1:11, each = 3)
importance.H$n <- rep(c("Hydro", "Bulk", "Polarity"), 11)

gg.f <- ggplot(importance.f, aes(y = MeanDecreaseGini, x = n)) +
  geom_bar(stat = "identity", fill = "#003361") +
  facet_grid(cols = vars(grp)) +
  theme(axis.title = element_text(size = 18, face="bold"),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = "#003361") +
  ggtitle("Factor Encoding")

gg.bm <- ggplot(importance.bm, aes(y = MeanDecreaseGini, x = n)) +
  geom_bar(stat = "identity", fill = "#003361") +
  facet_grid(cols = vars(gr)) +
  theme(axis.title = element_text(size = 18, face="bold"),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = "#003361") +
  ggtitle("BLOMAP Encoding")


gg.H <- ggplot(importance.H, aes(y = MeanDecreaseGini, x = n)) +
  geom_bar(stat = "identity", fill = "#003361") +
  facet_grid(cols = vars(gr)) +
  theme(axis.title = element_text(size = 18, face="bold"),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = "#003361") +
  ggtitle("Chemical Properties Encoding (Hydro)")


ggarrange(gg.f, gg.bm, gg.H, nrow = 3)
