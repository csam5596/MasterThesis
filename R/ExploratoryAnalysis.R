library(ggplot2)
library(ggpubr)

data <- read.csv("data/Immunogenicity_Dataset.csv")
data$Peptide <- as.character(data$Peptide)
data$length <- nchar(data$Peptide)

table(data$Immunogenicity)

peptides <- data.frame(
  "peptide" = data$Peptide,
  "immunogenicity" = data$Immunogenicity
)

##############################
#### EXPLORATORY ANALYSIS ####
##############################

## Binding Affinity ##
binding <- read.csv("data/Immunogenicity_Dataset_Binding.csv")
binding$Strength <- ifelse(binding$Rank > 2, "NB", ifelse(binding$Rank < 0.5, "SB", "WB"))
binding$Strength <- factor(binding$Strength, levels = c("SB", "WB", "NB"))


gg.ecdf <- ggplot(binding, aes(Rank)) +
  stat_ecdf(geom = "step") + 
  ggtitle("Empirical CDF") + 
  xlab("%Rank") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"),
        title = element_text(size = 20, face ="bold"),
        plot.margin = unit(c(0,1,0,0), "lines"))

gg.freq <- ggplot(binding, aes(x = Strength, group = Strength)) +
  geom_bar(aes(y = ..count.., fill = "#003361")) + 
  geom_text(aes( label = scales::percent(..count../sum(..count..)),
                 y= ..count.. ), stat= "count", vjust = -0.3, size = 5) +
  labs(y = "Count", fill="#003361") +
  theme(legend.position = "none", 
        axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        title = element_text(size = 20, face ="bold"),
        plot.margin = unit(c(0,0,0,1), "lines")) +
  ggtitle("Binding Strength Distribution") +
  scale_fill_manual(values = c("#003361", "#003361", "#003361"))

ggarrange(gg.ecdf, gg.freq)

## Distribution ##
ggplot(data, aes(x= Immunogenicity,  group=length)) + 
  geom_bar(aes(y = ..count.., fill = factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..count.. ), stat= "count", vjust = -0.3, size = 5) +
  labs(y = "Count", fill="Immunogenicity") +
  facet_grid(~length) +
  theme(legend.position = "none", axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        title = element_text(size = 20, face ="bold"),
        strip.text.x = element_text(size = 14)) +
  ggtitle("Data Distribution") +
  scale_fill_manual(values = c("#003361", "#f39200"))

addmargins(table(data$Immunogenicity, data$length))
