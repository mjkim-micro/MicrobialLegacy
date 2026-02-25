#### clinical data #####
### load data and packages ####
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)
library(ggpubr)
library(microbiome)
library(reshape2)
library(knitr)
library(cowplot)
library(ggsignif)
library(gridExtra)
library(car)
library(agricolae)
library(Maaslin2)
library(compositions)
library(pROC)
library(gplots)

### Maternal ####

### Exp01 ####
cl.exp1 <- read.table("/home1/minz11/knuh_kim/submission/rworks/in/clinic_exp01.tsv", sep = "\t", stringsAsFactors = F, comment.char = "", fill=T, quote = "", header = T, row.names = "SampleID")

str(cl.exp1)
cl.exp1$DSS <- as.character(cl.exp1$DSS)

cols <- c("0" = "#AFAFAF", "1" = "#91B0FF", "1.5" = "#0072B2")

### fig1B-1D ####
cl3_2_dapi.exp1 <- cl.exp1 %>%
  group_by(Age, DSS) %>%
  summarise(mean_value = mean(cl3_2_dapi, na.rm = TRUE),
            sd_value = sd(cl3_2_dapi, na.rm = TRUE),
            sem_value = sd(cl3_2_dapi, na.rm = T)/sqrt(n())) %>%
  ungroup()

fc.exp1 %>% 
  drop_na() %>% 
  ggplot(aes(x = Age, y = mean_value, fill = DSS, group = DSS)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = mean_value - sem_value,
                    ymax = mean_value + sem_value), width = 0.2, position = position_dodge(.9)) +
  geom_jitter(aes(x = Age, y = FC, group = DSS), data = subset.data.frame(cl.exp1, Age != "W5"), position = position_jitterdodge(jitter.width = .1, dodge.width = .9),size = 3) +
  scale_fill_manual(values = c("0" = "#AFAFAF", "1" = "#91B0FF", "1.5" = "#0072B2")) +
  labs(x = "Age", y = "Fecal calprotectin", title = "Exp01") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(0,14000)) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = .5), aspect.ratio = 1.25, axis.text = element_text(color = "black"))

# Wilcoxon rank-sum 테스트 수행
library(rstatix)

wilcox_results <- cl.exp1 %>%
  #filter(Age != "W5") %>%
  group_by(Age) %>%
  do(broom::tidy(pairwise.wilcox.test(.$infr, .$DSS, p.adjust.method = "bonferroni"))) %>%
  ungroup() %>%
  mutate(
    p.adj.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

print(wilcox_results)

# dunn_results <- cl.exp1 %>%
  # filter(Age != "W5") %>% 
  # group_by(Age) %>%
  # dunn_test(infr ~ DSS, p.adjust.method = "bonferroni") %>%
  # add_xy_position(x = "Age")

# p-value를 바탕으로 asterisk 할당
# dunn_results <- dunn_results %>%
  # mutate(asterisk = case_when(
    # p.adj < 0.001 ~ "***",
    # p.adj < 0.01 ~ "**",
    # p.adj < 0.05 ~ "*",
    # TRUE ~ "ns"
#  ))

# print(dunn_results)


### Exp02 ####
cl.exp2 <- read.table("/home1/minz11/knuh_kim/submission/rworks/in/clinic_exp02.tsv", sep = "\t", stringsAsFactors = F, comment.char = "", fill=T, quote = "", header = T, row.names = "SampleID")

str(cl.exp2)
cl.exp2$DSS <- as.character(cl.exp2$DSS)

### fig3B-3E ####
claudin3.exp2 <- cl.exp2 %>%
  group_by(DSS) %>%
  summarise(mean_value = mean(claudin3, na.rm = TRUE),
            sd_value = sd(claudin3, na.rm = TRUE),
            sem_value = sd(claudin3, na.rm = T)/sqrt(n())) %>%
  ungroup()

claudin3.exp2 %>% 
  drop_na() %>% 
  ggplot(aes(x = DSS, y = mean_value, fill = DSS, group = DSS)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = mean_value - sem_value,
                    ymax = mean_value + sem_value), width = 0.2, position = position_dodge(.9)) +
  geom_jitter(aes(x = DSS, y = claudin3, group = DSS), data = cl.exp2, position = position_jitterdodge(jitter.width = .1, dodge.width = .9),size = 3) +
  scale_fill_manual(values = c("0" = "#AFAFAF", "1" = "#91B0FF", "1.5" = "#0072B2")) +
  labs(x = "Group", y = "Claudin3", title = "Exp02") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.4)) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = .5), aspect.ratio = 1.25, axis.text = element_text(color = "black"))

# Wilcoxon rank-sum 테스트 수행
library(rstatix)

wilcox_results2 <- cl.exp2 %>%
  do(broom::tidy(pairwise.wilcox.test(.$claudin3, .$DSS, p.adjust.method = "bonferroni"))) %>%
  ungroup() %>%
  mutate(
    p.adj.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

print(wilcox_results2)
