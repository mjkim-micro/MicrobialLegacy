#### Exp03-04 #####
#### load packages and data #####
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
library(mixOmics)
library(gplots)

#### Exp03 ####
#### beta diversity ####
### PCoA ####
set.seed(1234)
plot_ordination(subset_samples(phy.03, Age == "W3"), ordinate(subset_samples(phy.03, Age == "W3"), "PCoA", "bray"), color = "Mother") +
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size = 1) +
  geom_hline(yintercept = 0, color = "grey", linetype="dashed", size = 1) +
  geom_point(aes(color = Mother), size = 3) +
  theme(aspect.ratio = 1) +
  theme(panel.background = element_rect(fill="white", color='black', linetype='solid', size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.title = element_text(size = 12, color = "black"), axis.text = element_text(size = 12, color = "black"), plot.title = element_text(hjust = .5)) +
  scale_color_manual(values = c("HC_HC" = "#AFAFAF", "UH-UH" = "#0072B2", "HC_UH" = "#FFB000", "UH_HC" = "#FE6100")) +
  stat_ellipse(aes(group = Mother)) +
  labs(x = "PCoA1  [6.5%]",y = "PCoA2  [5.6%]", title = "Exp 03")

# stats
set.seed(1234)
bc = phyloseq::distance(subset_samples(phy.03, Age == "W3"), method = "bray")
adonis2(bc ~ subset.data.frame(meta, Exp_No == "Exp03" & Age == "W3")$Mother)
pairwiseAdonis::pairwise.adonis(bc, subset.data.frame(meta, Exp_No == "Exp03" & Age == "W3")$Mother)

#### Exp04 ####
#### beta diversity ####
### PCoA ####
set.seed(1234)
plot_ordination(subset_samples(phy.03, Age %in% c("W10")), ordinate(subset_samples(phy.03, Age %in% c("W10")), "PCoA", "bray"), color = "Mother") +
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size = 1) +
  geom_hline(yintercept = 0, color = "grey", linetype="dashed", size = 1) +
  geom_point(aes(color = Mother, shape = Age), size = 3) +
  theme(aspect.ratio = 1) +
  theme(panel.background = element_rect(fill="white", color='black', linetype='solid', size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.title = element_text(size = 12, color = "black"), axis.text = element_text(size = 12, color = "black"), plot.title = element_text(hjust = .5)) +
  scale_color_manual(values = c("HC_HC" = "#AFAFAF", "UH_UH" = "#0072B2", "HC_UH" = "#FFB000", "UH_HC" = "#FE6100", "UH-UH" = "#0072B2")) +
  stat_ellipse(aes(group = Mother)) +
  labs(x = "PCoA1  [6.9%]",y = "PCoA2  [6.2%]", title = "Exp 04-W10")

# stats
set.seed(1234)
bc = phyloseq::distance(subset_samples(phy.03, Age %in% c("W10")), method = "bray")
adonis2(bc ~ subset.data.frame(meta, Exp_No == "Exp04" & Age == "W10")$Group3)
pairwiseAdonis::pairwise.adonis(bc, subset.data.frame(meta, Exp_No == "Exp04" & Age == "W10")$Group3)

### hierarchical clustering ####
phy.03@sam_data$Group2 <- paste0(phy.03@sam_data$Group, "_", phy.03@sam_data$DSS)
phy03.group <- merge_samples(subset_samples(phy.03, Age %in% c("W3")), "Group2")
phy03.group.rel <- transform_sample_counts(phy03.group, function(x) 100*x / sum(x))

phy03.group.otu <- as(otu_table(phy03.group.rel), "matrix")
phy03.group.dist <- vegdist(phy03.group.otu, method = "bray")

exp03.hclust <- hclust(phy03.group.dist, method = "ward.D2")
tree.03 <- as.phylo(exp03.hclust)

ggtree(tree.03) +
  geom_tiplab(size = 4, hjust = -0.05) +
  theme_tree2() +
  labs(title = "Hierarchical Clustering of Groups") +
  theme(plot.margin = margin(10, 100, 10, 10))

### distance comparison ####
bc.4 = phyloseq::distance(subset_samples(phy.03, Exp_No == "Exp04"), method = "bray")
bc.m.4 <- melt(as.matrix(bc.4))
bc.m.4[1:3,]

bc.m.4 <- bc.m.4 %>% filter(as.character(Var1) != as.character(Var2)) %>% 
  mutate_if(is.factor,as.character)

sd <- meta[,c("Group", "Mother", "DSS", "Age", "Group2")]
sd$Var1 <- rownames(sd)
bc.m2.4 <- left_join(bc.m.4, sd, by = "Var1")

colnames(sd) <- c("Group", "Mother", "DSS", "Age", "Group2", "Var2")
bc.m3.4 <- left_join(bc.m2.4, sd, by = "Var2")
colnames(bc.m3)

test <- bc.m3.4 %>% 
  subset(Mother.x == Mother.y & Age.x != Age.y)

test$concat <- paste0(test$Age.x, "_", test$Age.y)
test[1:3,]

# Calculate mean and standard deviation
dm.exp4.sum <- test %>%
  group_by(Mother.x, concat) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            sem_value = sd(value, na.rm = T)/sqrt(n())) %>%
  ungroup()

dm.exp4.sum <- dm.exp4.sum %>% 
  filter(concat %in% c("W8_W9", "W9_W10", "W8_W10"))

### stat
# Calculate Kruskal-Wallis for each Day and Parameter
test <- test %>% 
  filter(concat %in% c("W8_W9", "W9_W10", "W8_W10"))

set.seed(1234)
dunn.test(subset.data.frame(test, concat %in% c("W9_W10"))$value, subset.data.frame(test, concat %in% c("W9_W10"))$Mother.x, method="bonferroni")

test %>%
  filter(concat %in% c("W8_W9", "W9_W10", "W8_W10")) %>%
  ggplot(aes(fill = factor(Mother.x, levels = c("HC_HC", "UH_UH", "HC_UH", "UH_HC")), y = value, x = factor(concat, levels = c("W8_W9", "W9_W10", "W8_W10")))) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("HC_HC" = "#AFAFAF", "UH_UH" = "#0072B2", "HC_UH" = "#FFB000", "UH_HC" = "#FE6100")) +
  theme_bw() +
  theme(aspect.ratio = .75) +
  labs(x = "Comparison", y = "Bray-Curtis distance", fill = "Group")  

### dimdi exp04 ####
dimdi.out %>% 
  filter(Exp_No %in% c("Exp04")) %>% 
  filter(Group %in% c("Con", "Cross")) %>% 
  filter(DSS != "1") %>% 
  ggplot(aes(x = Age, y = dimdi, fill = Group2)) +
  geom_boxplot() +
  #geom_jitter(aes(x = Age, y = dimdi, group = DSS), position = position_jitterdodge(jitter.width = .1, dodge.width = .9),size = 2) +
  #scale_fill_manual(values = c("0" = "#AFAFAF", "1" = "#91B0FF", "1.5" = "#0072B2")) +
  #scale_fill_manual(values = c("Con_0" = "#AFAFAF", "Con_1.5" = "#0072B2", "Cross_0" = "#FFB000", "Cross_1.5" = "#FE6100")) +
  scale_fill_manual(values = c("Con_0" = "#AFAFAF", "Con_1.5" = "#0072B2", "Cross_0" = "#FFB000", "Cross_1.5" = "#FE6100")) +
  labs(x = "Age", y = "DiMDI") +
  theme_classic() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = .5), aspect.ratio = 1) +
  ylim(c(0,3.8))

set.seed(1234)
dimdi.wilcox <- dimdi.out %>%
  filter(Exp_No %in% c("Exp04")) %>%
  filter(Group %in% c("Con", "Cross")) %>% 
  group_by(Age) %>%
  do(broom::tidy(pairwise.wilcox.test(.$dimdi, .$Group2, p.adjust.method = "fdr"))) %>%
  ungroup() %>%
  mutate(
    p.adj.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

figs5.dimdi <-  dimdi.out %>%
  filter(Exp_No %in% c("Exp04")) %>%
  filter(Group %in% c("Con", "Cross"))
  
kruskal.test(dimdi ~ Group3, data = subset.data.frame(figs5.dimdi, Age == "W10"))
dunn.test(subset.data.frame(figs5.dimdi, Age == "W10")$dimdi, subset.data.frame(figs5.dimdi, Age == "W10")$Group3, method="bonferroni")

### hierarchical clustering ####
phy.04@sam_data$Group2 <- paste0(phy.03@sam_data$Group, "_", phy.03@sam_data$DSS)
phy04.group <- merge_samples(subset_samples(phy.03, Age %in% c("W10")), "Group2")
phy04.group.rel <- transform_sample_counts(phy04.group, function(x) 100*x / sum(x))

phy04.group.otu <- as(otu_table(phy04.group.rel), "matrix")
phy04.group.dist <- vegdist(phy04.group.otu, method = "bray")

exp04.hclust <- hclust(phy04.group.dist, method = "ward.D2")
tree.04 <- as.phylo(exp04.hclust)

ggtree(tree.04) +
  geom_tiplab(size = 4, hjust = -0.05) +
  theme_tree2() +
  labs(title = "Hierarchical Clustering of Groups") +
  theme(plot.margin = margin(10, 100, 10, 10))
