### revision #####
#### stacked barplot at the genus lv #####
genus <- as.data.frame(otu_table(aggregate_rare(microbiome::transform(phy, "compositional"), level = "Genus", detection = 1/100, prevalence = 15/100)))
genus <- as.data.frame(t(genus))
genus <- merge(genus, meta, by=0)
row.names(genus) <- genus$Row.names
genus <- genus[,-1]

colnames(genus)
genus.df <- genus %>% 
  pivot_longer(`[Eubacterium]_xylanophilum_group`:Unknown, names_to = "Genus", values_to = "Value")

ggplot(subset.data.frame(genus.df, Exp_No == "Exp06"), aes(x=SampleID, y=Value, fill=reorder(Genus, Value), group=reorder(Genus, Value))) + 
  geom_bar(stat = "identity") +
  facet_wrap(~Group3+Age, scales = "free", ncol = 6) +
  theme(panel.background = element_rect(fill="white", color='black', linetype='solid', size=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 10, color = "black"), axis.ticks = element_line(color = "black", size = 1), legend.position = "right") +
  scale_fill_manual(values = c("#E25B6D", "#FBC02D", "#8E24AA", "#F57C00", "#7B6FA2", "#13B879", "#7D6608", "#0097A7", "#00796B", "#689F38", "#8a2be2", "#ffebcd", "#5f8ea0", "#e25ac9", "#7fbdd3", "#a0c8ff", "#f55f20", "#b0744a", "#983f4a", "#ff9999", "#f88379", "#c154c1")) +
  theme(aspect.ratio = 2) #genus

# mean abundance ####
genus.avg <- genus.df %>%
  filter(Exp_No == "Exp06") %>%
  group_by(Group3, Age, Genus) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE),
            SE = sd(Value, na.rm = TRUE) / sqrt(n())) %>%
  ungroup()

# plot
ggplot(genus.avg, aes(x=Age, y=Mean_Value, fill=reorder(Genus, Mean_Value))) + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Group3, scales = "free_x", ncol = 3) +  
  theme(panel.background = element_rect(fill="white", color='black', 
                                        linetype='solid', size=1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 15, color = "black"), 
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(color = "black", size = 1), 
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#E25B6D", "#FBC02D", "#8E24AA", "#F57C00", 
                               "#7B6FA2", "#13B879", "#7D6608", "#0097A7", 
                               "#00796B", "#689F38", "#8a2be2", "#ffebcd", 
                               "#5f8ea0", "#e25ac9", "#7fbdd3", "#a0c8ff", 
                               "#f55f20", "#b0744a", "#983f4a", "#ff9999", 
                               "#f88379", "#c154c1")) +
  labs(x = "Age", y = "Relative Abundance", fill = "Genus") +
  theme(aspect.ratio = 1)

### MaAsLin2 ####
library(Maaslin2)

phy.genus <- tax_glom(phy, taxrank = "Genus", NArm = TRUE)
phy.genus.filt <- filter_taxa(phy.genus, 
                                  function(x) sum(x > 0) >= (0.1 * length(x)), 
                                  TRUE)

# 3. Feature table 준비 (taxa in rows, samples in columns)
g.maaslin <- as.data.frame(otu_table(subset_samples(phy.genus.filt, Group == "Con" & Age == "W10")))
g.tax <- as.data.frame(tax_table(phy.genus.filt))

if ("Genus" %in% colnames(g.tax)) {
  genus_names <- ifelse(is.na(g.tax$Genus) | g.tax$Genus == "", 
                        paste0("Unknown_", g.tax$Family),
                        g.tax$Genus)
  
  genus_names <- make.unique(genus_names, sep = "_")
  
  rownames(g.maaslin) <- genus_names
}

fit_genus <- Maaslin2(
  input_data = g.maaslin,
  input_metadata = subset.data.frame(meta, Group == "Con" & Age == "W10"),
  output = "suscp_w10_maaslin",
  fixed_effects = c("Group3"),
  reference = c("Group3,Con_0"),
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  min_abundance = 0.0001,
  min_prevalence = 0.1,
  max_significance = 0.25,
  cores = 4
)

### volcano plot between W8-W9 to response R1-C6 ####
maas.c6 <- read.delim("suscp_2p_maaslin/all_results.tsv", stringsAsFactors = FALSE)
# 3. Significance threshold
maas.c6 <- maas.c6 %>%
  mutate(
    Significant = case_when(
      qval < 0.05 ~ "q < 0.05",
      qval < 0.25 ~ "q < 0.25",
      TRUE ~ "Not Significant"
    ),
    Significant = factor(Significant, levels = c("q < 0.05", "q < 0.25", "Not Significant"))
  )

# 4. Top significant genera
top_genera <- maas.c6 %>%
  filter(qval < 0.25) %>%
  arrange(qval) %>%
  head(15) %>%
  pull(feature)

maas.c6 <- maas.c6 %>%
  mutate(Label = ifelse(feature %in% top_genera, feature, ""))

# 5. Volcano plot
ggplot(maas.c6, aes(x = coef, y = -log10(qval))) +
  geom_point(aes(color = Significant), alpha = 0.6, size = 2.5) +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed", color = "grey40", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey20", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", size = 0.5) +
  geom_text_repel(aes(label = Label),
                  size = 3.5,
                  max.overlaps = 20,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.3,
                  fontface = "italic",
                  min.segment.length = 0) +
  scale_color_manual(values = c("q < 0.05" = "#E25B6D",
                                "q < 0.25" = "#FBC02D", 
                                "Not Significant" = "grey70")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.position = "right",
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        aspect.ratio = 1) +
  labs(x = "Coefficient (Effect Size)",
       y = "-log10(q-value)",
       color = "Significance",
       title = "Differential genera - pre vs post-DSS 0%") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_continuous(expand = expansion(mult = 0.1))

### heatmap plot between offspring groups(Exp01) to response R1-C2 ####
libary(reshape2)

# 1. maaslin2 results for each week
week_folders <- c("suscp_w8_maaslin", "suscp_w9_maaslin", "suscp_w10_maaslin")
week_labels <- c("W8", "W9", "W10")

all_results <- data.frame()

for (i in 1:length(week_folders)) {
  results <- read.delim(paste0(week_folders[i], "/significant_results.tsv"), 
                        stringsAsFactors = FALSE)
  
  results$Week <- week_labels[i]
  all_results <- rbind(all_results, results)
}


# 1. significant genera
sig_genera <- all_results %>%
  filter(qval < 0.25) %>%
  pull(feature) %>%
  unique()

# 2. Coefficient 
coef_data <- all_results %>%
  filter(feature %in% sig_genera) %>%
  select(feature, Week, coef)

# 3. dcast 
coef_matrix <- dcast(coef_data, feature ~ Week, value.var = "coef", fill = 0)

# 4. Rownames set
rownames(coef_matrix) <- coef_matrix$feature
coef_matrix$feature <- NULL

# 5. Matrix
coef_matrix <- as.matrix(coef_matrix)

# 확인
print("Coefficient matrix:")
print(coef_matrix)
print(class(coef_matrix))
print(str(coef_matrix))

# 6. Significance matrix
sig_matrix <- matrix("", nrow = nrow(coef_matrix), ncol = ncol(coef_matrix))
rownames(sig_matrix) <- rownames(coef_matrix)
colnames(sig_matrix) <- colnames(coef_matrix)

for (i in 1:nrow(sig_info)) {
  genus <- sig_info$feature[i]
  week <- sig_info$Week[i]
  sig_mark <- sig_info$sig[i]
  
  if (genus %in% rownames(sig_matrix) && week %in% colnames(sig_matrix)) {
    sig_matrix[genus, week] <- sig_mark
  }
}

# 7. Heatmap
library(pheatmap)

pheatmap(coef_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-max(abs(coef_matrix)), max(abs(coef_matrix)), 
                      length.out = 101),
         display_numbers = sig_matrix,
         number_color = "black",
         fontsize_number = 14,
         fontsize_row = 11,
         fontsize_col = 12,
         angle_col = 0,
         border_color = "grey60",
         cellwidth = 50,
         cellheight = 30,
         main = "Differential Abundance Across Weeks (Exp01)",
         width = 7,
         height = 7)

# 1. maaslin2 results
week_folders <- c("suscp_w8_maaslin", "suscp_w9_maaslin", "suscp_w10_maaslin")
week_labels <- c("W8", "W9", "W10")

all_results <- data.frame()

for (i in 1:length(week_folders)) {
  results <- read.delim(paste0(week_folders[i], "/significant_results.tsv"), 
                        stringsAsFactors = FALSE)
  results$Week <- week_labels[i]
  all_results <- rbind(all_results, results)
}

# 1. Phyloseq data
phy_genus <- tax_glom(phy, taxrank = "Genus", NArm = TRUE)
phy_comp <- microbiome::transform(phy_genus, "compositional")

# 2. OTU table
otu_data <- as.data.frame(otu_table(phy_comp))
if (!taxa_are_rows(phy_comp)) {
  otu_data <- t(otu_data)
}

tax_data <- as.data.frame(tax_table(phy_comp))
genus_names <- tax_data$Genus
genus_names[is.na(genus_names)] <- "Unknown"
genus_names <- make.unique(genus_names, sep = "_")
rownames(otu_data) <- genus_names

# 3. significant genera
sig_genera <- all_results %>%
  pull(feature) %>%
  unique()

sig_genera <- c("Turicibacter", "Bifidobacterium", "Rikenellaceae_RC9_gut_group", "Candidatus_Arthromitus", "Helicobacter", "Tannerellaceae", "Bacteroides", "Gastranaerophilales", "Alistipes", "Coriobacteriaceae_UCG-002", "Ruminococcaceae", "uncultured_1", "[Eubacterium]_siraeum_group", "[Eubacterium]_coprostanoligenes_group", "Oscillibacter", "Faecalibaculum", "A2", "uncultured", "[Eubacterium]_xylanophilum_group", "Muribaculaceae", "Romboutsia", "Clostridia_UCG-014", "Desulfovibrio")

#sig_genera <- c("Bifidobacterium", "Escherichia-Shigella", "Romboutsia", "Streptococcus", "Ligilactobacillus", "Staphylococcus", "Staphylococcaceae", "Desulfovibrio", "uncultured_1", "uncultured_3", "Mammalia", "[Eubacterium]_ventriosum_group", "Faecalibaculum", "Clostridia_UCG-014", "Ruminococcaceae", "[Eubacterium]_coprostanoligenes_group", "Enterorhabdus", "Clostridium_sensu_stricto_1", "Akkermansia", "uncultured_2", "Family_XIII_UCG-001", "Clostridia_vadinBB60_group", "Odoribacter", "[Eubacterium]_xylanophilum_group", "Candidatus_Arthromitus", "Incertae_Sedis", "Parasutterella", "Bilophila", "Muribaculaceae", "Colidextribacter", "Roseburia", "Lactobacillus", "Mucispirillum", "HT002")
otu_sig <- otu_data[rownames(otu_data) %in% sig_genera, , drop = FALSE]

# 4. W1, W3, W5
samples_w <- meta %>%
  filter(Age %in% c("W8", "W9", "W10")) %>%
  filter(Group == "Con") %>% 
  rownames()

otu_sig_w <- otu_sig[, samples_w, drop = FALSE]

# 5. calculate mean for each group
meta_w <- meta[samples_w, ]

# Group + Age
otu_sig_w <- as.data.frame(t(otu_sig_w))
otu.sig.df <- merge(otu_sig_w, meta_w, by = 0)
rownames(otu.sig.df) <- otu.sig.df$Row.names
otu.sig.df <- otu.sig.df[,-1]
colnames(otu.sig.df)
otu.sig.ddf <- otu.sig.df %>% 
  #pivot_longer(names_to = "Genus", values_to = "Relabs", Mammalia:Roseburia)
  pivot_longer(names_to = "Genus", values_to = "Relabs", uncultured:A2)

maas.c1 <- otu.sig.ddf %>%
  group_by(Age, Genus, DSS) %>%
  summarise(mean_value = mean(Relabs, na.rm = TRUE),
            sd_value = sd(Relabs, na.rm = TRUE),
            sem_value = sd(Relabs, na.rm = T)/sqrt(n())) %>%
  ungroup()

maas.c3 <- otu.sig.ddf %>%
  group_by(Age, Genus, DSS) %>%
  summarise(mean_value = mean(Relabs, na.rm = TRUE),
            sd_value = sd(Relabs, na.rm = TRUE),
            sem_value = sd(Relabs, na.rm = T)/sqrt(n())) %>%
  ungroup()

# 6. Wide format
maas.c1 <- maas.c1 %>%
  mutate(Group3 = paste(Age, DSS, sep = "_"))

maas.c1.df <- maas.c1 %>%
  select(Genus, Group3, mean_value) %>%
  pivot_wider(names_from = Group3, values_from = mean_value) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

maas.c3 <- maas.c3 %>%
  mutate(Group3 = paste(Age, DSS, sep = "_"))

maas.c3.df <- maas.c3 %>%
  select(Genus, Group3, mean_value) %>%
  pivot_wider(names_from = Group3, values_from = mean_value) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

# 7. Column annotation
col_anno <- data.frame(
  Week = gsub("_.*", "", colnames(abundance_wide)),
  Group = gsub(".*_", "", colnames(abundance_wide)),
  row.names = colnames(abundance_wide)
)

# 8. Abundance Heatmap
hm.cols <- colorRampPalette(c("#3869a8", "#f4f4f4", "#fb4545"))

heatmap.2(as.matrix(maas.c3.df), scale = "row", col = hm.cols, trace = "none", density.info = "none", Colv = "T",Rowv = "T", dendrogram = "none", notecex = 1, notecol = "black", cexRow = .8)

### coefficient barplot ####
# 2. significant genera
sig.genera.coef <- all_results %>%
  filter(qval < 0.25)

### re-analyzing to respond R1-C6 ####
# 1. Phyloseq data
phy_genus <- tax_glom(phy, taxrank = "Genus", NArm = TRUE)
phy_comp <- microbiome::transform(phy_genus, "compositional")

# 2. OTU table
otu_data <- as.data.frame(otu_table(phy_comp))
if (!taxa_are_rows(phy_comp)) {
  otu_data <- t(otu_data)
}

tax_data <- as.data.frame(tax_table(phy_comp))
genus_names <- tax_data$Genus
genus_names[is.na(genus_names)] <- "Unknown"
genus_names <- make.unique(genus_names, sep = "_")
rownames(otu_data) <- genus_names

# 3. significant genera
sig_genera.2.df <- as.data.frame(read_tsv("exp02_dss_maaslin/significant_results.tsv"))
rownames(sig_genera.2.df) <- sig_genera.2.df$feature

sig_genera.2 <- c("Lactobacillus", "Gastranaerophilales", "Parasutterella", "[Eubacterium]_ventriosum_group", "Ligilactobacillus", "[Eubacterium]_coprostanoligenes_group", "Ruminococcaceae", "uncultured_1", "Faecalibaculum", "Alistipes", "UCG-005", "[Eubacterium]_xylanophilum_group", "Roseburia", "Mucispirillum", "Ruminococcus")

otu_sig <- otu_data[rownames(otu_data) %in% sig_genera.2, , drop = FALSE]

# 4. W1, W3, W5
samples_w <- meta %>%
  filter(Age %in% c("W8", "W9", "W10")) %>%
  filter(Exp_No == "Exp02") %>% 
  rownames()

otu_sig_w <- otu_sig[, samples_w, drop = FALSE]

# 5. calculate mean for each group
meta_w <- meta[samples_w, ]

# Group + Age별 평균
otu_sig_w <- as.data.frame(t(otu_sig_w))
otu.sig.df <- merge(otu_sig_w, meta_w, by = 0)
rownames(otu.sig.df) <- otu.sig.df$Row.names
otu.sig.df <- otu.sig.df[,-1]
colnames(otu.sig.df)
otu.sig.ddf <- otu.sig.df %>% 
  #pivot_longer(names_to = "Genus", values_to = "Relabs", Mammalia:Roseburia)
  pivot_longer(names_to = "Genus", values_to = "Relabs", uncultured_1:Roseburia)

maas.c6 <- otu.sig.ddf %>%
  group_by(Age, Genus, DSS) %>%
  summarise(mean_value = mean(Relabs, na.rm = TRUE),
            sd_value = sd(Relabs, na.rm = TRUE),
            sem_value = sd(Relabs, na.rm = T)/sqrt(n())) %>%
  ungroup()

# 6. Wide format
maas.c6 <- maas.c6 %>%
  mutate(Group3 = paste(Age, DSS, sep = "_"))

maas.c6.df <- maas.c6 %>%
  select(Genus, Group3, mean_value) %>%
  pivot_wider(names_from = Group3, values_from = mean_value) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

# 7. Column annotation
col_anno <- data.frame(
  Week = gsub("_.*", "", colnames(abundance_wide)),
  Group = gsub(".*_", "", colnames(abundance_wide)),
  row.names = colnames(abundance_wide)
)

# 8. Abundance Heatmap
hm.cols <- colorRampPalette(c("#3869a8", "#f4f4f4", "#fb4545"))

heatmap.2(as.matrix(maas.c6.df), scale = "row", col = hm.cols, trace = "none", density.info = "none", Colv = "T",Rowv = "T", dendrogram = "none", notecex = 1, notecol = "black", cexRow = .8)


### quantification of Fig.4i ####
# Distance
bc.fig4 <- phyloseq::distance(subset_samples(phy, Exp_No == "Exp05"), method = "bray")
bc.fig4.m <- melt(as.matrix(bc.fig4))
bc.fig4.m[1:3,]
meta.fig4 <- data.frame(sample_data(subset_samples(phy, Exp_No == "Exp05")))
