setwd("~/Collaborations_DKFZ/Amagdela project")

library(esquisse)
library(rlang)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(tidyverse)
library(BiocManager)
library(clusterProfiler)
library(GGally)
library(ggrepel)
library(naniar)
library(limma)
library(pheatmap)
library(readr)
library(DEqMS)
library(enrichplot)
library(msigdbr)
library(ggpubr)
library(stringr)
library(tibble)
library(ggrepel)
library(EnhancedVolcano)
library(fgsea)
library(openxlsx)
library(RColorBrewer)
library(tidyr)
library(broom)
library(enrichR)
library(fmsb)
library(corrplot)
library(factoextra)
library(corrplot)
library(reshape2)
library(conflicted)

Ama_Proteo <- read.csv ("~/Collaborations_DKFZ/Amagdela project/Amygdala_Final20Jan2024-Azmal.csv", header = TRUE, na.strings= TRUE)

colnames(Ama_Proteo)
colnames(Ama_Proteo) <- str_replace_all(colnames(Ama_Proteo), "\\.",replacement = "_")


Ama_Proteo <- Ama_Proteo %>% dplyr::mutate(Protein_Group = str_split(Protein_Group,";",simplify = TRUE)[,1])

Ama_Proteo <- Ama_Proteo %>%
  rename_with(~ str_replace(., "X20231213_FS1_DIA_Amag_Control_Rep", "Amag_Control_Rep"), starts_with("X20231213_FS1_DIA_Amag_Control_Rep"))

Ama_Proteo <- Ama_Proteo %>%
  rename_with(~ str_replace(., "X20231213_FS1_DIA_Amag_Fear_Rep", "Amag_Fear_Rep"), starts_with("X20231213_FS1_DIA_Amag_Fear_Rep"))

Ama_Proteo_Fear <- Ama_Proteo %>% dplyr::select (Protein_Group, Genes, matches("Fear_Rep"))
Ama_Proteo_Fear [3:7]<- log2(Ama_Proteo_Fear[3:7])
Ama_Proteo_con <- Ama_Proteo %>% dplyr::select (Protein_Group, Genes, matches("Amag_Control"))
Ama_Proteo_con [3:7] <- log2(Ama_Proteo_con[3:7])
FCAma_Proteo <- Ama_Proteo_Fear[3:7] - Ama_Proteo_con[3:7]
FCAma_Proteo$Protein_Group <- Ama_Proteo$Protein_Group
FCAma_Proteo$Genes <- Ama_Proteo$Genes
##-------------------------Differential exprssion----------------------------------------------------------------------
# preselect for  valid values # take care to adjust the is.na statement to really match your specific columns
FCAma_Proteo_main <-  FCAma_Proteo %>% 
  dplyr::select(Protein_Group, Genes, matches("Fear")) %>%
  filter(rowSums(is.na(FCAma_Proteo[3:7])) < 3)

FCAma_Proteo_main1 <- FCAma_Proteo_main[-1:-2]
rownames(FCAma_Proteo_main1) <- FCAma_Proteo_main$Genes


#limma analysis
fit1 <- lmFit(FCAma_Proteo_main1)# Select column replicates for Limma analysis
fit2 <- eBayes(fit1, robust = TRUE) # Perform Bayes statistics

#output proteinGroups_limma and DEqMS
table <- topTable(fit = fit2, number = Inf, genelist = FCAma_Proteo_main[,c("Genes")], coef = 1, adjust.method = "BH", sort.by = "logFC",resort.by = "logFC", p.value = 1,lfc = 0) # Limma output table
table_FCAma_Proteo <- rownames_to_column(table, var="Genes")

##extract top 10 up/down regulated proteins
top10_up1 <- table_FCAma_Proteo %>%
  filter(P.Value< 0.05) %>%
  arrange(desc(logFC))%>%
  dplyr::slice(1:10) 

top10_down1 <- table_FCAma_Proteo %>%
  filter(P.Value< 0.05) %>%
  arrange(logFC) %>%
  dplyr::slice(1:20)

EnhancedVolcano(table,
                lab = rownames(table),
                x = 'logFC',
                y = 'P.Value',
                title = 'Amygdala Stress',
                pCutoff = 10e-3,
                FCcutoff = 1,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)


##--------------------------------------------------------------------------------------
sortedTable <- table_FCAma_Proteo[order(table_FCAma_Proteo$logFC), ]
sortedTable <- rename(sortedTable, GeneID = ID)
ranks <- setNames(object = sortedTable$logFC, nm = sortedTable$GeneID)

as.data.frame(ranks)
as.matrix(ranks)

plot(ranks)
geneSets <- gmtPathways("~/Collaborations_DKFZ/Amagdela project/GSEAdb/m8.all.v2023.2.Mm.symbols.gmt")


CanonicalPWgeneSets <- gmtPathways("~/Collaborations_DKFZ/Amagdela project/GSEAdb/Canonicalpathways.v2023.2.Mm.symbols.gmt")
View(CanonicalPWgeneSets)


# Perform fgsea using fgseaMultilevel
fgseaRes_CanonicalPWgeneSets <- fgsea(pathways = CanonicalPWgeneSets, 
                  stats = ranks,
                  minSize = 10, 
                  maxSize = 1000)

# Above one is best parameters
fgsea_CanonicalPWgeneSets <- fgsea(pathways = CanonicalPWgeneSets, 
                   stats = ranks,
                   eps      = 0.0,
                   minSize  = 15,
                   maxSize  = 500)
# Save dataframe as an Excel file
write.xlsx(fgseaRes, "fgsea_CanonicalPWgeneSets.xlsx")
# -------------------------------------------------------------------

CCsubsetofGOgeneSets <- gmtPathways("~/Collaborations_DKFZ/Amagdela project/GSEAdb/CCsubsetofGO.v2023.2.Mm.symbols.gmt")

# Perform fgsea using fgseaMultilevel
fgsea_CCsubsetofGOgeneSets <- fgsea(pathways = CCsubsetofGOgeneSets, 
                                   stats = ranks,
                                   eps      = 0.0,
                                   minSize  = 15,
                                   maxSize  = 500)


# Save dataframe as an Excel file
write.xlsx(fgseaRes, "fgsea_CCsubsetofGOgeneSets.xlsx")
# -------------------------------------------------------------------

ontologygenesets <- gmtPathways("~/Collaborations_DKFZ/Amagdela project/GSEAdb/ontologygenesets.v2023.2.Mm.symbols.gmt")

# Perform fgsea using fgseaMultilevel
fgsea_ontologygenesets <- fgsea(pathways = ontologygenesets, 
                                    stats = ranks,
                                    eps      = 0.0,
                                    minSize  = 15,
                                    maxSize  = 500)

# Save dataframe as an Excel file
write.xlsx(fgseaRes, "fgsea_ontologygenesets.xlsx")
# -------------------------------------------------------------------
regulatorytargetgenesets <- gmtPathways("~/Collaborations_DKFZ/Amagdela project/GSEAdb/regulatorytargetgenesets.v2023.2.Mm.symbols.gmt")

# Perform fgsea using fgseaMultilevel
fgsea_regulatorytargetgenesets <- fgsea(pathways = regulatorytargetgenesets, 
                                stats = ranks,
                                eps      = 0.0,
                                minSize  = 15,
                                maxSize  = 500)


# Save dataframe as an Excel file
write.xlsx(fgseaRes, "fgsea_regulatorytargetgenesets.xlsx")


# ------------------------search the specific gene in the list -------------------------------------------

# Searching of specific in the list is possible for identifiying the pattern
# using the the grepl function, as below

Kcnj <- Ama_Proteo %>%
  filter(grepl("kcnj", Genes, ignore.case = TRUE))

ras <- Ama_Proteo %>% 
  filter(grepl("ras", Genes, ignore.case = T))


GPC1 <-  Ama_Proteo %>% 
  filter(grepl("protein coupled", First_Protein_Description, ignore.case = T))

df_GPCR <- GPC1
df <- df_GPCR %>% select("Genes", 5:14)


write.xlsx(df, "df.xlsx")

df <- Kcnj
colnames(df)

df <- df %>% select("Genes", 5:14)

conflict_scout()

df_long <- df %>%
  tidyr::pivot_longer(
    cols = -Genes,
    names_to = c("condition", "replicate"),
    names_pattern = "^(.*)_Rep(\\d+)$"
  ) %>% mutate(value = log2(value))

# df_long <- df_long %>% 
#   mutate(value = log2(value))
# 
# 
# df_long <- df %>%
#   pivot_longer(
#     cols = -Genes,
#     names_to = c("condition", "replicate"),
#     names_pattern = "^(.*)_Rep(\\d+)$"
#   ) %>%
#   mutate(value = log2(value))


# Enhanced ggplot code for better visualization
ggplot(df_long, aes(x = condition, y = value, fill = condition)) +
  geom_violin(trim = FALSE, scale = "area", color = "gray40", alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 1, color = "red", shape = 18) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +
  facet_wrap(~ Genes, scales = "free_y", strip.position = "top") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    strip.text.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Condition",
    y = "Log2 Expression",
    fill = "Condition",
    title = "Gene Expression by Condition"
  ) +
  scale_fill_brewer(palette = "Set1")  # Change the color palette

# This is for box plot
ggplot(df_long, aes(x = condition, y = value, fill = condition)) +
  geom_boxplot(trim = FALSE, scale = "area", color = "gray40", alpha = 0.4) +
  stat_summary(fun = mean, geom = "point", size = 1, color = "red", shape = 18) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +
  facet_wrap(~ Genes, scales = "free_y", strip.position = "top") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    strip.text.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Condition",
    y = "Log2 Expression",
    fill = "Condition",
    title = "Gene Expression of G-protein coupled receptor"
  ) +
  scale_fill_brewer(palette = "Set2")  # Change the color palette



# Calculate mean values for each combination of Genes and condition
statistics <- df_long %>%
  group_by(Genes, condition) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')

# Create the plot
ggplot(df_long, aes(x = condition, y = value, fill = condition)) +
  geom_violin(trim = FALSE, scale = "area", color = "gray40", alpha = 0.4) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +
  geom_text(data = statistics, aes(label = sprintf("Mean: %.2f", mean_value), y = mean_value + 0.1), 
            position = position_dodge(width = 0.75), vjust = 0, color = "blue") +
  facet_wrap(~ Genes, scales = "free_y", strip.position = "top") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    strip.text.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Condition",
    y = "Log2 Expression",
    fill = "Condition",
    title = "Gene Expression by Condition"
  ) +
  scale_fill_brewer(palette = "Set1")

# Perform t-tests for each gene
results <- df_long %>%
  group_by(Genes) %>%
  do(tidy(t.test(value ~ condition, data = .))) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))  # Adjusting p-values using Benjamini-Hochberg

# Filter for significant results to display (filter will create the problem for ploting, uneven)
# significant_results <- results %>%
#    dplyr::filter(p.adj < 0.05) %>%
#   mutate(signif_label = if_else(p.adj < 0.001, "***", if_else(p.adj < 0.01, "**", "*")))

# Apply labels to all results, including non-significant ones
GPCR_Sig_table_significant_results <- results %>%
  mutate(signif_label = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE          ~ "ns"  # Label for non-significant results
  )) 

write.xlsx(GPCR_Sig_table_significant_results, "GPCR_Sig_table_significant_results.xlsx")


# Create the plot
ggplot(df_long, aes(x = condition, y = value, fill = condition)) +
  geom_violin(trim = FALSE, scale = "area", color = "gray40", alpha = 0.4) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +
  # geom_text(data = significant_results, aes(label = signif_label, y = Inf), vjust = -0.5, color = "red") +
  facet_wrap(~ Genes, scales = "free_y", strip.position = "top") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    strip.text.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Condition",
    y = "Log2 Expression",
    fill = "Condition",
    title = "Gene Expression by Condition"
  ) +
  scale_fill_brewer(palette = "Set1")




# Create the plot

ggplot(df_long, aes(x = condition, y = value, fill = condition)) +
  geom_violin(trim = FALSE, scale = "area", color = "gray40", alpha = 0.4) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +
  geom_text(data = significant_results, aes(label = signif_label, x = Inf, y = Inf), 
            hjust = 1.1, vjust = -0.5, color = "red", inherit.aes = FALSE) +
  facet_wrap(~ Genes, scales = "free_y", strip.position = "top") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    strip.text.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Condition",
    y = "Log2 Expression",
    fill = "Condition",
    title = "Gene Expression by Condition"
  ) +
  scale_fill_brewer(palette = "Set1")



# Calculate a position for the text based on the data
max_values <- df_long %>%
  group_by(Genes) %>%
  summarise(max_value = max(value, na.rm = TRUE))

# # Merge this with your significant_results
# significant_results <- significant_results %>%
#   left_join(max_values, by = "Genes") %>%
#   mutate(y_pos = max_value + (0.1 * max_value))  # Adjust position as needed

# Merge this with your significant_results
significant_results <- significant_results %>%
  left_join(max_values, by = "Genes") %>%
  mutate(y_pos = max_value )  # Adjust position as needed + (0.1 * max_value)


# Update plot code
 ggplot(df_long, aes(x = condition, y = value, fill = condition)) +
  geom_violin(trim = FALSE, scale = "area", color = "gray40", alpha = 0.4) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +
  # geom_text(data = significant_results, aes(x = 1, y = y_pos, label = signif_label), 
  #           hjust = -5, vjust = 0, color = "black", size = 6, inherit.aes = FALSE) +
  facet_wrap(~ Genes, scales = "free_y", strip.position = "top") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    strip.text.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Condition",
    y = "Log2 Expression",
    fill = "Condition",
    title = "Kcnj Family Protein Expression by Condition"
  ) +
  scale_fill_brewer(palette = "Set1")

cs
ggsave("cs.tiff", plot = cs, height = 20, width = 30)

write.xlsx(significant_results, "significant_results.xlsx")

dim(significant_results)

# -------------------------------------------------------------------------

# Define the term you are looking for
term <- "HEDGEHOG"

# Filter the DataFrame where the 'pathway' column contains the term 'synaptic'
filtered_dfGO <- fgseaRes_CanonicalPWgeneSets %>%
  filter(grepl(term, pathway, ignore.case = TRUE))  # Case-insensitive partial match search in 'pathway'

ggplot(filtered_dfGO, aes(x = reorder(pathway, NES), y = NES)) +
  geom_point(aes(size = size, color = padj), alpha = 0.6) + # padj as an example, replace with any relevant metric
  scale_size_continuous(range = c(5, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  xlab("Pathway") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top 10 Pathways") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(filtered_dfGO, aes(x=reorder(pathway, NES), y=NES)) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Pathway") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top 30 Pathways")

# Filter to keep only rows where NES is positive
filtered_dfGO <- filtered_dfGO %>% 
  filter(NES > 0)

# Now create the plot
ggplot(filtered_dfGO, aes(x=reorder(pathway, NES), y=NES)) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Pathway") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top Pathways with Positive NES")


ggplot(filtered_dfGO, aes(x = reorder(pathway, NES), y = NES, size = size, color = pval)) +
  geom_point(alpha = 0.9) +  # Use points, with some transparency
  scale_color_gradient(low = "black", high = "blue") +  # Color gradient from high to low P-values
  scale_size(range = c(4, 20), name = "Size") +
  coord_flip() +  # Flip coordinates for better label reading
  labs(x = "Pathway", 
       y = "Normalized Enrichment Score (NES)", 
       size = "Size", 
       color = "P-value") +
  ggtitle("Visualization of Pathways") +
  theme_minimal() +  # Minimalist theme
  theme(legend.position = "right")  # Adjust legend position


colnames(fgsea_CP)

fgsea_CP <- fgsea_CanonicalPWgeneSets %>% 
  filter(pval < 0.05, NES >0)

fgsea_CP_N <- fgsea_CanonicalPWgeneSets %>% 
  filter(pval < 0.05, NES <0)

fgsea_GO_Po <- fgsea_ontologygenesets %>% 
  dplyr::filter(padj < 0.05, NES >0)

fgsea_GO_Neg <- fgsea_ontologygenesets %>% 
  dplyr::filter(pval < 0.05, NES <0)

UB_GSEA <- rbind(fgsea_GO_Po, fgsea_GO_Neg)


ggplot(fgsea_CP, aes(x = reorder(pathway, NES), y = NES, size = size, color = pval)) +
  geom_point(alpha = 0.9) +  # Use points, with some transparency
  scale_color_gradient(low = "orange", high = "blue") +  # Color gradient from high to low P-values
  scale_size(range = c(4, 10), name = "Size") +
  coord_flip() +  # Flip coordinates for better label reading
  labs(x = "Pathway", 
       y = "Normalized Enrichment Score (NES)", 
       size = "Size", 
       color = "P-value") +
  ggtitle("Visualization of Pathways") +
  theme_minimal() +  # Minimalist theme
  theme(legend.position = "right")  # Adjust legend position

# -----------------------------------------------------------------------
# Create a Heatmap from specific gene ontology term

rownames(CC_IF_values)<- CC_IF_values$Genes

library(ComplexHeatmap)

IF_Genes <- CCsubsetofGOgeneSets [["GOCC_INTERMEDIATE_FILAMENT"]]

Total_IFexp1 <- Ama_Proteo %>% filter(Genes %in% IF_Genes) %>% 
  mutate(across(.cols = 5:14, .fns = ~ log2(. + 1))) 
  
# column_to_rownames("Genes") %>% 
  # select(colnames(Ama_Proteo)[grepl("Rep",colnames(Ama_Proteo))]) %>% 
  # # filter( !if_any(everything(), is.na) ) %>% 
  # # as.matrix() %>% 
  # log2()

# Assuming Total_IFexp1 is your dataframe
Total_IFexp1[, 5:14] <- log2(Total_IFexp1[, 5:14] + 1)  # Applying log2 transformation directly


# Total_IFexp[ rowSums(is.na(Total_IFexp)) == 0, ]

# Need to convert it to matrix before plotting as heatmap
Heatmap(Total_IFexp)

# Total_IFexp <- Total_IFexp %>% as.data.frame()

write.xlsx(Total_IFexp1, "Total_IFexp1.xlsx", na = "NA", rowNames = FALSE)

library(openxlsx)

write.xlsx(Total_IFexp1, file = "Total_IFexp1.xlsx", na.string = "NA", rowNames = FALSE)




# This is for GOBP_ASTROCYTE_DEVELOPMENT
AstrodDifferen_Genes <- ontologygenesets [["GOBP_ASTROCYTE_DIFFERENTIATION"]]

AstroDiff_Total_IFexp1 <- Ama_Proteo %>% filter(Genes %in% AstrodDifferen_Genes) %>% 
  mutate(across(.cols = 5:14, .fns = ~ log2(. + 1))) 

AstroDiff_Total_IFexp1 <- AstroDiff_Total_IFexp1 %>% select(3, 5:14)
AstroDev_Total_IFexp1 <- AstroDev_Total_IFexp1 %>%  select(3, 5:14)
AstroPROJE_Total_IFexp1 <- AstroPROJE_Total_IFexp1 %>%  select(3, 5:14)


write.xlsx(AstroDiff_Total_IFexp1, "AstroDiff_Total_IFexp1.xlsx", na = "NA", rowNames = FALSE)
write.xlsx(AstroDev_Total_IFexp1, "AstroDev_Total_IFexp1.xlsx", na = "NA", rowNames = FALSE)
write.xlsx(AstroPROJE_Total_IFexp1, "AstroPROJE_Total_IFexp1.xlsx", na = "NA", rowNames = FALSE)


# column_to_rownames("Genes") %>% 
# select(colnames(Ama_Proteo)[grepl("Rep",colnames(Ama_Proteo))]) %>% 
# # filter( !if_any(everything(), is.na) ) %>% 
# # as.matrix() %>% 
# log2()

# Assuming Total_IFexp1 is your dataframe
Total_IFexp1[, 5:14] <- log2(Total_IFexp1[, 5:14] + 1)  # Applying log2 transformation directly


# Total_IFexp[ rowSums(is.na(Total_IFexp)) == 0, ]

# Need to convert it to matrix before plotting as heatmap
Heatmap(Total_IFexp)
write.xlsx(Total_IFexp1, file = "Total_IFexp1.xlsx", na.string = "NA", rowNames = FALSE)

# ----------------------------------------------------------------------

CC_IF <- fgsea_CCsubsetofGOgeneSets %>%
  filter(pathway == "GOCC_INTERMEDIATE_FILAMENT") %>% 
  select(leadingEdge) %>% 
  unlist() %>% unname() 

fgsea_CCsubsetofGOgeneSets %>%
  filter(pathway == "GOCC_KERATIN_FILAMENT") %>% 
  select(leadingEdge) %>% 
  unlist() %>% unname() #%>% length()

GOCC_KERATIN_FILAMENT

CC_IF_values <- Ama_Proteo %>% 
  filter(Genes%in% CC_IF)

# ----------------------------------------------------------
fgsea_CCsubsetofGOgeneSets %>% mutate(Term_MeanFCProteins =
lapply(fgsea_CCsubsetofGOgeneSets$leadingEdge,FUN = function(x){
  table %>% dplyr::filter(ID %in% x) %>% select(logFC) %>%  
    unlist() %>%  unname() %>% mean(na.rm = TRUE) 
}) %>% unlist()  
) -> df_test



df_test %>% ggplot(aes(Term_MeanFCProteins,-log10(padj)) )+
  geom_point() 



fgsea_CCsubsetofGOgeneSets %>% mutate(
  Term_MeanFCProteins =
    lapply(leadingEdge,
      FUN = function(x) {
        table %>% dplyr::filter(ID %in% x) %>% select(logFC) %>%
          unlist() %>%  unname() %>% mean(na.rm = TRUE)
      }
    ) %>% unlist()
) -> df_test

GoVol <- rownames_to_column(df_test, var="pathway")

df_test %>% ggplot(aes(Term_MeanFCProteins, -log10(padj))) +
  geom_point() 


fgsea_ontologygenesets %>% mutate(
  Term_MeanFCProteins =
    lapply(leadingEdge,
           FUN = function(x) {
             table %>% dplyr::filter(ID %in% x) %>% select(logFC) %>%
               unlist() %>%  unname() %>% mean(na.rm = TRUE)
           }
    ) %>% unlist()
) -> df_test



df_test %>% ggplot(aes(Term_MeanFCProteins, -log10(pval))) +
  geom_point() 


EnhancedVolcano(df_test,
                lab = df_test$pathway,
                x = 'Term_MeanFCProteins',
                y = 'pval',
                title = 'Amygdala Stress',
                pCutoff = 10e-3,
                FCcutoff = 0.8,
                pointSize = 4.0,
                labSize = 3.0,
                colAlpha = 0.7,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)


# this would be for canonical pathways
fgseaRes_CanonicalPWgeneSets %>% mutate(
  Term_MeanFCProteins =
    lapply(leadingEdge,
           FUN = function(x) {
             table %>% dplyr::filter(ID %in% x) %>% select(logFC) %>%
               unlist() %>%  unname() %>% mean(na.rm = TRUE)
           }
    ) %>% unlist()
) -> df_testCano



df_testCano %>% ggplot(aes(Term_MeanFCProteins, -log10(pval))) +
  geom_point() 


EnhancedVolcano(df_testCano,
                lab = df_testCano$pathway,
                x = 'Term_MeanFCProteins',
                y = 'pval',
                title = 'Amygdala Stress',
                pCutoff = 10e-2,
                FCcutoff = 0.8,
                pointSize = 4.0,
                labSize = 2.0,
                colAlpha = 0.7,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)

# This is for regulatory targets
# this would be for canonical pathways
fgsea_regulatorytargetgenesets %>% mutate(
  Term_MeanFCProteins =
    lapply(leadingEdge,
           FUN = function(x) {
             table %>% dplyr::filter(ID %in% x) %>% select(logFC) %>%
               unlist() %>%  unname() %>% mean(na.rm = TRUE)
           }
    ) %>% unlist()
) -> df_testRegu



df_testRegu %>% ggplot(aes(Term_MeanFCProteins, -log10(pval))) +
  geom_point() 


EnhancedVolcano(df_testRegu,
                lab = df_testRegu$pathway,
                x = 'Term_MeanFCProteins',
                y = 'pval',
                title = 'Amygdala Stress',
                pCutoff = 10e-2,
                FCcutoff = 0.8,
                pointSize = 4.0,
                labSize = 2.0,
                colAlpha = 0.7,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)


# ----------------------------------------------------------------------------
ggplot(fgsea_CP_N, aes(x = reorder(pathway, NES), y = NES, size = size, color = pval)) +
  geom_point(alpha = 0.9) +  # Use points, with some transparency
  scale_color_gradient(low = "orange", high = "blue") +  # Color gradient from high to low P-values
  scale_size(range = c(4, 10), name = "Size") +
  coord_flip() +  # Flip coordinates for better label reading
  labs(x = "Pathway", 
       y = "Normalized Enrichment Score (NES)", 
       size = "Size", 
       color = "P-value") +
  ggtitle("Visualization of Pathways") +
  theme_minimal() +  # Minimalist theme
  theme(legend.position = "right")  # Adjust legend position


BiocManager::install("enrichplot")
library(enrichplot)


# Plot results for a specific gene set
gseaplot2(fgseaRes = fgsea_res, geneSetID = "ZHANG_UTERUS_C1_REGENERATIVE_UP", title = "Gene Set 1 Enrichment Plot")



# To view top pathways
head(fgseaRes)

view(fgseaRes)

# Plotting the top pathways

plotGseaTable(fgseaRes[order(fgseaRes$padj), ][1:10, ])


head(fgseaRes[order(pval), ])


plotEnrichment(geneSets[["ZHANG_UTERUS_C0_SECRETORY_STROMAL3_NPPC_HIGH_CELL"]],
               ranks) + labs(title="ZHANG_UTERUS_C0_SECRETORY_STROMAL3_NPPC_HIGH_CELL")

plotEnrichment(geneSets[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
               ranks) + labs(title="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

plotEnrichment(geneSets[["HALLMARK_MTORC1_SIGNALING"]],
               ranks) + labs(title="HALLMARK_MTORC1_SIGNALING")


plotEnrichment(geneSets[["HALLMARK_MYC_TARGETS_V1"]],
               ranks) + labs(title="HALLMARK_MYC_TARGETS_V1")

# Sort fgsea results by adjusted p-value and select the top 10
top_pathways <- fgseaRes[order(fgseaRes$padj), ][1:10, ]

# Basic GSEA plot
ggplot(top_pathways, aes(x=reorder(pathway, NES), y=NES)) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Pathway") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top 30 Pathways")



ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES)) +
  geom_point(aes(size = padj, color = padj), alpha = 0.6) + # padj as an example, replace with any relevant metric
  scale_size_continuous(range = c(3, 20)) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  xlab("Pathway") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top 10 Pathways") +
  theme_minimal() +
  theme(legend.position = "bottom")


ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES)) +
  geom_point(aes(size = size, color = size), alpha = 0.6) + # padj as an example, replace with any relevant metric
  scale_size_continuous(range = c(3, 20)) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  xlab("Pathway") +
  ylab("Normalized Enrichment Score (NES)") +
  ggtitle("Top 10 Pathways") +
  theme_minimal() +
  theme(legend.position = "bottom")
# ---------------------------------------------------------------------------------

geneSetscgp <- gmtPathways("~/Amaygdala_proteome/m2.cgp.v2023.2.Mm.symbols.gmt")

# Perform fgsea using fgseaMultilevel
fgseaRescgp <- fgsea(pathways = geneSetscgp, 
                     stats = ranks,
                     minSize = 10, 
                     maxSize = 1000)

plotEnrichment(geneSetscgp[["MARSON_BOUND_BY_FOXP3_UNSTIMULATED"]],
               ranks) + labs(title="MARSON_BOUND_BY_FOXP3_UNSTIMULATED")

# ---------------------------------------------------------------------------------
geneSetsrtg <- gmtPathways("~/Amaygdala_proteome/m3.all.v2023.2.Mm.symbols.gmt")

# Perform fgsea using fgseaMultilevel
fgseaRescgp <- fgsea(pathways = geneSetscgp, 
                     stats = ranks,
                     minSize = 10, 
                     maxSize = 1000)

plotEnrichment(geneSetscgp[["MARSON_BOUND_BY_FOXP3_UNSTIMULATED"]],
               ranks) + labs(title="MARSON_BOUND_BY_FOXP3_UNSTIMULATED")


# ---------------------------------------------------------------------------------
significant_genes <- table_nras$Genes

databases <- listEnrichrDbs()
results <- enrichr(significant_genes, databases = c('ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', 'ARCHS4_TFs_Coexp'))

head(results$ENCODE_and_ChEA_Consensus_TFs_from_ChIP_X)
head(results$ARCHS4_TFs_Coexp)

str(results)
view(results)
str(results$ENCODE_and_ChEA_Consensus_TFs_from_ChIP_X)
str(results[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]])
top_TFs <- head(results[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]][order(-results[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]$Combined.Score), ], 10)

# Create the bar plot using ggplot2
TF<- ggplot(top_TFs, aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
  geom_bar(stat = "identity", aes(fill = Term)) +
  coord_flip() +
  xlab("Transcription Factors") +
  ylab("Combined Score") +
  ggtitle("Top 10 Transcription Factors Based on Combined Score") +
  theme_minimal() +
  theme(legend.position = "none")

TF

zz<- TF_point <- ggplot(top_TFs, aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
  geom_point(aes(color = Term), size = 4) +
  coord_flip() +
  xlab("Transcription Factors") +
  ylab("Combined Score") +
  ggtitle("Top 10 Transcription Factors Based on Combined Score") +
  theme_minimal() +
  theme(legend.position = "none")
zz


cc<-TF_lollipop <- ggplot(top_TFs, aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
  geom_segment( aes(x=reorder(Term, -Combined.Score), xend=reorder(Term, -Combined.Score), y=0, yend=Combined.Score), color="grey") +
  geom_point(aes(color = Term), size = 4) +
  coord_flip() +
  xlab("Transcription Factors") +
  ylab("Combined Score") +
  ggtitle("Top 10 Transcription Factors Based on Combined Score") +
  theme_minimal() +
  theme(legend.position = "none")
cc

# Generate the fancy point plot
TF_fancy <- ggplot(top_TFs, aes(x = reorder(Term, -Combined.Score), y = Combined.Score)) +
  geom_segment(aes(x=reorder(Term, -Combined.Score), xend=reorder(Term, -Combined.Score), y=0, yend=Combined.Score),
               color="grey50", size=0.5, alpha=0.7) +
  geom_point(aes(color = Combined.Score), size = 6, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_text_repel(aes(label = Term), box.padding = 0.5, point.padding = 0.5, nudge_y = 10) +
  coord_flip() +
  xlab("Transcription Factors") +
  ylab("Combined Score") +
  ggtitle("Top 10 Transcription Factors Based on Combined Score") +
  theme_void() +
  theme(legend.position = "bottom")

# Show the plot
TF_fancy

# Make sure to select only the top 10 TFs and the relevant columns for the radar chart
radar_data <- top_TFs[, c("Term", "Combined.Score", "Adjusted.P.value")]
colnames(radar_data) <- c("TF", "Combined Score", "Adjusted P-value")

# Set the max and min for each variable
radar_data <- rbind(rep(1, 3), rep(0, 3), radar_data)

# Create the radar chart
radarchart(radar_data)

# Create the bubble chart
bubble_chart <- ggplot(top_TFs, aes(x = reorder(Term, -Combined.Score), y = Combined.Score, size = Odds.Ratio, color = Adjusted.P.value)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(5, 30)) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  xlab("Transcription Factors") +
  ylab("Combined Score") +
  ggtitle("Top 10 Transcription Factors") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Show the plot
bubble_chart

# -----------------------------------------------------------------------------
# correlation plot

replicates <- Ama_Proteo %>% 
  select(contains("Rep")) %>% 
  mutate(across(everything(), ~log2(. + 1)))

# Calculate the correlation matrix for all log2 transformed replicates
cor_replicates <- cor(replicates, use = "pairwise.complete.obs")

# Melt the correlation matrix for plotting
cor_replicates_melted <- melt(cor_replicates)

# Plotting the heatmap
ggplot(cor_replicates_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name="Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
  labs(title = "Correlation Matrix for All Log2 Transformed Replicates", x = '', y = '')

# Using ggpairs to create a plot matrix
ggpairs(replicates, 
        upper = list(continuous = wrap("cor", size = 3)), 
        lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
        diag = list(continuous = wrap("barDiag")))

# Calculate the correlation matrix for all log2 transformed replicates
cor_replicates <- cor(replicates, use = "pairwise.complete.obs")

# Plot the correlation matrix using corrplot
corrplot(cor_replicates, method = "color", type = "upper", order = "hclust",
         addCoef.col = "black", # Add correlation coefficients
         tl.col = "black", tl.srt = 45, # Text label color and rotation
         diag = FALSE, # Do not show the diagonal
         col = colorRampPalette(c("red", "white", "blue"))(200), # Custom color palette for full range
         cl.lim = c(0, 1)) # Set color scale limits from -1 to 1


# Convert correlations to absolute values to range from 0 to 1
abs_cor_replicates <- abs(cor_replicates)

# Plot the correlation matrix using corrplot
corrplot(abs_cor_replicates, method = "color", type = "upper", order = "hclust",
         addCoef.col = "black", # Add correlation coefficients
         tl.col = "black", tl.srt = 45, # Text label color and rotation
         diag = FALSE, # Do not show the diagonal
         col = colorRampPalette(c("white", "blue"))(200), # Custom color palette from 0 to 1
         cl.lim = c(0, 1)) # Set color scale limits from 0 to 1




# --------------perform PCA plot------------------------------------
# For easy PCA plotting and results extraction
# Check for and remove rows with any NA values
replicates_clean <- replicates %>%
  na.omit()

# Now perform PCA on the cleaned data
pca_result <- prcomp(replicates_clean, scale. = TRUE, center = TRUE)

# Optionally visualize the PCA results using factoextra (if installed)
# Scree plot to see variance explained by each principal component
fviz_eig(pca_result)

# Biplot of the first two principal components
fviz_pca_biplot(pca_result, 
                geom = c("point", "text"),  # show points and variable names
                addEllipses = TRUE,         # add confidence ellipses
                label = "var",              # label variables (loadings)
                repel = TRUE)               # avoid text overlapping


fviz_pca_var(pca_result,
             geom = "arrow",    # Use arrows to represent loadings (variables)
             label = "var",     # Labels the variables
             col.var = "contrib",  # Color by contribution to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),  # Customize color gradient
             repel = TRUE,       # Avoid text overlapping
             arrow.size = 0.5,   # Adjust arrow size
             main = "PCA Variables Plot")  # Main title of the plot


# Using factoextra to plot variable loadings in a rectangular 2D plot
fviz_pca_var(pca_result,
             geom = "arrow",  # Use arrows to represent loadings (variables)
             label = "var",   # Labels the variables
             col.var = "contrib",  # Color by contribution to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),  # Customize color gradient
             repel = TRUE,    # Avoid text overlapping
             arrow.size = 0.5, # Adjust arrow size
             xlim = c(-1, 1),  # Set x-axis limits
             ylim = c(-1, 1)) + # Set y-axis limits
  ggtitle("PCA Variables Plot") +  # Set the title
  theme_minimal() +  # Use a minimal theme
  theme(axis.text = element_text(size = 12),  # Customize text size
        axis.title = element_text(size = 14))  # Customize axis title size