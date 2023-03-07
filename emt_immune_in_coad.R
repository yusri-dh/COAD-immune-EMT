# load library ------------------------------------------------------------

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(GSEABase)
library(edgeR)
library(singscore)
library(pheatmap)
library(NbClust)
library(psych)
library(corrplot)
library(ggplot2)
library(ggpubr)
library(randomForestSRC)
library(rstatix)

working_dir = "./"
setwd(working_dir)
genome <- "hg38"
#function to loads an RData file, and returns it
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# load data -------------------------------------------------------------
# we obtained the COAD_rna_hg38.rda from coad_data_download.R script
rna_all <- loadRData("./data/COAD_rna_hg38.rda")


# filter low expression gene
prop_expressed = rowMeans(assay(rna_all) >= 5)
keep = prop_expressed > 0.5
rna_deg = rna_all[keep, ]

assay(rna_deg, 'logTPM') = log1p(assay(rna_deg, 'tpm_unstrand'))
all_genes = rowData(rna_deg)$gene_name

# filter duplicated gene name
dups = unique(all_genes[duplicated(all_genes)])
rna_deg = rna_deg[!all_genes %in% dups,]

ensmbl_id_name_convert = rowData(rna_deg)[, c("gene_id", "gene_name")] %>%
  as_tibble() %>%
  separate(
    col = gene_id,
    remove = FALSE,
    sep = "\\.",
    into = c("id", "id2")
  ) %>%
  mutate(id2 = NULL)

rm(dups, keep)

# Which samples are Primary Tumor
samples.primary.tumour <-
  rna_deg$barcode[rna_deg$shortLetterCode == "TP"]

# which samples are solid tissue normal
samples.solid.tissue.normal <-
  rna_deg$barcode[rna_deg$shortLetterCode == "NT"]

# Normal vs tumor DEG analysis-------------------------------------------------------------
dataPrep <- TCGAanalyze_Preprocessing(object = rna_deg,
                                      cor.cut = 0.6)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.25)


dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[, samples.solid.tissue.normal],
  mat2 = dataFilt[, samples.primary.tumour],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT",
  pipeline = "edgeR"
)

rm(dataPrep,
   dataNorm,
   dataFilt,
   samples.primary.tumour,
   samples.solid.tissue.normal)
gc()

# Transcriptome signature analysis----------------------------------------------
# Load gene sets----------------------------------------------------------------
# Immune cell signatures gene set from Angelova et.al.(2015)
path_angelova = "./data/gene_set_angelova_2015.csv"
angelova_df = read_csv(path_angelova, col_names = c("gene_name", "set_name")) %>%
  mutate(set_name = str_replace_all(set_name, " ", "_"))


# M1 and M2 macrophages cell signatures gene set from Aran et.al.(2017)
path_xcell = "./data/gene_set_xcell_2017.csv"
xcell_df = read_csv(path_xcell)
m1_genes = xcell_df[c(
  "M1_BLUEPRINT_1",
  "M1_BLUEPRINT_2",
  "M1_BLUEPRINT_3",
  "M1_FANTOM_1",
  "M1_FANTOM_2",
  "M1_FANTOM_3"
)] %>% unlist() %>% unique() %>% as.character() %>% na.omit() %>%
  as_tibble() %>%
  mutate(set_name = "M1") %>%
  rename(gene_name = value)
m2_genes = xcell_df[c(
  "M2_BLUEPRINT_1",
  "M2_BLUEPRINT_2",
  "M2_BLUEPRINT_3",
  "M2_HPCA_1",
  "M2_HPCA_2",
  "M2_HPCA_3"
)] %>% unlist() %>% unique() %>% as.character() %>% na.omit() %>%
  as_tibble() %>%
  mutate(set_name = "M2") %>%
  rename(gene_name = value)
xcell_df = rbind(m1_genes, m2_genes)

# The immune subtypes gene signatures from Thorsson 2018
path_thorsson = "./data/gene_set_thorsson_2018.csv"
thorsson_df = read_csv(path_thorsson, col_names = c("gene_name", "set_name")) %>%
  mutate(set_name = str_replace_all(set_name, " ", "_"))

# The Epithelial-mesenchymal gene signatures from Thiery et.al. 2014
thiery_path = './data/Thiery_EMTsignatures.txt'
thiery_df = read.table(thiery_path, header = T) %>%
  dplyr::select(officialSymbol, epiMes_tumor) %>%
  rename(gene_name = officialSymbol) %>%
  drop_na(epiMes_tumor) %>%
  mutate(set_name = str_replace_all(epiMes_tumor, c("epi" = "Epithelial",
                                                    "mes" = "Mesenchymal"))) %>%
  dplyr::select(-epiMes_tumor)

combined_df = bind_rows(angelova_df, thiery_df, thorsson_df, xcell_df)
rm(
  angelova_df,
  thiery_df,
  thorsson_df,
  path_angelova,
  path_thorsson,
  thiery_path,
  path_xcell,
  xcell_df
)


gene_set_name <- combined_df$set_name %>% unique
gsc_list = list()
for (i in gene_set_name) {
  print(i)
  genes = combined_df %>%
    filter(set_name == i) %>%
    pull(gene_name) %>%
    unique()
  gsc_list[[i]] = GeneSet(genes, setName = i, geneIdType = SymbolIdentifier())
  rm(genes)
  rm(i)
}
gsc <- GeneSetCollection(gsc_list)

# calculate Single sample scoring of molecular phenotypes
rna = rna_deg[, rna_deg$sample_type == "Primary Tumor"]
rownames(rna) = rowData(rna)$gene_name

eranks = rankGenes(assay(rna, 'logTPM'))

coad_multiscore =  multiScore(eranks, gsc)
coad_scores = coad_multiscore$Scores

# Clustering the samples-------------------------------------------------------
# Calculate optimal number of cluster (K). We found that the optimal number of cluster = 3
# res_nblust <- NbClust::NbClust(as_tibble(t(coad_scores)),
#                                 method ="ward.D2")

distance_mat <- dist(t(coad_scores), method = 'euclidean')
Hierar_cl <- hclust(distance_mat, method = "ward.D2")
fit <- cutree(Hierar_cl, k = 3)
ordered_sample = Hierar_cl$labels[Hierar_cl$order]
group = factor(fit,
               levels = c("2", "1", "3"),
               labels = c("high", "medium", "low"))
rna$group = group

coad_scores_df = t(coad_scores) %>%
  as.data.frame() %>%
  mutate(group = group) %>%
  arrange(by = Mesenchymal)

# Hierarchical clustering visualization
pheatmap(
  coad_scores,
  main = "Immune and EMT signature scores",
  annotation_col = coad_scores_df["group"],
  cluster_rows = T,
  cluster_cols = T,
  cutree_cols = 3,
  cutree_rows = 3,
  clustering_method = "ward.D2",
  # treeheight_col=2,
  show_colnames = FALSE,
  fontsize_row = 5
)

# Compare the means of immune signature  between EMT groups--------------------------------
# create temporary dataframe for analysis
coad_scores_df_long <- t(coad_scores) %>%
  scale() %>%
  as.data.frame() %>%
  mutate(group = group) %>%
  gather(gene_set, score,-group)

# student t-test with bonferroni correction
stat.test <- coad_scores_df_long %>%
  group_by(gene_set) %>%
  t_test(score ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_y_position(step.increase = 0.25,
                 comparisons = list(c("medium", "high"), c("low", "medium"), c("low", "high")))

# visualization
f <- ggboxplot(
  coad_scores_df_long,
  x = "group",
  y = "score",
  order = c("low", "medium", "high"),
  color = "group",
  ylab = "Standardized Signature Score",
  ggtheme = theme_gray(),
  repel = T
) +
  stat_pvalue_manual(stat.test,
                     label = "p.adj.signif",
                     tip.length = 0.01,
                     size = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme(strip.text.x = element_text(size = 8), axis.title.x = element_blank())
facet(f,
      facet.by = "gene_set",
      nrow = 5,
      ncol = 8)

rm(coad_scores_df_long)

# Correlation analysis intervariables------------------------------------------
# create temporary dataframe for correlation analysis
coad_scores_df_long <- t(coad_scores) %>%
  scale() %>%
  as.data.frame() %>%
  mutate(group = group) %>%
  gather(gene_set, score, -Mesenchymal,-group)

cor.test.gene_set <- coad_scores_df_long %>%
  group_by(gene_set) %>%
  cor_test(Mesenchymal, score) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()  %>%
  select(gene_set, cor, p, conf.low, conf.high, p.adj, p.adj.signif) %>%
  arrange(desc(cor))

coad_scores_df_long$gene_set <-
  factor(coad_scores_df_long$gene_set, levels = cor.test.gene_set$gene_set)


f <- (
  ggplot(coad_scores_df_long, aes(x = Mesenchymal, y = score)) +
    geom_point(alpha = 0.5, size = 0.3) +
    geom_smooth(method = lm) +
    stat_cor(
      method = "pearson",
      aes(label = ..r.label..),
      label.y = -5,
      size = 3
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.05))) +
    ylab("Standardized Signature Score") +
    xlab("Mesenchymal signature scores")
) %>% facet(facet.by = "gene_set", ncol = 8)
f

rm(coad_scores_df_long)

# Immunomodulator correlation analysis --------------------------------------------
# load immunomodulator list from Thorsson et al. 2018.
immunomodulator_df = read.csv("./data/immunomodulator_genes.csv") %>%
  mutate(HGNC.Symbol2 = str_replace(HGNC.Symbol, "-", ".")) %>%
  mutate(label = coalesce(Immune.Checkpoint, Super.Category)) %>%
  arrange(label)

immunomodulator_df = immunomodulator_df[(immunomodulator_df$HGNC.Symbol %in% rownames(rna)), ]
inhibitory_genes = immunomodulator_df %>%
  filter(Immune.Checkpoint == "Inhibitory") %>%
  pull(HGNC.Symbol)
stimulatory_genes = immunomodulator_df %>%
  filter(Immune.Checkpoint == "Stimulatory") %>%
  pull(HGNC.Symbol)
apc_genes = immunomodulator_df %>%
  filter(Super.Category == "Antigen presentation") %>%
  pull(HGNC.Symbol) %>%
  str_replace("HLA-", "HLA.")

x = assay(rna, "logTPM") %>%
  t() %>%
  scale() %>%
  as.data.frame()
y_numeric = coad_scores_df[rownames(x), ] %>%
  mutate(group = NULL) %>%
  scale() %>%
  as.data.frame()
y_group = coad_scores_df[rownames(x), ] %>%
  pull(group)

if (identical(rownames(x), rownames(y_numeric))) {
  emt_df = data.frame(x, y_numeric)
  emt_df["group"] = y_group
}

df_long = emt_df[, c('Mesenchymal',
                     'group',
                     inhibitory_genes,
                     stimulatory_genes,
                     apc_genes)] %>%
  gather(measurement, value,-Mesenchymal, -group) %>%
  mutate(measurement = str_replace(measurement, "HLA.", "HLA-")) # synchronize the gene name with "HLA-"

cor.test.immunomod <- df_long %>%
  group_by(measurement) %>%
  cor_test(Mesenchymal, value) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  left_join(immunomodulator_df[, c("HGNC.Symbol", "label")], by = c("measurement" =
                                                                      "HGNC.Symbol"))

# select genes with adjusted P < 0.05 and correlation > 0.3
significant_cor_genes <- cor.test.immunomod %>%
  arrange(desc(cor)) %>%
  filter(p.adj < 0.05) %>%
  filter(abs(cor) > 0.3) %>%
  mutate(measurement = str_replace(measurement, "HLA-", "HLA.")) %>%
  pull(measurement)

# Identify important genes using multivariate random forest--------------------

x = subset(assay(rna, "logTPM"), subset = rownames(rna) %in% dataDEGs$gene_name) %>%
  t() %>%
  scale() %>%
  as.data.frame()
y_numeric = coad_scores_df[rownames(x), ] %>%
  mutate(group = NULL) %>%
  scale() %>%
  as.data.frame()

if (identical(rownames(x), rownames(y_numeric))) {
  emt_df = data.frame(x, y_numeric)
}



# function to run multivariate random forest
run_multivariate_rf <- function(gene_set, x_predictor, y_target) {
  y = y_target[, c("Mesenchymal", gene_set)]
  f = get.mv.formula(c("Mesenchymal", gene_set))
  mreg <- rfsrc(
    f,
    data.frame(x_predictor, y),
    ntree = 5000,
    importance = "permute",
    nodesize = 10,
    splitrule = "mahalanobis",
    nsplit = 10
  )
  imp <- get.mv.vimp(mreg,  standardize = TRUE) %>%
    data.frame()
  return(list("rf" = mreg, "importance_df" = imp))
}


immune_process_sets = c(
  "Lymphocyte_infiltration",
  "CSF1_response",
  "IFNG_response" ,
  "TGFB_response",
  "Wound_healing"
)

immune_process_rf = run_multivariate_rf(immune_process_sets, x, y_numeric)

# extract variable importance and store it in the dataframe format
immune_process_imp = immune_process_rf$importance_df %>%
  mutate(arithmetic_mean = NULL, multi = NULL) %>%
  #scale() %>%
  as.data.frame()
mean_imp = rowMeans(immune_process_imp) # calculate mean importance score
mult_imp = apply(immune_process_imp, 1, prod) # calculate multiplicative score
immune_process_imp[["mean_importance"]] = mean_imp
immune_process_imp[["importance_product"]] = mult_imp
