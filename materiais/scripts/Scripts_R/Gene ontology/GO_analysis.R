############################################################
# SCRIPT: ANALISE DE GENE ONTOLOGY (GO)
# ADAPTADO PARA O ARQUIVO: degs.csv
#
# OBJETIVO:
#   - Ler a tabela de DEGs
#   - Usar a coluna Symbol como identificador
#   - Separar genes upregulated e downregulated
#   - Rodar enriquecimento GO com clusterProfiler
#   - Exportar tabelas e gráficos
############################################################


##############################
# PACOTES NECESSARIOS
##############################

# Se precisar instalar:
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "AnnotationDbi",
#                        "org.Hs.eg.db", "enrichplot", "AnnotationHub"))
# install.packages(c("dplyr", "readr", "ggplot2"))

library(dplyr)
library(readr)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

# Como seus genes parecem ser humanos (ex.: FCRL1, TRDC, IGLV3-10),
# vamos usar o banco de anotação humano.
library(org.Hs.eg.db)

# Caso não seja humano, biscar por Db mais apropriado para sua espécie.
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, "OrgDb")
query(ah, c("Escherichia", "OrgDb"))


##############################
# DEFINIR DIRETORIO DE TRABALHO
##############################

setwd("/Users/felipemesquita/Downloads/") 

##############################
# LER A TABELA DE DEGs
##############################

# Vamos ler o arquivo mantendo a estrutura original.
# Como há colunas com vírgula decimal, não vamos depender delas agora.
# Para GO, o mais importante será Symbol, sig, DEG_up e DEG_down.

deg <- read_csv("degs.csv", show_col_types = FALSE)

# Verificar estrutura
glimpse(deg)
head(deg)


##############################
# CONFERIR COLUNAS IMPORTANTES
##############################

# Colunas esperadas no seu arquivo:
# Symbol
# adj.P.Val
# DEG_up
# DEG_down
# sig

colnames(deg)


##############################
# LIMPAR A COLUNA DE GENES
##############################

# Remover genes vazios ou NA
deg <- deg %>%
  filter(!is.na(Symbol), Symbol != "")

# Remover duplicatas simples por símbolo, se existirem
deg <- deg %>%
  distinct(Symbol, .keep_all = TRUE)


##############################
# DEFINIR OS GRUPOS DE INTERESSE
##############################

# Aqui vamos confiar nas colunas já prontas da sua tabela:
# DEG_up == TRUE
# DEG_down == TRUE
# sig == "significant"

deg_sig <- deg %>%
  filter(sig == "significant")

deg_up <- deg_sig %>%
  filter(DEG_up == TRUE)

deg_down <- deg_sig %>%
  filter(DEG_down == TRUE)

# Universo/background:
# idealmente usar todos os genes testados na análise
background_genes <- unique(deg$Symbol)

cat("Genes totais no arquivo:", length(background_genes), "\n")
cat("Genes significativos:", nrow(deg_sig), "\n")
cat("Genes upregulated:", nrow(deg_up), "\n")
cat("Genes downregulated:", nrow(deg_down), "\n")


##############################
# CONVERTER SYMBOL -> ENTREZID
##############################

# O clusterProfiler pode trabalhar com diferentes tipos de IDs,
# mas uma estratégia clássica é converter os símbolos para ENTREZID.

deg_up_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(deg_up$Symbol),
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "SYMBOL"
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

deg_down_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(deg_down$Symbol),
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "SYMBOL"
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

background_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(background_genes),
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "SYMBOL"
) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

genes_up_entrez <- unique(deg_up_map$ENTREZID)
genes_down_entrez <- unique(deg_down_map$ENTREZID)
universe_entrez <- unique(background_map$ENTREZID)

cat("Genes UP com ENTREZID:", length(genes_up_entrez), "\n")
cat("Genes DOWN com ENTREZID:", length(genes_down_entrez), "\n")
cat("Genes no universo com ENTREZID:", length(universe_entrez), "\n")


##############################
# ENRIQUECIMENTO GO - UPREGULATED
##############################

go_up_bp <- enrichGO(
  gene          = genes_up_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.20,
  minGSSize     = 5,
  maxGSSize     = 500,
  readable      = TRUE
)

go_up_cc <- enrichGO(
  gene          = genes_up_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.20,
  universe      = universe_entrez,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)

go_up_mf <- enrichGO(
  gene          = genes_up_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.20,
  universe      = universe_entrez,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)


##############################
# ENRIQUECIMENTO GO - DOWNREGULATED
##############################

go_down_bp <- enrichGO(
  gene          = genes_down_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.20,
  universe      = universe_entrez,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)

go_down_cc <- enrichGO(
  gene          = genes_down_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.20,
  universe      = universe_entrez,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)

go_down_mf <- enrichGO(
  gene          = genes_down_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.10,
  universe      = universe_entrez,
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)


##############################
# EXPORTAR TABELAS
##############################
outdir <- "/Users/felipemesquita/Downloads/"

write_csv(as.data.frame(go_up_bp),   file.path(outdir, "GO_up_BP.csv"))
write_csv(as.data.frame(go_up_cc),   file.path(outdir, "GO_up_CC.csv"))
write_csv(as.data.frame(go_up_mf),   file.path(outdir, "GO_up_MF.csv"))

write_csv(as.data.frame(go_down_bp), file.path(outdir, "GO_down_BP.csv"))
write_csv(as.data.frame(go_down_cc), file.path(outdir, "GO_down_CC.csv"))
write_csv(as.data.frame(go_down_mf), file.path(outdir, "GO_down_MF.csv"))

##############################
# VISUALIZAR TOP TERMOS
##############################

head(as.data.frame(go_up_bp))
head(as.data.frame(go_up_cc))
head(as.data.frame(go_up_mf))

head(as.data.frame(go_down_bp))
head(as.data.frame(go_down_cc))
head(as.data.frame(go_down_mf))


##############################
# DOTPLOTS
#############################


dotplot(go_up_bp, showCategory = 10) +
  ggtitle("GO Biological Process - Upregulated genes")
  


dotplot(go_down_bp, showCategory = 15) +
  ggtitle("GO Biological Process - Downregulated genes")

##############################
# BARPLOTS
##############################


barplot(go_up_bp, showCategory = 10) +
  ggtitle("GO Biological Process - Upregulated genes")



##############################
# CNETPLOT
##############################


cnetplot(go_up_bp, showCategory = 5)


##############################
# EMAPPLOT
##############################

# Calcula similaridade entre termos antes de usar emapplot

go_up_bp_sim <- pairwise_termsim(go_up_bp)

emapplot(go_up_bp_sim, showCategory = 15)

