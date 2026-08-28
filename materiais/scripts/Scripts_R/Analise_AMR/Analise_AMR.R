############################################################
## Programa de Pós-Graduação em Microbiologia Médica
## Disciplina: Bioinformática
## Prof. Dr. Felipe P. Mesquita
## Análise terciária de genes de resistência antimicrobiana
## Formato de entrada: tabela TSV do AMRFinderPlus
############################################################

##############################
## 1. INSTALAR E CARREGAR PACOTES
##############################

# Execute a instalação somente uma vez, se necessário:
# install.packages("tidyverse")

library(tidyverse)

##############################
## 2. DEFINIR ARQUIVO E LIMIARES
##############################

# Coloque o script e a tabela na mesma pasta.
# No RStudio, abra o projeto ou defina essa pasta como diretório de trabalho.
arquivo_amr <- "amr_resultados_exemplo.tsv"

# Limiares didáticos. Eles devem ser ajustados conforme a pergunta,
# o organismo, a ferramenta e o protocolo empregado.
cobertura_minima <- 90
identidade_minima <- 90

##############################
## 3. IMPORTAR A TABELA
##############################

# O AMRFinderPlus produz uma tabela separada por tabulação (TSV).
# check.names = FALSE preserva nomes como "% Identity to reference".
amr_bruto <- read.delim(
  arquivo_amr,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dim(amr_bruto)
head(amr_bruto)
colnames(amr_bruto)

##############################
## 4. CONFERIR COLUNAS NECESSÁRIAS
##############################

colunas_necessarias <- c(
  "Name", "Element symbol", "Element name", "Type", "Class",
  "Subclass", "Method", "% Coverage of reference",
  "% Identity to reference", "Contig id", "Start", "Stop"
)

colunas_ausentes <- setdiff(colunas_necessarias, colnames(amr_bruto))

if (length(colunas_ausentes) > 0) {
  stop(
    "A tabela não contém as seguintes colunas: ",
    paste(colunas_ausentes, collapse = ", ")
  )
}

##############################
## 5. SELECIONAR E RENOMEAR COLUNAS
##############################

amr <- amr_bruto %>%
  transmute(
    amostra = Name,
    gene = `Element symbol`,
    descricao = `Element name`,
    tipo = Type,
    classe = Class,
    subclasse = Subclass,
    metodo = Method,
    cobertura = as.numeric(`% Coverage of reference`),
    identidade = as.numeric(`% Identity to reference`),
    contig = `Contig id`,
    inicio = as.numeric(Start),
    fim = as.numeric(Stop)
  )

glimpse(amr)

##############################
## 6. FILTRAR ACHADOS DE AMR
##############################

# Mantemos somente elementos classificados como AMR e que atingem
# os limiares definidos de cobertura e identidade.
amr_filtrado <- amr %>%
  filter(
    tipo == "AMR",
    cobertura >= cobertura_minima,
    identidade >= identidade_minima,
    !is.na(gene),
    gene != ""
  ) %>%
  distinct(amostra, gene, contig, inicio, fim, .keep_all = TRUE)

cat("Achados antes do filtro:", nrow(amr), "\n")
cat("Achados após o filtro:", nrow(amr_filtrado), "\n")
cat("Número de amostras:", n_distinct(amr_filtrado$amostra), "\n")
cat("Genes distintos:", n_distinct(amr_filtrado$gene), "\n")

##############################
## 7. RESUMOS DESCRITIVOS
##############################

resumo_amostras <- amr_filtrado %>%
  count(amostra, name = "numero_de_genes") %>%
  arrange(desc(numero_de_genes))

resumo_genes <- amr_filtrado %>%
  count(gene, classe, subclasse, name = "numero_de_amostras") %>%
  arrange(desc(numero_de_amostras), gene)

resumo_classes <- amr_filtrado %>%
  count(amostra, classe, name = "numero_de_genes")

resumo_amostras
resumo_genes
resumo_classes

##############################
## 8. GRÁFICO: GENES POR AMOSTRA
##############################

grafico_amostras <- ggplot(
  resumo_amostras,
  aes(x = reorder(amostra, numero_de_genes), y = numero_de_genes)
) +
  geom_col(fill = "#167d83", width = 0.72) +
  coord_flip() +
  labs(
    title = "Genes de resistência detectados por amostra",
    x = "Amostra",
    y = "Número de genes"
  ) +
  theme_minimal(base_size = 12)

grafico_amostras

##############################
## 9. GRÁFICO: CLASSES POR AMOSTRA
##############################

grafico_classes <- ggplot(
  resumo_classes,
  aes(x = amostra, y = numero_de_genes, fill = classe)
) +
  geom_col(width = 0.75) +
  labs(
    title = "Perfil de classes de resistência",
    x = "Amostra",
    y = "Número de genes",
    fill = "Classe"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grafico_classes

##############################
## 10. MATRIZ DE PRESENÇA E AUSÊNCIA
##############################

# Se um gene aparece uma ou mais vezes na amostra, recebe valor 1.
presenca_long <- amr_filtrado %>%
  distinct(amostra, gene) %>%
  mutate(presenca = 1L)

matriz_presenca <- presenca_long %>%
  pivot_wider(
    names_from = gene,
    values_from = presenca,
    values_fill = 0
  ) %>%
  arrange(amostra)

matriz_presenca

##############################
## 11. HEATMAP DE GENES POR AMOSTRA
##############################

grafico_heatmap <- ggplot(
  presenca_long,
  aes(x = gene, y = amostra, fill = factor(presenca))
) +
  geom_tile(color = "white", linewidth = 0.7) +
  scale_fill_manual(
    values = c("1" = "#167d83"),
    labels = c("1" = "Presente"),
    name = "Detecção"
  ) +
  labs(
    title = "Presença de genes de resistência",
    x = "Gene",
    y = "Amostra"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

grafico_heatmap

##############################
## 12. EXPORTAR RESULTADOS
##############################

dir.create("resultados_amr", showWarnings = FALSE)

write_csv(amr_filtrado, "resultados_amr/achados_amr_filtrados.csv")
write_csv(resumo_amostras, "resultados_amr/resumo_por_amostra.csv")
write_csv(resumo_genes, "resultados_amr/resumo_por_gene.csv")
write_csv(matriz_presenca, "resultados_amr/matriz_presenca_ausencia.csv")

ggsave(
  "resultados_amr/genes_por_amostra.png",
  grafico_amostras,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "resultados_amr/classes_por_amostra.png",
  grafico_classes,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  "resultados_amr/heatmap_genes_amostras.png",
  grafico_heatmap,
  width = 10,
  height = 6,
  dpi = 300
)

############################################################
## INTERPRETAÇÃO
############################################################

# A detecção genotípica de um determinante de resistência não confirma,
# isoladamente, o fenótipo de resistência. A interpretação deve considerar
# espécie, integridade e expressão do gene, mecanismo, método de detecção,
# versão do banco de dados e teste de sensibilidade aos antimicrobianos.

cat("Análise concluída. Resultados salvos em: resultados_amr/\n")
