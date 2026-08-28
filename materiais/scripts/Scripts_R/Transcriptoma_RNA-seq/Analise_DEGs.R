############################################################
## Programa de Pós-Graduação em Microbiologia Médica
## Disciplina: Bioinformática
## Prof. Dr. Felipe P. Mesquita
# Expressão diferencial com Limma-Voom (RNA-seq)
############################################################

############################
## Instalar os pacotes que serão usados nesse script. Cada pacote desse carrega funções diferentes.
## Para ver essas funções você pode consultar no help ou no doc do pacote no bioconductor.
############################

BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("ggrepel")
BiocManager::install("tidyverse")

############################
## Para que ela esteja funcional, você precisa carregar ela com a função "library"
############################

library(tidyverse) # Conjunto de pacotes para manipulação e visualização de dados.
library(edgeR) # Usado para organizar counts, filtrar baixa expressão e normalizar por TMM.
library(limma) # Usado para voom, modelo linear e estatística moderada.
library(ggrepel)  # Usado para inserir nomes de genes no volcano plot sem sobreposição excessiva.


############################
## subir a tabela de counts. Aqui eu escolhi uma aleatória dentro do GEO. Mas se for os seus dados melhor ainda.
############################
# função "setwd" para configurar o diretório de trabalho
# verifique o caminho da sua pasta de trabalho. Não esqueça de mudar.
setwd("/Users/felipemesquita/Downloads/") 

# Aqui você faz o upload do arquivo tsv dos counts. Função "read.delim". Ela lê arquivos tabulares separados por tabulação.
## comment.char = "#": Ignora linhas iniciadas com "#", comuns nos arquivos do featureCounts.
# Vale ressaltar que outras funções podem ser utilizadas, como "read_table" do pacote "readr"
# dei o nome da variável de fc, mas você pode dar o nome que quiser, mas lembre disso nos próximos comandos
fc <- read.delim("GSE297704_featureCounts.tsv", comment.char = "#", check.names = FALSE)

dim(fc) # Quantas colunas e linhas da matriz tem?
head(fc) # cabeçalho + primeiras linhas da matriz

############################################################
## Separar apenas os counts. Quero só as colunas dos counts.
############################################################

# Aqui eu estou atribuindo uma nova variável. Lembre que o sinal <- é de atribuição de funções/valores/vetores/etc a uma variável.
# quem decide o nome das variáveis é você!
# Como eu não queria essas colunas dentro da matriz, então removemos.
# A função "c" cria vetores. Vetores são sequencia de elementos do mesmo tipo. Aqui no caso é do tipo caracter.
annotation_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

#Aqui já aparece uma nova função nativa do R: "setdiff".
# retorna os elementos que estão no primeiro vetor, mas não estão no segundo.
# Vetor 1: "fc"
# Vetor 2: "annotation_columns"
## Portanto, count_columns recebe todas as colunas que NÃO são colunas de anotação. Essas colunas correspondem às amostras.
count_columns <- setdiff(colnames(fc), annotation_columns)

## O operador [ , ] seleciona linhas e colunas de um objeto. [linhas,colunas]
## fc[, count_columns]: seleciona todas as linhas e apenas as colunas dos counts.
counts <- fc[, count_columns]

## A função "rownames" define os nomes das linhas da matriz. Aqui cada linha passa a ser identificada pelo Geneid.
rownames(counts) <- fc$Geneid

## as.matrix() converte o data.frame em matriz.
## Pacotes como edgeR e limma trabalham melhor com matriz numérica para counts.
## Vale a pena lembrar dos outros "as.": as.numeric, as.vetor, as.caracter. Com eles, voce transforma o tipo de variável.
counts <- as.matrix(counts)

## storage.mode() define o tipo de dado armazenado na matriz. Contagens de RNA-seq devem ser números inteiros (Integer).
storage.mode(counts) <- "integer"


############################################################
## Criar metadados das amostras. Não precisaria se o df fosse do meu estudo pois já teria isso tabulado.
## Dá pra fazer no excel ou no nano do bash.
## Como é um estudo aleatório do GEO, tive que criar.
############################################################

## colnames(counts) recupera os nomes das amostras que estão nas colunas da matriz de counts.
## Vejam que o nome das funções são sugestivas: "colnames" - column names ou nome das colunas.
sample_names <- colnames(counts)

## data.frame() cria uma tabela de metadados.
## stringsAsFactors = FALSE: impede que textos sejam convertidos automaticamente em fatores.
## Essa tabela terá inicialmente uma coluna chamada sample, contendo o nome completo de cada amostra.
sample_info <- data.frame(sample = sample_names, stringsAsFactors = FALSE)

##COnfira ali do lado no enviroment como está o nome das amostras na variável "sample_info" (clica na tabelinha do lado)
# gsub() faz substituição de texto para limpar esse nome.
## pattern = "\\.sorted\\.bam$": procura nomes que terminam com ".sorted.bam".
## replacement = "": remove esse trecho do nome.
## x = sample_info$sample: indica onde será feita a substituição (coluna sample)
sample_info$sample_clean <- gsub(pattern = "\\.sorted\\.bam$", replacement = "", x = sample_info$sample)

# Aqui usamos de novo o gsub().
##pattern = "[-_][0-9]+-[0-9]+$"
#remove o sufixo final que identifica replicatas.
# Exemplo: "WT_MG_1-2" vira "WT_MG"
sample_info$condition <- gsub(pattern = "[-_][0-9]+-[0-9]+$", replacement = "", x = sample_info$sample_clean)

sample_info

table(sample_info$condition)

############################################################
## Definir contraste de interesse
############################################################
## Aqui começamos a ter definição da comparação que será feita nas nossas amostras.
## Na aula eu falei que é sempre necessário definir uma amostra teste e outra comparadora.
## Aqui o operador define os dois grupos que serão comparados.
## grupo_referencia: grupo usado como base da comparação.
## grupo_teste: grupo que será comparado contra a referência.
## A interpretação final será:
## grupo_teste versus grupo_referencia.
## Portanto:
## logFC > 0 significa maior expressão no grupo_teste.
## logFC < 0 significa maior expressão no grupo_referencia.
## aqui escolhi aleatoriamente dois grupos.
grupo_referencia <- "WT_MG"
grupo_teste <- "36-8D_MG"

############################################################
## Selecionar somente as amostras dos dois grupos que acabamos de definir.
############################################################
## %in% testa se cada elemento pertence a um conjunto.
## sample_info$condition %in% c(grupo_referencia, grupo_teste) retorna TRUE para as amostras pertencentes aos dois grupos de interesse e FALSE para as demais.
## Veja bem, eu só fiz isso pq esse estudo tinha 4 grupos diferentes e eu queria analisar só 2 para demonstrar.

amostras_selecionadas <- sample_info$condition %in% c(grupo_referencia, grupo_teste)

sample_info_2groups <- sample_info[amostras_selecionadas, ]

counts_2groups <- counts[, sample_info_2groups$sample]

dim(counts_2groups)

sample_info_2groups

############################################################
## Criar fator de condição
############################################################
## factor() transforma a variável condition em variável categórica.
## levels = c(grupo_referencia, grupo_teste): define a ordem dos grupos.
## O primeiro nível é a referência.
## O segundo nível é o grupo teste.
## Isso é fundamental para a interpretação do logFC no limma.

sample_info_2groups$condition <- factor(sample_info_2groups$condition, levels = c(grupo_referencia, grupo_teste))

sample_info_2groups$condition


############################################################
## Criar objeto DGEList do edgeR
############################################################
## DGEList() cria um objeto próprio do edgeR para dados de RNA-seq.
## Esse objeto armazena:
## - matriz de contagens;
## - tamanho das bibliotecas;
## - fatores de normalização;
## Ele é do tipo lista e por isso fica todo estranho ali no enviroment ->

dge <- DGEList(counts = counts_2groups)

dge$samples
#Repare que o norm.factors estão todos igual a 1. Isso significa que os reads não estão normalizados.
#A normalização a seguir com o calcNormFactors vai compensar amostras que foram sequenciadas com profundidades diferentes

############################################################
## Filtrar genes com baixa expressão
############################################################
## filterByExpr() remove genes com expressão muito baixa (poucos counts)
## y = dge: indica o objeto com as contagens. Dá uma olhada ali do lado onde está o dge ->
## group = sample_info_2groups$condition: informa a qual grupo cada amostra pertence.
keep <- filterByExpr(y = dge, group = sample_info_2groups$condition)

table(keep)
# 7 genes genes removidos por quantidade de counts.

## Filtra o objeto DGEList mantendo apenas os genes expressos.
## keep.lib.sizes = FALSE: recalcula os tamanhos das bibliotecas após a filtragem.
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

dim(dge_filtered)

############################################################
## Normalização TMM
############################################################
## calcNormFactors() calcula fatores de normalização.
## method = "TMM": usa o método Trimmed Mean of M-values.
## A normalização TMM corrige diferenças de composição entre bibliotecas,
## por exemplo quando uma amostra tem muitos reads concentrados em poucos genes altamente expressos.
## Veja que agora eu chamo para essa função o dge_filtered
dge_filtered <- calcNormFactors(object = dge_filtered, method = "TMM")

dge_filtered$samples
#Veja que agora o norm.factors mudou, indicando a compensação do tamanho da biblioteca.

############################################################
## Criar matriz de desenho experimental
############################################################
## model.matrix() cria a matriz de desenho experimental.
## ~ condition: indica que queremos modelar a expressão dos genes em função da condição experimental.
## Eii.. olha o dataframe sample_info_2groups para ver de onde vem esse "condition".
design <- model.matrix(~ condition, data = sample_info_2groups)

design

############################################################
##  Transformação voom (count para log2-CPM)
############################################################
## voom transforma os counts em log2-CPM (counts-per-million) e estima pesos observacionais para cada gene/amostra.
## counts = dge_filtered: usa o objeto filtrado e normalizado.
## design = design: informa o desenho experimental (Aquele que acabamos de criar aqui em cima)
## plot = TRUE: mostra o gráfico da tendência média-variância.
voom_data <- voom(counts = dge_filtered, design = design, plot = TRUE)

############################################################
## Ajustar modelo linear com limma
############################################################
## lmFit() ajusta um modelo linear para cada gene.
## object = voom_data: usa os dados transformados pelo voom.
## design = design: usa a matriz de desenho experimental.
## O resultado é um modelo estatístico gene a gene.
fit <- lmFit(object = voom_data, design = design)

############################################################
## Aplicar eBayes: estabiliza as estimativas de variância
############################################################
## eBayes aplica estatística Bayesiana.
## Essa etapa modera as estimativas de variância entre genes.
## O resultado melhora a estabilidade dos testes estatísticos.
## Sim, é a mesma variância de sempre (Calculada a partir da média e do desvio).
# Efetivamente, é aqui que a estatística comparativa é aplicada para ter o p-value
fit <- eBayes(fit)

############################################################
## Extrair resultados do contraste
############################################################
## topTable() extrai a tabela de resultados estatísticos que acabamos de fazer aqui em cima.
## fit = fit: usa o modelo ajustado após eBayes.
## coef = 2: usa o segundo coeficiente da matriz design.
## Neste modelo com dois grupos, o coeficiente 2 representa: grupo_teste versus grupo_referencia.
## number = Inf: retorna todos os genes, não apenas os primeiros.
## adjust.method = "BH": aplica correção de múltiplos testes pelo método Benjamini-Hochberg, gerando FDR.
## sort.by = "P": ordena os genes pelo p-valor.

resultado <- topTable(
  fit = fit,
  coef = 2,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

## rownames(resultado) contém os IDs dos genes.
## Aqui criamos uma coluna explícita chamada gene_id para usar na hora do volcano plot.
resultado$gene_id <- rownames(resultado)

head(resultado)

############################################################
## Classificar genes diferencialmente expressos
############################################################
## Define o limiar de significância estatística.
## alpha = 0.05 significa FDR menor que 5%.
# Aqui você pode ser mais rigoroso e colocar a = 0.01 e logfc_cutoff = 2, por exemplo.

alpha <- 0.05

logfc_cutoff <- 1

## Aqui criamos uma nova coluna chamada regulation e todas as "células vão ter o valor "not significant"
# Por enquanto...
resultado$regulation <- "Not significant"

## Classifica como Downregulated os genes com: FDR < 0.05 e logFC <= -1.
## Esses genes estão mais expressos no grupo_referencia ou reduzidos no grupo_teste.
resultado$regulation[resultado$adj.P.Val < alpha & resultado$logFC <= -logfc_cutoff] <- "Downregulated"

## Classifica como Upregulated os genes com: FDR < 0.05 e logFC >= 1.
## Esses genes estão mais expressos no grupo_teste em relação ao grupo_referencia.
resultado$regulation[resultado$adj.P.Val < alpha & resultado$logFC >= logfc_cutoff] <- "Upregulated"


table(resultado$regulation)


############################################################
## Preparar dados para volcano plot
############################################################

volcano_data <- resultado #aqui estamos só duplicando a matriz. uma cópia

## -log10(FDR) transforma o FDR para facilitar visualização.
## Quanto maior o valor de -log10(FDR), mais significativo é o gene.
## pmax() evita problema com FDR igual a zero.
## .Machine$double.xmin é o menor número positivo representável pelo R.
volcano_data$neg_log10_FDR <- -log10(pmax(volcano_data$adj.P.Val, .Machine$double.xmin))

## Cria uma coluna vazia para os rótulos dos genes no gráfico.
## Inicialmente nenhum gene será rotulado.
volcano_data$label <- NA

## Identifica quais genes são diferencialmente expressos.
## o sinal != significa "diferente de", ou seja, só vai ser considerado o que for diferente de not significant
genes_significativos <- volcano_data$regulation != "Not significant"

## which() retorna os índices/posições dos genes significativos.
indices_significativos <- which(genes_significativos)

# order() ordena os genes significativos pelo FDR.
## Os menores valores de FDR aparecem primeiro
indices_ordenados <- indices_significativos[order(volcano_data$adj.P.Val[indices_significativos])]

## head() seleciona os primeiros 20 genes mais significativos.
## Aqui vai definir a label no volcano. quantos genes serão marcos com seu nome.
## depois faça o teste com outros numeros
top_labels <- head(indices_ordenados, 20)

# Adiciona o nome dos genes mais significativos na coluna label.
## Apenas esses genes serão rotulados no volcano plot.
volcano_data$label[top_labels] <- volcano_data$gene_id[top_labels]

############################################################
## Volcano plot
############################################################

volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_FDR, color = regulation)) +
geom_point(alpha = 0.75, size = 1.8) +
geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed") +
geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
geom_text_repel(aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "steelblue", "Not significant" = "grey70")) +
labs(title = paste0(grupo_teste, " vs ", grupo_referencia), x = "log2 Fold Change", y = "-log10 FDR", color = "Regulation") +
theme_bw(base_size = 14)

volcano_plot

############################################################
## Salvar tabela apenas com DEGs
############################################################
## Seleciona apenas genes classificados como Upregulated ou Downregulated.
## Genes "Not significant" são removidos.
degs <- resultado[resultado$regulation != "Not significant",]

degs <- degs[order(degs$adj.P.Val), ]

head(degs)

# é isso meu povo!
# Guardem esse script para quando precisarem analisar genes diferencialmente expressos.
# Dica: Pesquise também como fazer analise de vias com gene ontology e GSEA.

