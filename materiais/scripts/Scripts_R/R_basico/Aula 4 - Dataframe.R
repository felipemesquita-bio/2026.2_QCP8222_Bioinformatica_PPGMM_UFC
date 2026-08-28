# DATAFRAME
# Objetivo: trabalhar com base de dados

#Definir area de trabalho
setwd("/Users/felipemesquita/Documents/Bioinformática UFC/scripts_basicos_R/")

#importanto os dados
## txt
#tabela <- read.table("tabela.txt", 
#                     header = TRUE, sep = "")

## csv
#tabela <- read.csv("tabela.csv", 
#                   header = TRUE, sep = ";")

library(readxl)
df <- read_excel("DataFrame.xlsx")

#Analisando um dataframe
View(df)
dim(df) # Dimensões (linhas / colunas)
str(df) # Mais informações sobre cada linha
names(df) # Nomes das colunas
head(df) # 5 primeiras linhas
tail(df) # 5 ultimas linhas
summary(df)

#Selecionando variaveis
df
df[3]
df$PAS

Col1 <- df[3]
Col2 <- df$PAS

class(df$PAS)
class(Col2)
class(Col1)

# Modificando df

# Excluindo variavel
df$PAS <- NULL
df

# Alterando o tipo de variavel dentro do df
df$Glasgow
summary(df$Glasgow)
df$Glasgow <- as.factor(df$Glasgow)
df$Glasgow
summary(df$Glasgow)



###### Checando colunas específicas



####### Media, Mediana, Max, Min (variaveis numericas)
table(df$Óbito) # tabela de frequencias
mean(df$PAS)
median(df$PAS)
var(df$PAS) # Variância
sd(df$PAS) # Desvio padrão

summary(df)


####### Dplyr 
install.packages("dplyr")
library("dplyr")

# arrange and desc : Ordenar os dados de acordo com alguma variável 

df_obito <- arrange(df, Óbito)
View(df_obito)

df_desc <- arrange(df, desc(PAS))
View(df_desc)


# Select : selecionando colunas
#selecionar determinadas colunas dentro de um data.frame

df_select <- select(df, Paciente, PAS, Óbito)

df_select <-  select(df, -PAS)

df_select <- df %>%
  select(Paciente, PAS, Óbito)


# Filter : filtrando casos

df_alive <- df %>%
  filter(Óbito == "N")

df_high <- df %>%
  filter(PAS > 100)

# operadores & (and) e | (or)

df_or <- df %>%
  filter(Glasgow > 5 | Glasgow < 5)
View(df_or)

df_and <- df %>%
  filter(Óbito == "N" & PAS > 70)
View(df_and)

# Mutate (criando novas colunas)

df$novacoluna <- df$Glasgow + df$PAS

df_mutate <- df %>%
  mutate(ratio = PAS / Glasgow,
         fd = `Dosagem Proteína Y`*100)

View(df_mutate)


####### Plots 

# boxplot
boxplot(df$`Dosagem Proteína Y`)
boxplot(df$`Dosagem Proteína Y` ~ df$Óbito)

#histograma
hist(df$PAS)

# correlacao
plot(df$`Dosagem Proteína Y` ~ df$PAS)

####### Finalizando 
rm(list = ls()) # limpar lista do environment


#Quer ver mais como fazer mais plots?
# Consultar pacote ggplot2 