# OBJETIVO: primeiro contato com a linguagem R
#colocar enviroment em GRID

# Navegação e exemplos
#criando variáveis
#para executar clique sobre a linha ou selecione o conjunto de linhas e pressionar crtl + enter

#operações sem atribuir variavel não SALVA
5+5
4*9


a <- 10
a

b
b = 10
b
B

#utilizaremos <- na criação de variáveis e = nas atribuições de funções (possuem a mesma funções, mas com o = com algumas particularidades)
#para padronizar utilizaremos sempre <-
#variável <- valor

a <- 10
b <- 5

c <- a + b

#Criamos uma string, variavel não-numérica (sempre entre aspas)
a <- "Emerson"
b <- "Celina"

c <- a + b
c

#Funçoes (Ver packages)
#Funções estão smpre entre parenteses
#Fazem parte de um pacote
#Pacotes contem muitas funçoes e precisam estar instaladas
#Nao sendo um pacote padrao, ele precisa ser chamado

c <- c(a,b)
c

#Help
?c

a <- c(10, 5, 15, 20)
a
a[1]
a[3]

?summary
summary(a)
summary(c)

#Funçao de um pacote nao padrao (ver rdocumentation.org)
?str_c

install.packages("stringr")
library(stringr)

?str_c

Nome <- Emerson
Sobrenome <- Silva

NomeCompleto <- str_c(Nome, Sobrenome)
NomeCompleto

NomeCompleto <- str_c(Nome, " ", Sobrenome)
NomeCompleto











