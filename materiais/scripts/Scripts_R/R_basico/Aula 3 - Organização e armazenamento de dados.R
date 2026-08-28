
#OBJETIVO: entender como o R organiza e armazena os dados

#Armazenamento numérico

Salario <- 2200.89
Horas <- 220

SH <- Salario/Hora
SHi <- as.integer(Salario/Horas)
SHi * 5
SHr <- round(Salario/Horas)

Numeros1 <- c(Salario,Horas)

#Armazenamento de caracteres

Nome_1 <- "Eduardo Abreu"
Nome_2 <- "Amanda Lopes"
Idade <- "25"
Idade + 5

Nomes <- Nome_1 + Nome_2
Nomes <- c(Nome_1,Nome_2)

Nomes
Nomes[1]
Nomes[2]

#Armazenamento de fatores
Cargahoraria <- c(220,220,150,100,100)
summary(Cargahoraria)

Cargahoraria <- as.factor(Cargahoraria)
summary(Cargahoraria)

Cargahoraria <- as.factor(c(220,220,150,100,100))
summary(Cargahoraria)


#Armazenamento Logico

L1 <- Salario > Horas
L1

L2 <- Salario < Horas
L2

Logico <- TRUE
Logico1 <- "TRUE"

Logico2 <- c(1,TRUE,3)
Logico2


# VETORES - estrutura basica de dados
# Uma sequencia de dados do mesmo tipo

# Vetor de caracteres
is.vector(Names)
mode(Names)

# Vetor numerico
is.vector(Horas)
mode(Horas)

# Vetor Logico
is.vector(L1)
mode(L1)

# Vetor é uma estrutura de um dado com um unico tipo. 
# P.Ex: L1 só tem um dado com uma sequencia logica; Horas só tem variavel numerica

# LISTAS
# é como se fosse um Vetor com tipos de dados diferentes

a <- c(1,2,3,4,5)
a

b <- c(1,"2",3,4,5)
b

is.list(a)
is.list(b)

is.vector(a)
is.vector(b)

b <- list(10,"2",8)
is.list(b)
mode(b)
str(b)

e <- list(c(10,6,51,5),"2",8)
e
str(e)
e[[1]][1]


# MATRIZES - duas dimesões de um tipo de dado (semelhante aos vetores)
m <- matrix(1:9, nrow = 3)
m
mode(m)
class(m)

#filtrar dados na matriz
m[1,3]

#modificar dados na matriz
m[1,3] <- "a"

mode(m)
class(m)

