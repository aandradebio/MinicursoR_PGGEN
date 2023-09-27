### Bem vindos ao Minicurso de Exploração e Visualização de dados biológicos em R
### Msc. Amanda Araújo Serrão de Andrade - Laboratório de Bioinformática (LNCC/UFRJ)
### Duvidas podem ser encaminhadas por email: aandradebio@gmail.com

### Temas da Aula 2: Estatística descritiva e Pré-processamento dos dados tabulares e sequências de DNA

### Vamos carregar os pacotes utilizados na aula de hoje 
library(data.table) # install.packages("data.table")
library(ggplot2) # install.packages("ggplot2")
library(Biostrings)

# Estatística descritiva são métodos estatísticos usados para resumir, organizar e descrever um conjunto de dados. 
# A etapa de pré-processamento de dados é uma parte crítica da análise de dados biológicos que envolve a limpeza, organização e transformação dos dados brutos em uma forma adequada para análises subsequentes. 

## Sessão 1: Dados tabulares

taxonomy <- read.table(file = "tax.csv", header = TRUE, sep = ",", dec = ".") # Correto

# Sessão 1.1: Estatística Descritiva

dim(taxonomy)
print(taxonomy)
head(taxonomy)
str(taxonomy)

# Identificação de linhas duplicadas
duplicatas <- taxonomy[duplicated(taxonomy), ]

# Identificação de linhas duplicadas em uma coluna específica
duplicatas <- taxonomy[duplicated(taxonomy$species), ]
duplicatas <- taxonomy[duplicated(taxonomy$family), ]

# Identificação de duplicatas com base em várias colunas
duplicatas <- taxonomy[duplicated(taxonomy[c('species', 'family')]), ]
duplicatas <- taxonomy[duplicated(taxonomy[c('order', 'family')]), ]

# Análise das variáveis categóricas

# Contagens das variaveis categóricas

# Este comando calcula a contagem de quantas vezes cada espécie de vírus aparece no conjunto de dados. É útil para determinar quais espécies estão mais representadas e podem ser de interesse para análises posteriores.
table(taxonomy$species) # Contagem das espécies

# Quantas espécies pertencem a cada família de vírus. Isso pode ajudar a identificar quais famílias são mais diversas em termos de espécies
table(taxonomy$family) # Contagem do número de espécies por família
table(taxonomy$genus) # Contagem do número de espécies por gênero
table(taxonomy$order) # Contagem do número de espécies por ordem

# Este comando cria uma tabela de contingência cruzada que mostra como as espécies estão distribuídas entre diferentes famílias de vírus. Pode ajudar a identificar relações entre espécies e famílias.
table(taxonomy$species, taxonomy$family)

prop.table(table(taxonomy$species)) # Contagem das espécies
prop.table(table(taxonomy$family)) # Contagem do número de espécies por família
prop.table(table(taxonomy$genus)) # Contagem do número de espécies por gênero
prop.table(table(taxonomy$order)) # Contagem do número de espécies por ordem

# Esses comandos identificam o nome da espécie, família, gênero e ordem mais frequentes no conjunto de dados. Isso pode ser útil para destacar as categorias mais predominantes.
names(sort(table(taxonomy$species), decreasing = TRUE)[1]) # Nomes das categorias mais frequentes
names(sort(table(taxonomy$family), decreasing = TRUE)[1])
names(sort(table(taxonomy$genus), decreasing = TRUE)[1])
names(sort(table(taxonomy$order), decreasing = TRUE)[1])

# Análise da variável numérica

summary(taxonomy) # Resumo Estatístico de cada coluna
summary(taxonomy$num_contigs) # Resumo Estatístico da coluna numérica

# A média é a soma de todos os valores dividida pelo número de observações. 
# A mediana é o valor que divide os dados ao meio. Ambas são medidas de tendência central.

mean(taxonomy$num_contigs) # Calculo da média da coluna numérica
median(taxonomy$num_contigs) # Calculo da mediana da coluna numérica

# O desvio padrão mede a dispersão dos dados em torno da média. 
# A variância é o quadrado do desvio padrão.

sd(taxonomy$num_contigs) # Calculo do desvio padrão da coluna numérica
var(taxonomy$num_contigs) # Calculo da variância da coluna numérica

# Percentis são valores que dividem os dados em percentagens específicas. Os quartis são percentis que dividem os dados em quatro partes iguais. 

quantile(taxonomy$num_contigs, probs = 0.25) # Quartil 25
quantile(taxonomy$num_contigs, probs = 0.75) # Quartil 75

# Contagens 

table(taxonomy$num_contigs)

# Comentários sobre a coluna numérica: 
# Os vírus apresentam no mínimo 101 contigs por espécie, têm um valor máximo de 867.381 contigs e uma média de 43.657 contigs.
# O primeiro quartil (25%) está em 2.240 contigs, o que significa que 25% dos valores estão abaixo desse ponto.
# A mediana (50%) está em 6.073 contigs, o que é próximo da média, sugerindo uma distribuição aproximadamente simétrica.
# O terceiro quartil (75%) está em 23.610 contigs, o que significa que 25% dos valores estão acima desse ponto.

# Perguntas:
# Qual é a espécie com o maior numéro de contigs?
taxonomy$species[which.max(taxonomy$num_contigs)]
# Quais famílias apresentam o maior numero de espécies?
names(sort(table(taxonomy$family), decreasing = TRUE)[1])
# No total, a família Flaviviridae apresenta quantas espécies? E quantos contigs?
sum(taxonomy$family == "Flaviviridae")
sum(taxonomy$num_contigs[taxonomy$family == "Flaviviridae"])

# No total, o gênero Orbivírus apresenta quantos contigs?
sum(taxonomy$num_contigs[taxonomy$genus == "Orbivirus"])

# Quais espécies apresentam os 5 maiores números de contigs?
head(sort(taxonomy$num_contigs), 5)
taxonomy$species[order(taxonomy$num_contigs)][1:5]

# E os 5 maiores?
tail(sort(taxonomy$num_contigs), 5)
taxonomy$species[order(taxonomy$num_contigs)][nrow(taxonomy):(nrow(taxonomy)-5)]

# Tarefa de casa: 
# 1) Apliquem esses comandos e treinem a interpretação de seus dados. 
# 2) Elaborem novas perguntas e proponham adaptem os comandos mostrados para respondê-las.

# Sessão 1.2: Pré-processamento 

# Tratamento de dados ausentes
# Determine quais dados estão faltando em seu conjunto de dados e por que estão ausentes. 

sum(is.na(taxonomy$num_contigs)) # Identificar NA na coluna de contigs
sapply(taxonomy, function(x) sum(is.na(x))) # Contar o número de NA por coluna

# Lidar com os NAs
# Remoção
taxonomy2 <- na.omit(taxonomy) # Remove todas as linhas com pelo menos 1 NA
taxonomy2 <- taxonomy[ , colSums(is.na(taxonomy)) == 0] # Remover colunas com pelo menos 1 NA
taxonomy2 <- taxonomy[complete.cases(taxonomy) | rowSums(is.na(taxonomy)) <= 3, ] # Remove as linhas com mais de 3 NA

# Preenchimento de NA
taxonomy[is.na(taxonomy)] <- 0 # Preencher todos os NAs com um valor específico (por exemplo, 0)

# Imputação de dados numéricos
# Use métodos estatísticos, como média, mediana ou moda, para substituir NAs em colunas numéricas
taxonomy$num_contigs[is.na(taxonomy$num_contigs)] <- mean(taxonomy$num_contigs, na.rm = TRUE)

# Modelagem de machine learning: Use modelos de machine learning para prever valores ausentes com base em outras variáveis do conjunto de dados. A abordagem exata depende do modelo escolhido e do contexto.
# Tratar NAs como uma categoria separada
# Análises separadas para observações com e sem NAs, você pode criar subconjuntos de dados.

# Amostragem dos dados
set.seed(123)  # Define uma semente para a reprodutibilidade
amostra <- taxonomy[sample(nrow(taxonomy), size = 10), ] # Amostra aleatória

# Ordenação e filtragem

taxonomy <- taxonomy[order(taxonomy$num_contigs), ]  # Ordena por número de contigs
taxonomy2 <- taxonomy[taxonomy$num_contigs > 1000, ]  # Filtra as espécies com mais de 1000 contigs
taxonomy2 <- taxonomy[order(taxonomy$num_contigs, decreasing = TRUE)[1:10], ] # Filtra apenas as 10 espécies com o maior numero de contigs

# Filtro condicional
# Filtragem condicional com operador &
filtro_condicional <- subset(taxonomy, family == "Flaviviridae" & num_contigs >= 1000)

# Filtragem condicional com operadores OU (||): 
# Filtrar todas as espécies que pertencem à família "Flaviviridae" OU têm mais de 1000 contigs:
filtro_condicional <- subset(taxonomy, family == "Flaviviridae" | num_contigs > 1000)

# Filtragem condicional com operador NÃO
filtro_condicional <- subset(taxonomy, family != "Flaviviridae")

# Filtragem condicional combinando operadores: 
# Filtrar todas as espécies que pertencem à família "Flaviviridae" E NÃO têm menos de 500 contigs
filtro_condicional <- subset(taxonomy, family == "Flaviviridae" | family == "Reoviridae" & !(num_contigs < 500))

# Filtragem por padrões de texto
filtro_padrão_texto <- subset(taxonomy, grepl("orthobunyavirus", species, ignore.case = TRUE))

# Normalização
# Padroniza os dados para que eles estejam na mesma escala, especialmente em experimentos de expressão gênica, onde a normalização é crucial para comparar diferentes amostras.
# Log-transformação: Em alguns casos, aplicar uma transformação logarítmica aos dados pode estabilizar a variância e tornar a distribuição mais próxima de uma distribuição normal.
# Escalonamento: Redimensione os dados para um intervalo específico, como [0, 1], se necessário.

# Min-Max Scaling:
# Esta é uma técnica comum que dimensiona os valores de uma variável para um intervalo específico, como [0, 1].
scaled_data <- (taxonomy$num_contigs - min(taxonomy$num_contigs)) / (max(taxonomy$num_contigs) - min(taxonomy$num_contigs))
taxonomy$scaled_data <- (taxonomy$num_contigs - min(taxonomy$num_contigs)) / (max(taxonomy$num_contigs) - min(taxonomy$num_contigs))

#Z-Score Standardization
# Esta técnica padroniza os valores para que a média seja 0 e o desvio padrão seja 1.
taxonomy$standardized_data <- (taxonomy$num_contigs - mean(taxonomy$num_contigs)) / sd(taxonomy$num_contigs)

# Log Transformation
# Em alguns casos, aplicar o logaritmo aos dados pode ajudar a normalizá-los
taxonomy$log_transformed_data <- log(taxonomy$num_contigs)

## Sessão 2: Sequências de DNA

virus <- readDNAStringSet("virus.fasta.gz")
virus

# Sessão 2.1: Estatística Descritiva e Expĺoração dos dados

head(virus)
length(virus) # Número de sequências
names(virus) # Headers de cada sequência
reverseComplement(virus)
virus[1] # Acessar sequências individualmente
virus[1:3]

# Contagem de bases nucleotídicas
alphabetFrequency(virus)
alphabetFrequency(virus, baseOnly=TRUE, as.prob=TRUE)

letterFrequency(virus, "A", as.prob=TRUE)

# Comprimento das sequências
width(virus)

# Contagem de motivos específicos
motif <- DNAString("ATGC")
vcountPattern(motif, virus)

# Selecionar apenas certos nucleotídeos das sequências
subseq(virus, start=3, end=3) # Seleção apenas do 3 nucleotídeo
subseq(virus, start=1, end=100) # Seleção dos 100 primeiros
subseq(virus, start=10, end=100) # Seleção do nucleotídeo 10 ao 100

# Análise do conteúdo GC
vcountPattern("CG", virus) / width(virus) # Por sequência em %
vcountPattern("GC", virus)

# Sessão 2.2: Pré-processamento 
# Remoção de sequências com % de bases não informativas (!= de A,C,T,G)
df <- cbind(virus@ranges@NAMES, alphabetFrequency(virus, baseOnly=TRUE, as.prob=T))
df <- as.data.frame(df)
colnames(df) <- c("RefSeqIds", "A", "C", "G", "T", "other")
subsequences <- subset(df, df$other <= 0.10) # Manter apenas os abaixo de de 10% 
subset <- virus[subsequences$RefSeqIds]

# Remover sequências a partir de um tamanho mínimo
tamanho_minimo <- 1000  # Substitua pelo tamanho desejado
# Filtrar sequências com tamanho maior ou igual ao tamanho mínimo
virus_filtrado <- virus[nchar(virus) >= tamanho_minimo]

# Remover sequências a partir do valor de GC
# Defina o conteúdo GC mínimo desejado (em porcentagem)
gc_minimo <- 40  # Substitua pelo conteúdo GC desejado

# Função para calcular o conteúdo GC de uma sequência
calculate_gc_content <- function(sequence) {
    gc_count <- sum(unlist(strsplit(as.character(sequence), "")) %in% c("G", "C"))
    sequence_length <- nchar(sequence)
    gc_percentage <- (gc_count / sequence_length) * 100
    return(gc_percentage)
}

# Filtrar sequências com conteúdo GC maior ou igual ao valor mínimo
virus_filtrado <- virus[sapply(virus, calculate_gc_content) >= gc_minimo]

# Visualização
# Histograma de comprimento de sequências
seq_lengths <- lengths(virus)
summary(seq_lengths)
ggplot(data = data.frame(Length = seq_lengths), aes(x = Length)) +
    geom_histogram(binwidth = 100, fill = "blue", color = "black") +
    labs(title = "Histograma de Comprimento de Sequências", x = "Comprimento", y = "Contagem")

# Gráfico de Dispersão de Comprimento vs. Conteúdo GC
seq_lengths <- lengths(virus)
gc_contents <- sapply(virus, calculate_gc_content)  # Use a função calculate_gc_content definida anteriormente
df <- data.frame(Length = seq_lengths, GC_Content = gc_contents)
ggplot(data = df, aes(x = Length, y = GC_Content)) +
    geom_point() +
    labs(title = "Dispersão de Comprimento vs. Conteúdo GC", x = "Comprimento", y = "Conteúdo GC")
