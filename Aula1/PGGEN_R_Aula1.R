### Bem vindos ao Minicurso de Exploração e Visualização de dados biológicos em R
### Msc. Amanda Araújo Serrão de Andrade - Laboratório de Bioinformática (LNCC/UFRJ)
### Duvidas podem ser encaminhadas por email: aandradebio@gmail.com

# Sessão 1: Introdução
# Nesta aula, vamos aprender sobre a importância dos pacotes R, como instalá-los e utilizá-los.
# Também veremos onde encontrar recursos adicionais.

## Pacotes
# Os pacotes são coleções de funções tematicamente relacionadas no R.
# Cada pacote possui documentação detalhada explicando suas funcionalidades.

## Links Úteis:
# CRAN (Comprehensive R Archive Network): Principal repositório de pacotes do R.
# Bioconductor: Pacotes de bioinformática pré-curados.
# GitHub: Repositório de pacotes de desenvolvedores.
# RDocumentation: Documentação detalhada para pacotes.
# rdrr.io: Mais informações sobre pacotes R.

## Recomendações:
# Manter este script .R e os demais arquivos utilizados durante esta prática no mesmo diretório
# Recomendo que refaçam meus comentários na medida que vamos avançando a aula.

# Define o diretório de trabalho como o diretório do script atual
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_dir)

### Sessão 2: Instalação de Pacotes
# Vamos começar instalando os pacotes necessários para o curso.
# Existem várias maneiras de instalar pacotes, incluindo o uso do RStudio ou CRAN.

# Instalando pacotes do tidyverse, que inclui várias bibliotecas úteis.
install.packages("tidyverse")

# Instalando o pacote "data.table" para manipulação eficiente de tabelas grandes.
install.packages("data.table")

# Instalando o pacote "seqinr" especializado em análise de sequências biológicas.
install.packages("seqinr")

# Instalando o pacote "ggplot2" para criação de gráficos personalizáveis.
install.packages("ggplot2")

# Instalando o pacote "RColorBrewer" para paletas de cores pré-definidas em gráficos.
install.packages("RColorBrewer")

# Instalando o pacote "ggtree" para visualização de árvores filogenéticas.
install.packages("ggtree")

# Instalando o pacote "ape" para análise filogenética.
install.packages("ape")

# Instalando o pacote "dplyr" para manipulação eficiente de dados.
install.packages("dplyr")

# Instalando o pacote "tidyr" para organização e pré-processamento de dados.
install.packages("tidyr")

# Instalando o pacote "cowplot" para criação de gráficos complexos.
install.packages("cowplot")

# Instalando o pacote "devtools" para obtenção e desenvolvimento de pacotes R.
install.packages("devtools")

# Instalando o pacote "ggpubr" para ferramentas adicionais de criação de gráficos com ggplot2.
install.packages("ggpubr")

# Instalando o pacote "readxl" para importação de dados de arquivos do Excel.
install.packages("readxl")

# Instalando o pacote "stringr" para manipulação de strings (textos).
install.packages("stringr")

# Instalando o pacote "gridExtra" para layouts personalizados de gráficos múltiplos.
install.packages("gridExtra")

# Sessão 3: Bioconductor
# Bioconductor é uma fonte de pacotes voltados para bioinformática.
# Vamos verificar e instalar o pacote "BiocManager" e, em seguida, instalar o pacote "Biostrings".
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::available() # Listar pacotes disponíveis
BiocManager::install(c("Biostrings"))

# Sessão 4: GitHub
# Às vezes, pacotes ou versões mais recentes podem ser encontrados no GitHub.
# Vamos usar o "devtools" para instalar o pacote "ggplot2" diretamente do GitHub.
devtools::install_github("hadley/ggplot2")

# Sessão 5: Arquivo ZIP
# Você também pode instalar pacotes a partir de arquivos ZIP locais.
install.packages("C:/caminho/para/o/arquivo/zipado/nome-do-pacote.zip", repos = NULL)

# O CRAN, Bioconductor e o Github podem apresentar diferentes versões dos mesmos pacotes. 

# Instalar pacotes se não estiverem instalados
required_packages <- c(
    "tidyverse", "data.table", "seqinr", "ggplot2", "RColorBrewer",
    "ggtree", "ape", "dplyr", "tidyr", "cowplot", "devtools", "ggpubr", "readxl"
)

# Verificar e instalar pacotes ausentes
for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
        install.packages(package, dependencies = TRUE)
    }
}

# Sessão 6: Carregar Pacotes
# Agora que instalamos os pacotes, vamos carregá-los para usá-los.
library(data.table) 
library(ggpubr) 
library(ape) 
library(ggplot2) 
library(seqinr) 
library(RColorBrewer) 
library(dplyr) 
library(tidyr) 
library(cowplot) 
library(readxl) 
library(devtools) 
library(Biostrings)

# Verificar as funções de um pacote

ls("package:tidyr")

# Sessão 7: Menu de Ajuda
# Você pode acessar a documentação de um pacote ou função usando o RStudio ou "?funcao".
help(package = 'seqinr') # Informações sobre o pacote seqinr
?read.fasta # Documentação da função read.fasta

# Sessão 8: Demo Datasets
# Vamos explorar alguns conjuntos de dados de demonstração.
data() # Listar todos os conjuntos de dados disponíveis
data(mtcars)
data(ToothGrowth)
data(iris)
?iris # Obter informações sobre o conjunto de dados "iris"

dim(iris) # Dimensões do data.frame
table(iris$Species) # Contagem de indivíduos por espécie
prop.table(table(iris$Species)) # Proporção em porcentagem

# Sessão 9: Importação de Dados
# Vamos aprender a importar dados de diferentes fontes.

## Tabelas
# A função read.table lê arquivos de tabela para criar um dataframe.
# A função fread do pacote data.table é mais eficiente para arquivos grandes.
# Outras funções como read.csv e read.csv2 são úteis para formatos específicos.

# Exemplos de importação de tabelas
table1 <- read.table(file.choose(), sep = ",", header = TRUE, dec = ".") # Correto
table2 <- read.table(file = "tax.csv", header = TRUE, sep = ",", dec = ".") # Correto
table3 <- read.table(file = "tax.tab", header = TRUE, sep = "\t", dec = ".") # Correto

table4 <- read.table(file = "tax.txt", sep = ",") # Incorreto (separador diferente)
table5 <- read.table(file = "tax.txt", sep = "\t") # Incorreto (separador diferente)
table6 <- read.table(file = "tax.txt", sep = "|", header = TRUE) # Correto

# Exemplos de importação de tabelas usando fread
# A função fread é parte do pacote data.table e é especialmente eficiente para importar grandes conjuntos de dados.

# Importação de um arquivo CSV com fread
table1 <- fread("tax.csv", header = TRUE) 
# Comentário: Aqui, estamos importando um arquivo CSV usando fread. O argumento header = TRUE indica que a primeira linha contém nomes de coluna.

# Importação de um arquivo TSV (separado por tab) com fread
table2 <- fread("tax.tab", header = TRUE, sep = "\t")
# Comentário: Neste caso, estamos importando um arquivo TSV (separado por tab) com fread. Usamos o argumento sep = "\t" para especificar o separador.

# Importação de um arquivo TXT com um separador personalizado usando fread
table3 <- fread("tax.txt", header = TRUE, sep = "|")
# Comentário: Aqui, estamos importando um arquivo TXT com um separador personalizado ("|"). Novamente, indicamos que a primeira linha contém nomes de coluna.

# Diferença entre fread e read.table:
# A principal diferença entre fread e read.table está na eficiência de importação de dados. fread é notavelmente mais rápido e eficiente, o que o torna ideal para grandes conjuntos de dados. Além disso, fread é mais flexível na detecção automática do formato do arquivo, como delimitadores e formatos de números.

# Use fread quando você precisa de alta eficiência, especialmente com grandes conjuntos de dados e quando deseja que o R detecte automaticamente o formato do arquivo.
# Use read.table quando estiver trabalhando com arquivos de texto simples e precisar de flexibilidade adicional na especificação de formatos.
# Use read.csv quando estiver trabalhando especificamente com arquivos CSV e quiser uma configuração padrão mais orientada para CSV.

# Importação de arquivos Excel
table7 <- read_excel("tax.xlsx")
table8 <- read_excel("tax.xlsx", sheet = "Sheet1") # Por nome
table9 <- read_excel("tax.xlsx", sheet = 1) # Pelo número

# Importação de dados de URLs
table_url1 <- read.csv('https://raw.githubusercontent.com/Statology/Miscellaneous/main/basketball_data.csv')
table_url2 <- fread('http://www.sthda.com/upload/boxplot_format.txt')

# Sequências de nucleotídeos/aminoácidos
# Importando a partir de um arquivo fasta

virusSequences <- readDNAStringSet(file="virus.fasta.gz") 
virusSequences <- readAAStringSet(file="virus.faa.gz") # Aminoácidos

# Obtenção de sequências por meio do R
virusIDS <- c("MH734972.1","NC_025353.1","AY898809.1","KJ399977.1","NC_043636.1","MG029271.1",
              "NC_055531.1","MW809639.1","KU173874.1","HQ630929.1")
virusSequences <- read.GenBank(virusIDS, seq.names = virusIDS, as.character = TRUE)
write.dna(virusSequences, 'virusSequences.fasta', format = "fasta")

### Exportação de dados
# A função write.table pode ser utilizada para salvar um dataframe/matrix em arquivo
# write.csv exporta o ponto "." como ponto decimal e a vírgula como separador
# write.csv2 exporta a vírgula como ponto decimal e o ; como separador

# Parâmetros: 
# x = o dataframe/matrix que vai ser salvo
# file = o arquivo resultante
# append = TRUE/FALSE. Adicionar a um arquivo existente (em vez de o substituir)
# sep = O caracter separador de campos. Ex: "\t"
# dec = O caracter utilizado para os pontos decimais
# row.names = TRUE/FALSE
# col.names = TRUE/FALSE

data(mtcars)

# Tabelas

write.table(iris, file="iris1.csv", sep=",", row.names=T, col.names = T)
write.table(iris, file="iris2.csv", sep=",", row.names=F, col.names = T)
write.table(iris, file="iris3.csv", sep=",", row.names=F, col.names = F)
write.table(iris, file="iris4.tab", sep="\t", row.names=T, col.names = T)
write.table(iris, file="iris5.tab", sep="\t", row.names=T, col.names = T, dec = ",")
write.table(mtcars, file="iris5.tab", sep="\t", row.names=F, col.names = F, dec = ",", append=T)

# Formato Excel
data(mtcars)
write.xlsx(iris, file="iris.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, 
           append = FALSE)
write.xlsx2(mtcars, "iris.xlsx", sheetName = "mtcars", col.names = TRUE, row.names = TRUE, 
            append = TRUE)

# Capturar a saída de um comando
output <- prop.table(table(iris$Species))
capture.output(output, file="output.txt")
output <- prop.table(table(mtcars$mpg))
capture.output(output, file="output.txt",append = T)

capture.output(prop.table(table(iris$Species)), file="output.txt",append = T)

# RDS
# Salva um objeto para um arquivo .rds
saveRDS(mtcars, "mtcars.rds")
# Importa o objeto novamente
mtcars <- readRDS("mtcars.rds")
data(iris)

# RData
# Não perder o processo. Pode ser carregado em diferentes computadores
# Salva múltiplos objetos em um arquivo .RData
save(mtcars, file = "Aula1.RData")
# Save multiple objects
save(mtcars, iris, file = "Aula3.RData")
# To load the data again
load("Aula1.RData")

save.image(file = "Aula1.RData")
load("Aula1.RData")

# Sessão 11: Boas Práticas
# Algumas dicas úteis para nomear variáveis e dados.
# Evitar espaços entre os nomes. Ex: Zika virus -> Zika_virus
# Aos dados faltantes devem ser atribuídos 0 ou NA (not available)
# Evitar símbolos ?, $, *, +, #, (, ), -, /, }, {, |, >, < 
# Evitar iniciar o nome das colunas com número. Ex: 123virus -> virus123 ou v123
# Os nomes das colunas devem ser únicos

# Limpeza de resultados intermediários, se necessário
rm(table1, table2, table3, table4, table5, table6, table7, table8, table9, table_url1, table_url2, virusSequences, virusIDS)

# Encerramento da aula
cat("Aula 1 concluída. Lembre-se de praticar e explorar mais!\n")