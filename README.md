# Análisis de trazos e identificación de potenciales eventos vicariantes en el ambiente R


This repository contains source code and documentation for the project "panBioGeo: panbiogeographic tools (and other things)."

Please note that although comments & pull requests are welcome, those not involving an actual software bug (which we always try to fix as promptly as we can) may not be addressed in a timely manner.

## Authors: Augusto Ferrari, José Ricardo I. Ribeiro & Diego Janisch Alvares


UNIDAD DE ACCION:

FASE PREPARATÓRIA

Função que remove todos os objetos do ambiente
> rm(list = ls())

INSTALANDO PACOTES NECESSÁRIOS: necessário apenas na primeira vez
> install.packages('maptools');install.packages('geiger')
> install.packages('phytools');install.packages('devtools') 
> install.packages('viridis');install.packages('mapdata')
> install.packages('rgdal');install.packages("letsR")
> install.packages('fossil');install.packages('vegan')
> install.packages("spdep");install.packages('dismo');install.packages("TreeSearch")
> install.packages('spatstat');install.packages('ctv');install.packages('Xming')
> install.packages("stringr"); install.packages('stringi')


CARREGANDO PACOTES PARA A PRÁTICA: necessário sempre que rodar a prática
> library(phytools)
> library(devtools)
> library(rgdal)
> library(maptools)
> library(geiger)
> library(raster)

## INDICANDO O DIRETÓRIO DE TRABALHO ##
recomendamos que todos os arquivos da prática sejam descarregados no mesmo diretório
> setwd(" <endereço> ") # direcionando o R para esse diretório

***
# BLOQUE 1

CARREGANDO O ARQUIVO COM OS REGISTROS GEOREFERENCIADOS DAS ESPÉCIES
> lycipta.coords <- read.table('Lycipta.csv', sep = ',', dec = '.', header = TRUE)

visualização da estrutura dos dados:
OBS: Para utilização de banco de dados próprio, organizar o arquivo com o seguinte cabeçalho e coordenadas em decimais: spp, Long, Lat
> head(lycipta.coords, 10)

visualizacao espacial dos dados (mapa)
CARREGANDO SHAPE DA AMÉRICA DO SUL COM CLASSIFICAÇÃO BIOGEOGRÁFICA
> asul <- rgdal::readOGR(dsn = 'Lowenberg_Neto_2014.shp')


PLOTANDO OS PONTOS
> plot(asul, axes = TRUE, col = "gray")

> points(lycipta.coords[,2:3], pch = 16, col = "red", cex = 1)

CARREGANDO A FILOGENIA CORRESPONDENTE AOS DADOS DE DISTRIBUIÇÃO

> lycipta.tree <- read.tree(file = 'Lycipta.tre')

VISUALIZAÇÃO DA FILOGENIA:
> lycipta.tree <- drop.tip(lycipta.tree, "ROOT")
> plotTree(lycipta.tree, node.numbers = TRUE, ftype = "i")
> tiplabels()

***
# BLOQUE 2

carregar a função que gera os traços individuais
> source('node_terminal.R')

ATENÇÃO: a função abaixo irá gerar os MST's para todas as espécies presentes no arquivo inicial.
Neste ponto definimos a resolução espacial (argumento "resol") que será mantida nas análises subsequentes (PAE_PCE)


> mst_all_taxa <- terminal_node(coordin = lycipta.coords, shape_file = asul,
                    sobrepo = FALSE, caption = TRUE, resol = c(10, 10), tree = lycipta.tree, xmin = min(lycipta.coords$Long) - 5,
                    xmax = max(lycipta.coords$Long), ymin = min(lycipta.coords$Lat) - 7, ymax = max(lycipta.coords$Lat) + 7)


***
# BLOQUE 3

Exemplo 1
Permite-se selecionar apenas um táxon/nó terminal: aqui é possível alterar tamanho de quadrícula como opção de visualização de MST's em outras escalas espaciais (i.e. resol = c(5, 5))

> mst_single_taxa1 <- terminal_node(taxon = c('E_L_imitator'), coordin = lycipta.coords, shape_file = asul, sobrepo = FALSE,
                        caption = TRUE, resol = c(2, 2), tree = lycipta.tree, xmin = min(lycipta.coords$Long) - 5,
                        xmax = max(lycipta.coords$Long), ymin = min(lycipta.coords$Lat) - 8, ymax = max(lycipta.coords$Lat) + 7)

> mst_single_taxa2 <- terminal_node(taxon = c('E_L_illotus'), coordin = lycipta.coords, shape_file = asul, sobrepo = FALSE,
                        caption = TRUE, resol = c(2, 2), tree = lycipta.tree, xmin = min(lycipta.coords$Long) - 5,
                        xmax = max(lycipta.coords$Long), ymin = min(lycipta.coords$Lat) - 7, ymax = max(lycipta.coords$Lat) + 7)

OBS.: os MST's individuais podem ser visualizados em conjunto, adicionando os nomes dos terminais ao argumento "taxon=", como exemplificado abaixo:

> mst_multi_taxa <- terminal_node(taxon = c('E_L_illotus','E_L_imitator'), coordin = lycipta.coords, shape_file = asul,
                      sobrepo = FALSE, caption = TRUE, resol = c(2, 2), tree = lycipta.tree, xmin = min(lycipta.coords$Long) - 5,
                      xmax = max(lycipta.coords$Long), ymin = min(lycipta.coords$Lat) - 7, ymax = max(lycipta.coords$Lat) + 7)



Exemplo 2
Permite-se selecionar apenas um nó: essa função irá sobrepor MST's das espécies que compoem um determinado nó
ATENÇÃO: aqui visualizaremos somente os nós internos; para nós terminais, utilizar o argumento "taxon" como acima no exemplo 1;
aqui tambem é possível alterar tamanho de quadrícula como opção de visualização de MST's em outras escalas espaciais (i.e. resol = c(5, 5))

nó 22:
> mst_single_node22 <- terminal_node(coordin = lycipta.coords, tree = lycipta.tree, shape_file = asul,
                         caption = TRUE, resol = c(5, 5), nodes = 22, seephylog = FALSE, sobrepo = FALSE,
                         xmin = min(lycipta.coords$Long) - 5, xmax = max(lycipta.coords$Long),
                         ymin = min(lycipta.coords$Lat) - 7, ymax = max(lycipta.coords$Lat) + 7)
                         
nó 18:
> mst_single_node18 <- terminal_node(coordin = lycipta.coords, tree = lycipta.tree, shape_file = asul,
                         caption = TRUE, resol = c(5, 5), nodes = 18, seephylog = FALSE, sobrepo = FALSE,
                         xmin = min(lycipta.coords$Long) - 5, xmax = max(lycipta.coords$Long),
                         ymin = min(lycipta.coords$Lat) - 7, ymax = max(lycipta.coords$Lat) + 7)
                         
nó 20:
> mst_single_node20 <- terminal_node(coordin = lycipta.coords, tree = lycipta.tree, shape_file = asul,
                         caption = TRUE, resol = c(5, 5), nodes = 20, seephylog = FALSE, sobrepo = FALSE,
                         xmin = min(lycipta.coords$Long) - 5, xmax = max(lycipta.coords$Long),
                         ymin = min(lycipta.coords$Lat) - 7, ymax = max(lycipta.coords$Lat) + 7)


***
# BLOQUE 4

realização da análise de PAE-PCE dentro do ambiente R

> source('pae_pce_pan.R')

-para verificarmos as sinapomorfias (homoplásicas e não homoplásicas): IR > 0 (argumento "oneGrid = FALSE")
> pae1 <- pae_pce(mat = mst_all_taxa, k = 1000, maxit = 100, oneGrid = FALSE, shapeFile = asul, resolut = c(10, 10),
                gridView = TRUE, labelGrid = TRUE, legendSpecies = TRUE, N = 5)

ao visualizar o objeto que contém o resultado, podemos relacionar quais espécies estão associadas a quais quadrículas do grid

> pae1

-para verificarmos as autapomorfias: IR = 0 e IC > 0 (argumento "oneGrid = TRUE" e nVoltas = número de iterações ocorridas durante a execução da pae_pce com o argumento "oneGrid = FALSE")
> pae2 <- pae_pce(mat = mst_all_taxa, k = 1000, maxit = 100, oneGrid = TRUE, shapeFile = asul, resolut = c(10, 10),
                gridView = TRUE, labelGrid = TRUE, legendSpecies = TRUE, nVoltas = 1)

-é possivel inserir alguns taxons que sustentam os traços generalizados na solução final...

> ex1 <- rgdal::readOGR(dsn = 'out/mst_E_L_longicornis.shp', verbose = FALSE)
> plot(ex1, cex = 1.1, pch = 21, col = 'blue', add = TRUE, lwd = 3, lty = 1)

> ex2 <- rgdal::readOGR(dsn = 'out/mst_E_L_illotus.shp', verbose = FALSE)
> plot(ex2, cex = 1.1, pch = 21, col = 'black', add = T, lwd = 3, lty = 1)

> ex3 <- rgdal::readOGR(dsn = 'out/mst_E_L_picticornis.shp', verbose = FALSE)
> plot(ex3, cex = 1.1, pch = 21, col = 'green', add = T, lwd = 3, lty = 1)

> ex4 <- rgdal::readOGR(dsn = 'out/mst_E_L_cribrarius.shp', verbose = FALSE)
> plot(ex4, cex = 1.1, pch = 21, col = 'yellow', add = T, lwd = 3, lty = 1)


***
## Esse passo é um passo alternativo ao passo 4:

# BLOQUE 5

realização da análise de PAE-PCE com o uso do TNT

carregar a função que gera a matriz de presença e ausência
> source('tnt_matrix.R')

a função gera a matriz de presença e ausência para rodar no TNT (o arquivo resultante ficará disponível no diretório de trabalho)
> tnt_matrix(mst_all_taxa)

## Obs: a matriz resultante terá como base o tamanho de quadrícula definido no BLOQUE 2


***
# BLOQUE 6

COMPARAÇÕES DE GRUPOS-IRMÃOS (TABELA 1)

-COMPARAÇÃO DO ENTRE AS ESPÉCIES DO GRUPO-IRMÃO I
> grupoirmao1 <- terminal_node(taxon = c('E_L_longicornis'), coordin = lycipta.coords, tree = lycipta.tree,
                   shape_file = asul, caption = TRUE, resol = c(5, 5), nodes = 14, seephylog = FALSE, sobrepo = FALSE,
                   xmin = min(lycipta.coords$Long) - 5, xmax = max(lycipta.coords$Long), ymin = min(lycipta.coords$Lat) - 7,
                   ymax = max(lycipta.coords$Lat) + 7)
  
-COMPARAÇÃO ENTRE AS ESPÉCIES DO GRUPO-IRMÃO II
> grupoirmao2 <- terminal_node(coordin = lycipta.coords, tree = lycipta.tree, shape_file = asul, caption = TRUE,
                   resol = c(5, 5), nodes = c(15, 23), seephylog = FALSE, sobrepo = FALSE, xmin = min(lycipta.coords$Long) - 5,
                   xmax = max(lycipta.coords$Long), ymin = min(lycipta.coords$Lat) - 7, ymax = max(lycipta.coords$Lat) + 7)

-COMPARAÇÃO ENTRE AS ESPÉCIES DO GRUPO-IRMÃO III
> grupoirmao3 <- terminal_node(coordin = lycipta.coords, tree = lycipta.tree, shape_file = asul, caption = TRUE,
                   resol = c(5, 5), nodes = c(16, 17), seephylog = FALSE, sobrepo = FALSE, xmin = min(lycipta.coords$Long) - 5,
                   xmax = max(lycipta.coords$Long), ymin = min(lycipta.coords$Lat) - 7, ymax = max(lycipta.coords$Lat) + 7)

-COMPARAÇÃO ENTRE AS ESPÉCIES DO GRUPO-IRMÃO IV
> grupoirmao4 <- terminal_node(taxon = c('E_L_triangulator','E_L_picticornis'), coordin = lycipta.coords, tree = lycipta.tree,
                   shape_file = asul, caption = TRUE, resol = c(5, 5), seephylog = FALSE, sobrepo = FALSE,
                   xmin =  min(lycipta.coords$Long) - 5, xmax = max(lycipta.coords$Long), ymin = min(lycipta.coords$Lat) - 7,
                   ymax = max(lycipta.coords$Lat) + 7)


***
# BLOQUE 7

para sabermos se existe ou não vicariância entre os MSTs dos nós:

> setwd(' <ENDEREÇO> ')

> source('<ENDEREÇO INDICANDO ONDE ESTÁ A FUNÇÃO is_allopatry.R> ')

-a função produz um arquivo .txt no diretório /out, chamado "vicarius.txt": devemos SEMPRE usar o resultado obtido com a função node_terminal de todas as espécies (os MSTs de todas as espécies da árvore)!

> is_allopatry(results = mst_all_taxa, phylogeny = lycipta.tree, coordin = lycipta.coords)
