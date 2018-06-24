#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
# UNIVERSIDADE FEDERAL DO PARANA                                                                     #
# PPG-ENTOMOLOGIA                                                                                    #
#                                                                                                    #
# Maio de 2017                                                                                       #
# Disciplina: "Desenho amostral e ferramentas estatisticas para o estudo de comunidades de animais"  #
# Professor: Sebastian Sendoya                                                                       #
#                                                                                                    #
# Aula: Analises de ordenaacao                                                                       #
# Professora: Luciana de Campos Franci (lucianafranci@gmail.com)                                     #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

#Antes de tudo, defina seu local de trabalho (working directory)


# Escalonamento Multidimensional Nao Metrico (NMDS) -----------------------

library(vegan)

set.seed(2) #indica que a semente da aleatorizacao
com_matriz <- matrix(sample(1:100,300, replace = T), nrow = 10, dimnames = list(paste("com_",1:10,sep=""),
                                                                                paste("sp",1:30,sep="")))

# NMDS1 -----------------------------------------------------------------------

#Mostrara cada iteracao do NMDS ate que uma solucao seja atingida (ou seja, o estresse foi minimizado apos 
#um numero de reconfiguracoes dos pontos em duas dimensoes)
com_NMDS <- metaMDS(com_matriz, k = 2, distance = "bray") #k indica o numero de dimensoes reduzidas
com_NMDS

#A funcao metaMDS automaticamente transforma os dados para raiz quadrada e calcula as distancias de Bray-Curtis
#para a matriz de comunidade por local

#Vamos testar qual valor de k (dimensoes) gera menor valor de stress
for(i in 1:5){
  print(metaMDS(com_matriz, distance = "bray", k = i, trace = F)$stress) #vai retornar os valores de stress
}

#Refazendo a NMDS usando k = 3
com_NMDS <- metaMDS(com_matriz, k = 3, distance = "bray") #k indica o numero de dimensoes reduzidas
com_NMDS

#extraindo os escores
scores(com_NMDS)

# NMDS2 -----------------------------------------------------------------------


#Analisando o gráfico de Shepard (grafico de estresse), o qual mostra a dispersao em torno da regressao entre
#a distancia inter-pontos na configuracao final (ou seja, a distancia entre cada par de comunidades) e 
#a dissimilaridade original
stressplot(com_NMDS)

#Grande dispersao em torno da linha sugere que as dissimilaridades originais nao estao bem preservadas
#no numero reduzido de dimensoes

# NMDS3 -----------------------------------------------------------------------


#Plotando a NMDS. O gráfico mostra tanto as comunidades ('locais', circulos abertos) e especies (cruzes vermelhas)
#mas sem identificar cada ponto
plot(com_NMDS, cex  = 1.5, display = "sites")

#Usando as funcoes ordiplot e orditorp  para adicionar rotulos ao grafico
ordiplot(com_NMDS, type = "n")
orditorp(com_NMDS, display = "species", col = "red", 
         air = 0.01) #air indica o espaco permitido entre os rotulos, se for menor que 1 permite sobreposicao
orditorp(com_NMDS, display = "sites", cex = 1.25, air = 0.01)

# NMDS4 -----------------------------------------------------------------------


#Dicas
#1.Supondo que as comunidades de 1 a 5 passaram por um tratamento e de 6-10 por um tratamento diferente. Pode-se
#desenhar poligonos ou elipse conectando as comunidades. Esse e um modo intuitivo de mostrar como as comunidades e especies
#se agrupam baseando-se em tratamentos
trat <- c(rep("Tratamento1",5), rep("Tratamento2",5)) #criando um objeto com os tratamentos
ordiplot(com_NMDS, type = "n") 
orditorp(com_NMDS, display = "sites", 
         col = c(rep("darkgreen",5), rep("blue",5)), air = 0.01, cex = 1.25) #colocando as comunidades
ordihull(com_NMDS, groups= trat, draw = "polygon", 
         col = c("lightgreen","lightblue"), label = F) #para desenhar elipses use ordiellipse
orditorp(com_NMDS, display = "species", 
         col = "red", air = 0.01) #colocando os pontos para as especies


# NMDS5 -----------------------------------------------------------------------


#2.Supondo que o tratamento e um gradiente como a elevacao de um terreno ou temperatura, 
#podemos plotar as linhas de contorno
cota <- runif(10,0.5,1.5) #criando um objeto com as cotas altimetricas
ordisurf(com_NMDS, cota, main = "", col = "black")
orditorp(com_NMDS, display = "species", col = "red", air = 0.1, cex = 1)
orditorp(com_NMDS, display = "sites", col = c(rep("green",5), rep("blue",5)), air = 0.01, cex = 1.25)


# Analise de componentes principais (PCA) --------------------------------------------------------------------

#Abrindo os dados 'iris'
data(iris)
head(iris)

#Transformando as colunas de 1 a 4 usando log para normalizar os dados e salvando num objeto
log_ir <- log(iris[, 1:4])

#Salvando as especies em um objeto
ir_species <- as.data.frame(iris[, 5])
colnames(ir_species)[1] <- "specie"

# PCA1 --------------------------------------------------------------------

#Abrindo o pacote 'vegan'
library(vegan)

#Fazendo a analise de PCA e salvando em um objeto
iris_pca <- rda(log_ir, scale = T) #scale = T indica que os dados serao padronizados
summary(iris_pca) #pedindo os resultados
#Saida:
#inertia = a soma da  variancia de todas as variaveis
#eigenvalues = a quantidade de variancia explicada por cada componente principal (medida de importancia de cada PC)
#species scores = autovetores
#sites scores = coordenadas dos locais

# PCA2 --------------------------------------------------------------------

#Para sabermos quanto cada eixo explica, precisamos dividir o autovalor de cada eixo pela inercia total
sum(apply(scale(log_ir), 2, var)) #tem que ser igual ao valor de inercia

eixo1_iris <- 2.9325/4; eixo1_iris
eixo2_iris <- 0.907/4; eixo2_iris
eixo3_iris <- 0.133/4; eixo3_iris
eixo4_iris <- 0.0275/4; eixo4_iris

#Plotando os autovalores de cada PC
screeplot(iris_pca, type = "l")

#Sabemos que a PCA foi a analise adequada quando no maximo ate os tres primeiros eixos explicam ate 70% da variacao dos dados

# PCA3 --------------------------------------------------------------------

#Plotando os eixos num biplot
colvec <- c("orange", "blue", "pink") #criando um vetor com as cores a serem usadas para as especies no grafico
levels(ir_species$specie) #as cores serao colocadas de acordo com a ordem dos niveis do vetor onde elas estao guardadas

#eixos 1 e 2
#scaling indica a escala dos dados. 
#scaling = 1 - Boa opcao quando o foco for interpretar as relacoes entre locais. O comprimento do vetor indica a importancia da variavel.
#scaling = 2 - Boa opcao quando o foco for interpretar as relacoes entre especies (ou as variaveis explanatorias como nesse exemplo)#scaling = 3 - escala simetrica tanto para especie quanto para locais
#scaling = 3 - Escalona locais e especies

biplot(iris_pca, scaling = 2, type = c("points", "points"), col = "black", choices = c(1,2), xlab = "PC1 (73.3%)",
       ylab = "PC2 (22.6%)") #plotando o biplot - choices indica quais eixos serao plotados, no caso, eixos 1 e 2.
with(ir_species, points(iris_pca, scaling = 2, display = "sites", col = colvec[specie],
                        pch = 21, bg = colvec[specie], choices = c(1,2)))
text(iris_pca, display = "species", cex = 0.9, col = "black", choices = c(1,2))
with(ir_species, legend("bottomright", legend = c("setosa", "versicolor", "virginica"), bty = "n", 
                        col = colvec, pch = 21, pt.bg = colvec))

#Colocando o circulo de contribuicao de equilibro - serve para indicar quais especies (ou medidas) sao as mais representativas, ou melhor, que mais contribuiram para os eixos
#Os vetores que ultrapassam o circulo sao os que mais contribuem, os que nao ultrapassam podem ser considerados redundantes aos que contribuem
d <- 2 #numero de dimensoes do grafico, no exemplo 2 dimensoes (eixos 1 e 2)
p <- 4 #numero de especies (no exemplo, sao 4 medidas)
cec <- sqrt(d/p) #calculando a raiz quadrada da razao entre dimensao e numero de especies/medidas
circ <- seq(0,2*pi,length=100)
coords <- t(rbind(sin(circ)*cec, cos(circ)*cec))
lines(coords, col="red")

# PCA4 --------------------------------------------------------------------


#eixos 2 e 4
biplot(iris_pca, scaling = 3, type = c("points", "points"), col = "black", choices = c(2,4), xlab = "PC2 (22%)",
       ylab = "PC4 (0.6%)")
with(ir_species, points(iris_pca, scaling = 3, display = "sites", col = colvec[specie],
                        pch = 21, bg = colvec[specie], choices = c(2,4)))
text(iris_pca, display = "species", cex = 0.9, col = "black", choices = c(2,4))
with(ir_species, legend("bottomleft", legend = c("setosa", "versicolor", "virginica"), bty = "n",
                        col = colvec, pch = 21, pt.bg = colvec))


# PCA5 --------------------------------------------------------------------


#Extraindo os eixos (scores) para usar em analises estatisticas como variaveis explanatorias
#Quantos eixos usar nas analises? Preferencialmente usar os eixos que, somados, expliquem ate 80%
#Nesse exemplo, o primeiro eixo explica 73%, e o segundo explica 22%, entao primeiro eixo ja e suficiente
ir_eixos <- round(as.data.frame(scores(iris_pca, choices = 1, display = "sites", scaling = 2)), digits = 6)
head(ir_eixos)#mostra apenas as 5 primeiras linhas do objeto
write.csv(ir_eixos, "eixos_iris.csv", sep = ";", row.names = F)

# Analise de Correspondencia Canonica (CCA) --------------------------------------------------------------------

library(vegan)

#Abrindo os dados das especies
data(dune)
str(dune)

dune_relativa <- decostand(dune, method = "total") #calculando a abundancia relativa de cada espécie por local
dune_relativa[1:5, 1:5] #vendo as 5 primeiras linhas e 5 primeiras colunas para ver se transformou corretamente

#Abrindos os dados ambientais
data(dune.env)
str(dune.env)


# CCA1 --------------------------------------------------------------------

#Rodando a CCA

dune_cca <- cca(dune_relativa ~ Moisture + Use + A1, data = dune.env, scale = TRUE)
summary(dune_cca)
print(dune_cca)


# CCA2 --------------------------------------------------------------------

#Testando se as variaveis sao redundantes usando o fator de inflacao da variancia (VIF)
#Valores maiores que 10 indicam redundancia

vif.cca(dune_cca)

# CCA3 --------------------------------------------------------------------

#Plotando a CCA em um triplot (plotando especies, locais, variaveis explanatorias)

plot(dune_cca, scaling = 2, choices = c(1,2), type = "n", xlab = "CCA1 (20.6%)", ylab = "CCA2 (11.5%)") #scaling = 2 porque nos interessa as especies
text(dune_cca, display = "cn")
points(dune_cca, pch = 21, col = "red", bg = "orange", cex = 1.5)
#para colocar os rotulos dos lugares e so usar text(dune_cca, display = "wa")
text(dune_cca, "species", col = "blue", cex = 1)

# Análise de redundancia (RDA) --------------------------------------------------------------------

library(vegan)
#Abrindo os dados das especies
data(dune)
str(dune)

#Abrindos os dados ambientais
data(dune.env)
str(dune.env)

# RDA1 --------------------------------------------------------------------

#Rodando a RDA

dune_rda <- rda(dune ~ Moisture + Use + A1, data = dune.env, scale = TRUE)
print(dune_rda)


# RDA2 --------------------------------------------------------------------

#Testando se as variaveis sao redundantes usando o fator de inflacao da variancia (VIF)
#Valores maiores que 10 indicam redundancia

vif.cca(dune_rda)


# RDA3 --------------------------------------------------------------------

#Plotando a RDA em um triplot (plotando especies, locais, variaveis explanatorias)

plot(dune_rda, scaling = 2, choices = c(1,2), type = "n")
text(dune_rda, display = "cn")
points(dune_rda, pch = 21, col = "red", bg = "orange", cex = 1.5)
#para colocar os rotulos dos lugares é só usar text(dune_cca, display = "wa")
text(dune_rda, "species", col = "blue", cex = 1)

# Comparando a CCA e a RDA ------------------------------------------------

#Testando a significancia dos eixos na RDA e na CCA para saer qual analise foi mais adequada

anova(dune_cca, by = "axis", step = 999)
anova(dune_rda, by = "axis", step = 999)
