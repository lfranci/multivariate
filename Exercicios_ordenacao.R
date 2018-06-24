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

#Exercicios

#Escalonamento Multidimensional Nao Metrico (NMDS)
# 1 -----------------------------------------------------------------------

#Faca a analise NMDS para a matriz de comunidades do arquivo 'dune'
data(dune)
#Interprete se o estresse foi alto ou nao

#  ------------------------------------------------------------------------

library(vegan)

for(i in 1:5){
  print(metaMDS(dune, distance = "bray", k = i, trace = F)$stress) #vai retornar os valores de stress
}

#5 dimensoes resultou em menor stress
dune_nmds <- metaMDS(dune, k = 5)
dune_nmds


# 2 -----------------------------------------------------------------------


#Faca o grafico de Shepard (estresse) e interprete o resultado

#  ------------------------------------------------------------------------


stressplot(dune_nmds)


# 3 -----------------------------------------------------------------------


#Se o estresse tiver sido baixo e se a dispersao dos pontos em torno do estresse plot nao for grande, plote a NMDS
#usando os nomes das especies e os locais
#Pelo contrario, refaca a NMDA usando um numero de dimensoes maior

#  ------------------------------------------------------------------------


plot(dune_nmds, cex  = 1.5)

ordiplot(dune_nmds, type = "n")
orditorp(dune_nmds, display = "species", col = "red", cex = 1, air = 0.01)
orditorp(dune_nmds, display = "sites", cex = 1.25, air = 0.01)


# 4 -----------------------------------------------------------------------


#Suponha que as regioes de 1 a 7 sao areas com amplitude de temperatura entre 0 e 7 graus
#enquanto as areas de 8-15 tem amplitude entre 5 e 18 graus, e as areas de 15-20, amplitude entre 17 e 32 graus.
#Faca uma grafico de NMDS indicando essa informacao

#  ------------------------------------------------------------------------


temp <- c(runif(7,0,7), runif(8,5,18), runif(5,17,32)) #criando um objeto com as cotas altimetricas

ordisurf(dune_nmds, temp, main = "", col = "gray", type = "n")
orditorp(dune_nmds, display = "species", col = "blue", air = 0.1, cex = 1)
orditorp(dune_nmds, display = "sites", col = "red", air = 0.01, cex = 1.25)


# 5 -----------------------------------------------------------------------


#Suponha que as areas de 1 a 7 estao presentes em um area A, areas de 8-15 em uma area B e areas de 15-20 em uma area C
#Plot um poligono indicando essas areas no grafico do NMDS

#  ------------------------------------------------------------------------


areas <- c(rep("Area_A",7), rep("Area_B",8), rep("Area_C",5)) #criando um objeto com as areas
ordiplot(dune_nmds, type = "n") 
orditorp(dune_nmds, display = "sites", col = c(rep("darkgreen",7), rep("blue",8), rep("red",5)), 
         air = 0.01, cex = 1.25) #colocando as comunidades
ordihull(dune_nmds, groups= areas, draw = "polygon", 
         col = c("lightgreen","lightblue", "lightpink"), label = F) #para desenhar elipses use ordiellipse
orditorp(dune_nmds, display = "species", col = "Black", air = 0.01) #colocando os pontos para as especies


# Analise de Componentes Principais (PCA) -----------------------------------------------------------------------

#Faca uma PCA dos dados 'dune'

data(dune)
library(vegan)
#Transformacao de Hellinger - usada para dados de abundancia, essa transformacao da peso baixo para variaveis com pouca
#contagem ou muitos zeros
dune_hell <- decostand(dune, "hellinger")


#  ------------------------------------------------------------------------

dune_pca <- rda(dune_hell, scale = T)
summary(dune_pca)

# 1 -----------------------------------------------------------------------

#Faca um biplot dos eixos 1 e 2 com scaling = 2, inclua o circulo de contruicao de equilibrio e compare com o resultado da NMDS do arquivo dune

#  ------------------------------------------------------------------------

biplot(dune_pca, scaling = 2, xlab = "PC1 (23,8%)", ylab = "PC2 (17,3%)", choices = c(1,2))
d <- 2 
p <- 30 #30 variaveis 
cec <- sqrt(d/p) #calculando a raiz quadrada da razao entre dimensao e numero de especies
circ <- seq(0,2*pi,length=100)
coords <- t(rbind(sin(circ)*cec, cos(circ)*cec))
lines(coords, col="blue")

# 2 -----------------------------------------------------------------------

#extraia os scores dos eixos e salve em um objeto

#  ------------------------------------------------------------------------

dune_eixos <- round(as.data.frame(scores(dune_pca, choices = 1:4, display = "sites", scaling = 2)), digits = 6)
head(dune_eixos)

# Analise de Correspondencia Canonica (CCA) -------------------------------
#Faca uma CCA dos dados varespec (valor de cobertura de especies) e varechem (variaveis ambientais)
data(varespec)
varespec_hel <- decostand(varespec, "hellinger")
data(varechem)

#  ------------------------------------------------------------------------
vare_cca <- cca(varespec_hel ~ N + P + Ca + Mg + S + Al + Fe + Mn + Zn + Mo + Baresoil + Humdepth + pH, data = varechem, scale = T)
summary(vare_cca)

# 1 -----------------------------------------------------------------------

#Usando o Valor de Inflacao da Variancia (VIF), verifique se ha variaveis redundantes. Se houver, retire-as e refaca o modelo

#  ------------------------------------------------------------------------

vif.cca(vare_cca) #S e Al sao redundantes, vamos retirar a com maior valor de VIF, ou seja, o Al
#Refazendo o modelo sem Al
vare_cca <- cca(varespec_hel ~ N + P + Ca + Mg + S + Fe + Mn + Zn + Mo + Baresoil + Humdepth + pH, data = varechem, scale = T)
#conferindo se retirando o Al os VIFs estao todos menores que 10
vif.cca(vare_cca)


# 2 -----------------------------------------------------------------------

#Construa um triplot para a cca feita

#  ------------------------------------------------------------------------

plot(vare_cca, scaling = 2, choices = c(1,2), type = "n")
text(vare_cca, display = "cn")
points(vare_cca, pch = 21, col = "red", bg = "orange", cex = 1.5)
text(vare_cca, "species", col = "blue", cex = 1)

# Analise de Redundancia (RDA) -----------------------------------------------------------------------

#Faca uma RDA dos dados varespec (valor de cobertura de especies) e varechem (variaveis ambientais)
data(varespec)
varespec_hel <- decostand(varespec, "hellinger")
data(varechem)

#  ------------------------------------------------------------------------
vare_rda <- rda(varespec_hel ~ N + P + Ca + Mg + S + Al + Fe + Mn + Zn + Mo + Baresoil + Humdepth + pH, data = varechem, scale = T)
summary(vare_rda)

# 1 -----------------------------------------------------------------------

#Usando o Valor de Inflacao da Variancia (VIF), verifique se ha variaveis redundantes. Se houver, retire-as e refaca o modelo

#  ------------------------------------------------------------------------

vif.cca(vare_rda) #S e Al sao redundantes, vamos retirar a com maior valor de VIF, ou seja, o Al
#Refazendo o modelo sem Al
vare_rda <- rda(varespec_hel ~ N + P + Ca + Mg + S + Fe + Mn + Zn + Mo + Baresoil + Humdepth + pH, data = varechem, scale = T)
#conferindo se retirando o Al os VIFs estao todos menores que 10
vif.cca(vare_rda)


# 2 -----------------------------------------------------------------------

#Construa um triplot para a cca feita

#  ------------------------------------------------------------------------

plot(vare_rda, scaling = 2, choices = c(1,2), type = "n")
text(vare_rda, display = "cn")
points(vare_rda, pch = 21, col = "red", bg = "orange", cex = 1.5)
text(vare_rda, "species", col = "blue", cex = 1)




# 3 -----------------------------------------------------------------------


#Usando o teste de Anova, compare os resultados da RDA e da CCA que voce acabou de fazer e indique qual o teste mais adequado
anova(vare_rda, by = "axis", step = 999)
anova(vare_cca, by = "axis", step = 999)
