
#1. Simulação dos Dados
#Criaremos um conjunto de dados simulados para um delineamento em blocos casualizados com 30 híbridos avaliados em 5 réplicas, introduzindo um desbalanceamento removendo 20% das observações de forma aleatória.


# Carregando os pacotes necessários
library(lme4)
library(lmerTest)
library(tidyverse)

set.seed(123)

# Definindo os números de híbridos e réplicas
n_hybrids <- 30
n_reps <- 5

# Criando os fatores para réplicas e híbridos
Replicate <- factor(rep(1:n_reps, each = n_hybrids))
Hybrid <- factor(rep(1:n_hybrids, times = n_reps))

# Montando o data frame inicial
data_sim <- data.frame(Replicate, Hybrid)

# Para simular um experimento com desbalanceamento, removemos 20% das observações
data_sim <- data_sim %>% sample_frac(0.8)

# Parâmetros para a simulação:
# Efeito fixo das réplicas (valores arbitrários)
rep_effect <- c(2, -1, 0, 1, -2)
# Variância do efeito aleatório dos híbridos: sigma_g² = 5
sigma_g <- sqrt(5)
# Efeito aleatório dos híbridos
random_hybrid <- rnorm(n_hybrids, mean = 0, sd = sigma_g)
# Variância residual: sigma_e² = 10
sigma_e <- sqrt(10)

# Simulando a resposta (yd) como soma do efeito da réplica, efeito do híbrido e erro residual
data_sim$response <- rep_effect[as.numeric(data_sim$Replicate)] +
  random_hybrid[as.numeric(data_sim$Hybrid)] +
  rnorm(nrow(data_sim), mean = 0, sd = sigma_e)

# Visualizando as primeiras linhas dos dados simulados
head(data_sim)
#2. Ajuste do Modelo Misto com lme4
#Utilizamos o lmer para ajustar o modelo misto, com efeito fixo das réplicas e efeito aleatório dos híbridos:
  

# Ajustando o modelo misto
model <- lmer(response ~ Replicate + (1 | Hybrid), data = data_sim)

# Sumário do modelo
summary(model)

#No sumário, os efeitos fixos (réplicas) e os componentes de variância dos efeitos aleatórios são apresentados.

#3. Extração dos Componentes de Variância e Cálculo da Herdabilidade
#Extraímos os componentes de variância para o efeito dos híbridos e o erro residual, e em seguida estimamos a herdabilidade em nível individual e para as médias dos híbridos.


# Extraindo os componentes de variância
var_comps <- as.data.frame(VarCorr(model))
var_comps

# Variância do efeito dos híbridos (σ²_g) e erro residual (σ²_e)
sigma2g <- as.numeric(var_comps[var_comps$grp == "Hybrid", "vcov"])
sigma2e <- sigma(model)^2  # sigma(model) retorna o erro residual

# Variância fenotípica individual
sigma2f <- sigma2g + sigma2e

# Herdabilidade em nível individual: H² = σ²_g / σ²_f
H2_ind <- sigma2g / sigma2f; H2_ind

# Calculando o número médio de repetições por híbrido (observando o desbalanceamento)
rep_per_hybrid <- table(data_sim$Hybrid)
mean_reps <- mean(rep_per_hybrid)

# Herdabilidade para as médias dos híbridos: H²_media = σ²_g / (σ²_g + σ²_e / r)
H2_mean <- sigma2g / (sigma2g + sigma2e / mean_reps); H2_mean

# Exibindo os resultados
data.frame(
  Parâmetro = c("Variância dos Híbridos (σ²_g)",
                "Variância Residual (σ²_e)",
                "Variância Fenotípica (σ²_f)",
                "Herdabilidade (individual)",
                "Herdabilidade (médias)"),
  Valor = c(round(sigma2g, 3),
            round(sigma2e, 3),
            round(sigma2f, 3),
            round(H2_ind, 3),
            round(H2_mean, 3))
)

#4. Obtenção dos BLUPs com Erros Padrão Usando merTools
#O pacote merTools permite estimar, por simulação, os erros padrão 
#associados aos BLUPs dos efeitos aleatórios. 
#Assim, podemos avaliar a confiabilidade dos BLUPs.


# Carregando o pacote merTools
library(merTools)

# Extraindo os BLUPs dos híbridos (efeitos aleatórios)
blups <- ranef(model)$Hybrid
blups <- data.frame(Hybrid = rownames(blups), BLUP = blups[,1])
rownames(blups) <- NULL


# Exibindo os 10 híbridos com maiores BLUPs e seus respectivos erros padrão e confiabilidade
sel = blups %>% arrange(desc(BLUP)) %>% head(10)
DS = mean(as.numeric(sel$BLUP)) - mean(as.numeric(blups$BLUP))
DS


GS = DS*H2_mean
GS
