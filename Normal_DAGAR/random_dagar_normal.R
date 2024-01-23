# librerias
library(rstan)
library(mvtnorm)
library(tidyverse)
library(loo)
library(sf)
library(reshape2)
library(ggrepel)
library(spdep)
library(igraph)

# configuraciones de stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# importamos funciones utiles
source("Codigos/Dagar/Implementacion/funciones_dagar.R")

# elegimos una semilla tal que se compruebe que la matriz sea DAG
semilla <- choose_seed(n = 50, p = 0.02)

# fijamos la semilla
set.seed(semilla)

# obtenemos un gráfico dagar
grafo <- sample_gnp(n=50, p=0.02, directed = TRUE, loops = FALSE)

# podemos visualizar el grafo anterior
plot(grafo, edge.arrow.size=0.01,vertex.size = 0.5, xlab = "Erdos-Renyi Model")

# obtenemos una matriz de adyacencia
M_adj <- as_adjacency_matrix(grafo) |> as.matrix()
M <- M_adj + t(M_adj)
M[M != 0] <- 1

# simulamos los datos
set.seed(1234)
# cantidad de regiones
N <- nrow(M)

z1 <- rnorm(N, 0, 1)
z2 <- rnorm(N, 0, 1)

b0 <- 4
b1 <- 2
b2 <- -1
sigma_u <- sqrt(1)
sigma_e <- sqrt(10)
rho <- 0.1

# como definir los w
w <- vector(mode  = "numeric", length = N)
w[1] <- rnorm(1, mean = 0, sd = sigma_u)

vec_var = (1 - rho^2) / (1 + (N_nei - rep(1,N)) * rho^2 )
b = rho / (1 + (N_nei - rep(1,N)) * rho^2 )

for(i in 2:N){
  w[i] <- rnorm(1,
                mean = sum(w[nei[(adj.ends[i-1]+1):adj.ends[i]]]) * b[i],
                sd = sqrt(vec_var[i]) * sigma_u)
}

y <- vector(mode = "numeric", length = N)
for(i in 1:N){
  y[i] <- rnorm(1,
                mean = b0 + b1 * z1[i] + b2 * z2[i] + w[i],
                sd = sigma_e)
}

# Ajustamos el modelo DAGAR
dagar_data <- list(N = N, y = y, z1 = z1, z2 = z2, M = M)

fit_dagar_grid <- stan(file = 'Codigos/Dagar/Implementacion/random_normal_dagar.stan', 
                       data = dagar_data, chains = 3,
                       iter = 20000, 
                       pars = c("beta0", "beta1", "beta2", "sigma2_e", "sigma2_u",
                                "rho", "prob"),
                       cores = 4)

hola <- rstan::extract(fit_dagar_grid,  permuted = FALSE)

hola[,1,1]

# guardamos los resultados del summary
(resumen <- fit_dagar_grid |> summary())

# tabla resumen del ajuste del modelo
table_fit(resumen, rho = ,prob = TRUE)

# realizamos un traceplor de las cadenas
traceplot(fit_dagar_grid)
