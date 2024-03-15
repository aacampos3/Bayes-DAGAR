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

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("Codigos/Dagar/Implementacion/funciones_dagar.R")


# Modelo DAGAR


# Matriz mapa -------------------------------------------------------------

vector_unos <- function(vector, n = 15){
  
  # si el vector tal que sea mayor que la cantidad de regiones arroja un "error"
  if(sum(vector > n)){
    return("Hay una número mayor que el largo de los datos")
  }
  
  resultado <- vector(mode = "numeric", length = n)
  
  for(v in vector){
    resultado[v] <- 1
  }
  return(resultado)
}

vector_unos(c(1,4), n =5)

# En primer lugar definiremos la matriz de pesos
M <- rbind(vector_unos(c(5)), # 1
       vector_unos(c(3)), # 2
       vector_unos(c(2, 4, 5)), # 3
       vector_unos(c(3, 5, 15)), # 4
       vector_unos(c(1, 3, 4, 6, 15)), # 5
       vector_unos(c(5, 7, 15)), # 6
       vector_unos(c(6)), # 7
       vector_unos(c(9, 14)), # 8
       vector_unos(c(8, 13, 14)), # 9
       vector_unos(c(11, 13)), # 10
       vector_unos(c(10, 12, 13)), # 11
       vector_unos(c(11)), # 12 a
       vector_unos(c(9, 10, 11)), # 12 b
       vector_unos(c(8, 9)), # 12 c
       vector_unos(c(4, 5, 6))  # 12 d
)



# Modelo DAGAR on grafo de vecindad ---------------------------------------

semilla <- choose_seed(n = 100, p = 0.015)

set.seed(semilla)

grafo <- sample_gnp(n=100, p=0.015, directed = TRUE, loops = FALSE)

plot(grafo, edge.arrow.size=0.01,vertex.size = 0.5, xlab = "Erdos-Renyi Model")

M_adj <- as_adjacency_matrix(grafo) |> as.matrix()

M <- M_adj + t(M_adj)

M[M != 0] <- 1


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

N_edges <- sum(M)/2 ### number of edges
nei <- neighbors(M) 
adj.ends <- adj_index(M)$adj_index 
N_nei <- adj_index(M)$N_nei

# como definir los w?
w <- vector(mode  = "numeric", length = N)
w[1] <- rnorm(1, mean = 0, sd = sigma_u)

vec_var = (1 - rho^2) / (1 + (N_nei - rep(1,N)) * rho^2 )
b = rho / (1 + (N_nei - rep(1,N)) * rho^2 )

for(i in 2:N){
  w[i] <- rnorm(1,
                mean = sum(w[nei[(adj.ends[i-1]+1):adj.ends[i]]]) * b[i],
                sd = sqrt(vec_var[i]) * sigma_u)
}

# datos
# añadir
y <- vector(mode = "numeric", length = N)
for(i in 1:N){
  y[i] <- rnorm(1,
                mean = b0 + b1 * z1[i] + b2 * z2[i] + w[i],
                sd = sigma_e)
}

dagar_data <- list(N = N, N_edges = N_edges, nei = nei, 
                   adjacency_ends = adj.ends, y = y, N_nei = N_nei, z1 = z1, 
                   z2 = z2, M = M)

fit_dagar_grid <- stan(file = 'Codigos/Dagar/Implementacion/normal_dagar.stan', 
                       data = dagar_data, chains = 3,
                       iter = 10, 
                       pars = c("beta0", "beta1", "beta2", "sigma2_e", "sigma2_u",
                                "rho"),
                       cores = 4)

dfdiresumen <- fit_dagar_grid |> summary()


# Grafo dirigido ----------------------------------------------------------


set.seed(1234)
graph <- graph_from_adjacency_matrix(M_adj, mode = "directed", diag = FALSE)
order <- topo_sort(graph)

A_reordenada <- M_adj[order, order]

(M2 <- t(A_reordenada))

N_edges <- sum(M2) ### number of edges
nei <- neighbors(M2) 
adj.ends <- adj_index(M2)$adj_index 
N_nei <- adj_index(M2)$N_nei

# como definir los w?
w <- vector(mode  = "numeric", length = N)
w[1] <- rnorm(1, mean = 0, sd = sigma_u)

vec_var = (1 - rho^2) / (1 + (N_nei - rep(1,N)) * rho^2 )
b = rho / (1 + (N_nei - rep(1,N)) * rho^2 )

for(i in 2:N){
  w[i] <- rnorm(1,
                mean = sum(w[nei[(adj.ends[i-1]+1):adj.ends[i]]]) * b[i],
                sd = sqrt(vec_var[i]) * sigma_u)
}

# datos
# añadir
y <- vector(mode = "numeric", length = N)
for(i in 1:N){
  y[i] <- rnorm(1,
                mean = b0 + b1 * z1[i] + b2 * z2[i] + w[i],
                sd = sigma_e)
}

dagar_data <- list(N = N, N_edges = N_edges, nei = nei, 
                   adjacency_ends = adj.ends, y = y, N_nei = N_nei, z1 = z1, z2 = z2)

fit_dagar_grid <- stan(file = 'Codigos/Dagar/Implementacion/normal_dagar.stan', 
                       data = dagar_data, chains = 3,
                       iter = 10000, 
                       pars = c("beta0", "beta1", "beta2", "sigma2_e", "sigma2_u",
                                "rho"),
                       cores = 4)

resumen <- fit_dagar_grid |> summary()

print_results(resumen)


# Grafo aleatorio ---------------------------------------------------------

semilla <- choose_seed(n = 50, p = 0.02)

set.seed(semilla)

grafo <- sample_gnp(n=50, p=0.02, directed = TRUE, loops = FALSE)

plot(grafo, edge.arrow.size=0.01,vertex.size = 0.5, xlab = "Erdos-Renyi Model")

M_adj <- as_adjacency_matrix(grafo) |> as.matrix()

M <- M_adj + t(M_adj)

M[M != 0] <- 1


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
rho <- 0.5

N_edges <- sum(M)/2 ### number of edges
nei <- neighbors(M) 
adj.ends <- adj_index(M)$adj_index 
N_nei <- adj_index(M)$N_nei

# como definir los w?
w <- vector(mode  = "numeric", length = N)
w[1] <- rnorm(1, mean = 0, sd = sigma_u)

vec_var = (1 - rho^2) / (1 + (N_nei - rep(1,N)) * rho^2 )
b = rho / (1 + (N_nei - rep(1,N)) * rho^2 )

for(i in 2:N){
  w[i] <- rnorm(1,
                mean = sum(w[nei[(adj.ends[i-1]+1):adj.ends[i]]]) * b[i],
                sd = sqrt(vec_var[i]) * sigma_u)
}

# datos
# añadir
y <- vector(mode = "numeric", length = N)
for(i in 1:N){
  y[i] <- rnorm(1,
                mean = b0 + b1 * z1[i] + b2 * z2[i] + w[i],
                sd = sigma_e)
}

dagar_data <- list(N = N, y = y, z1 = z1, z2 = z2, M = M)

fit_dagar_grid <- stan(file = 'Codigos/Dagar/Implementacion/random_normal_dagar.stan', 
                       data = dagar_data, chains = 3,
                       iter = 20000, 
                       pars = c("beta0", "beta1", "beta2", "sigma2_e", "sigma2_u",
                                "rho", "prob"),
                       cores = 4)

(resumen <- fit_dagar_grid |> summary())

traceplot(fit_dagar_grid)