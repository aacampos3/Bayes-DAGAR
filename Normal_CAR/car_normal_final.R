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
source("funciones_dagar.R")

# elegimos una semilla tal que se compruebe que la matriz sea DAG
semilla <- choose_seed(n = 60, p = 0.034, N = 1000000)

# fijamos la semilla
set.seed(semilla)

# obtenemos un gráfico dagar
grafo <- sample_gnp(n=60, p=0.034, directed = TRUE, loops = FALSE)

# podemos visualizar el grafo anterior
plot(grafo, edge.arrow.size=0.01,vertex.size = 0.5, xlab = "Erdos-Renyi Model")

# obtenemos una matriz de adyacencia
M_adj <- as_adjacency_matrix(grafo) |> as.matrix()
W <- M_adj + t(M_adj)
W[W != 0] <- 1

# contamos la cantidad de vecinos
apply(W, MARGIN = 2, sum)

# Funcion que transforma la matriz
pesos_matriz <- function(W){
  n <- ncol(W)
  Wt <- matrix(ncol = n, nrow = n)
  for(i in 1:n){
    n_i <- sum(W[i,])
    for(j in 1:n){
      Wt[i,j] <- W[i,j]/n_i
    }
  }
  return(Wt)
} 

Wt <- pesos_matriz(W)

# Definimos rho, sigma_e y sigma_u
rho <- 0.1
sigma_u <- 1
sigma_e <- sqrt(10)

# definimos el tamaño de la matriz
n <- 60

# matriz de varianzas covarianzas de los efectos aleatorios
W1 <- diag(as.numeric(W %*% matrix(rep(1, ncol(W)))), ncol = ncol(W), nrow = ncol(W))
# obtenemos la matriz de precision
S <- (1/sigma_u)^2 * (W1 - rho * W)


# simulamos los valores z1 y z2 de una distribucion Normal(0, 1) que corresponden a los
# "datos" del modelo
z1 <- rnorm(n)
z2 <- rnorm(n)

# definimos los valores de beta_0, beta_1 y beta_2 y t_e
beta_0 <- 4
beta_1 <- 2
beta_2 <- -1

# creamos matrices para guardar los resultados
filas <- 3750
matriz_beta0 <- matrix(ncol = 10, nrow = filas)
matriz_beta1 <- matrix(ncol = 10, nrow = filas)
matriz_beta2 <- matrix(ncol = 10, nrow = filas)
matriz_te <- matrix(ncol = 10, nrow = filas)
matriz_tu <- matrix(ncol = 10, nrow = filas)
matriz_rho <- matrix(ncol = 10, nrow = filas)
matriz_covergencia <- matrix(ncol = 10, nrow = 7)

rownames(matriz_covergencia) <- c("b0", "b1", "b2", "sigma_e", "tu", "rho", 
                                  "sigma2_e")

# Ciclo for ---------------------------------------------------------------

# fijamos la semilla
set.seed(semilla)

# simulamos el efecto aleatorio
u <- mvtnorm::rmvnorm(1, mean = rep(0, n), sigma = solve(S))

# tambien se puede guardar en un archivo RDS
#u <- readRDS(file = "Normal_CAR/efecto_aleatorio.rds")


for(i in 1:10){
  # obtenemos la realizaciones de u
  # simulamos y
  cat("Iteración:",i, "\r")
  
  y1 <- rnorm(n = n,
              mean = beta_0 + beta_1 * z1 + beta_2 * z2 + u,
              sd = sigma_e)
  
  # creamos una lista con los datos necesarios para el modelo en stan
  datos <- list(N = length(y1),
                y = y1,
                z1 = z1,
                z2 = z2,
                W = W,
                dWI = W1)
  
  # ejecutamos el modelo
  mod <- stan(file = 'Normal_CAR/normal_car.stan',
              data = datos,
              chains = 3, # tres cadenas
              iter = 8000, # 8000 iteraciones
              thin = 2, # thining 2
              # parametros a retornar
              pars = c("b0", "b1", "b2", "sigma_e", "sigma_u", "rho"))
  
  # resumen del modelo
  aux <- summary(mod)
  
  # guardamos la convergencia
  (convergencia <- aux$summary[,10])
  
  # extraemos las cadenas del modelo
  cadena <- rstan::extract(mod)
  
  # guardamos los resultados de cada resultado
  matriz_beta0[,i] <- cadena$b0
  matriz_beta1[,i] <- cadena$b1
  matriz_beta2[,i] <- cadena$b2
  matriz_te[,i] <- median(cadena$sigma_e)
  matriz_tu[,i] <- median(cadena$sigma_u)
  matriz_rho[,i] <- cadena$rho
  matriz_covergencia[,i] <- convergencia
}

# veamos los resultados
# muestra las medias a posteriori de las 10 iteraciones
matriz_beta0 |> apply(2, mean)
matriz_beta1 |> apply(2, mean)
matriz_beta2 |> apply(2, mean)
matriz_sigmae |> apply(2, mean)
matriz_sigmau |> apply(2, mean)
matriz_rho |> apply(2, mean)
matriz_covergencia # muestra la convergencia por iteracion y parametro
