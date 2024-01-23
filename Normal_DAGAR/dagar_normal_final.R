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
  (semilla <- choose_seed(n = 60, p = 0.034, N = 100000))
  
  # semilla: 16213
  
  # fijamos la semilla
  set.seed(semilla)
  
  # obtenemos un gráfico dagar
  grafo <- sample_gnp(n=60, p=0.034, directed = TRUE, loops = FALSE)
  
  # comprobamos si es dag
  is.dag(grafo)
  
  # podemos visualizar el grafo anterior
  plot(grafo, edge.arrow.size=0.01,vertex.size = 0.5, xlab = "Erdos-Renyi Model")
  
  # obtenemos una matriz de adyacencia
  M_adj <- as_adjacency_matrix(grafo) |> as.matrix()
  M <- M_adj + t(M_adj)
  M[M != 0] <- 1
  
  # datos necesarios para la implementacion
  N_edges <- sum(M)/2 ### number of edges
  nei <- neighbors(M) 
  adj.ends <- adj_index(M)$adj_index 
  N_nei <- adj_index(M)$N_nei
  
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
  rho <- 0.9
  
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
  
  filas <- 50
  matriz_beta0  <- matrix(ncol = 6, nrow = filas)
  matriz_beta1  <- matrix(ncol = 6, nrow = filas)
  matriz_beta2  <- matrix(ncol = 6, nrow = filas)
  matriz_sigmae <- matrix(ncol = 6, nrow = filas)
  matriz_sigmau <- matrix(ncol = 6, nrow = filas)
  matriz_rho    <- matrix(ncol = 6, nrow = filas)
  # matriz_rhat   <- matrix(ncol = 6, nrow = 8)
  col_select = c(1, 3, 4, 6, 8, 10)
  
  for(i in 1:50){
    # obtenemos la realizaciones de u
    # simulamos y
    cat("Iteración:",i, "\r")
    y <- rnorm(n = N,
               mean = b0 + b1 * z1 + b2 * z2 + w,
               sd = sigma_e)
    
    # Ajustamos el modelo DAGAR
    dagar_data <- list(N = N, N_edges = N_edges, nei = nei, 
                       adjacency_ends = adj.ends, y = y, N_nei = N_nei, z1 = z1, 
                       z2 = z2, M = M)
    
    fit_dagar <- stan(file = 'Codigos/Dagar/Implementacion/normal_dagar.stan', 
                           data = dagar_data, chains = 3,
                           iter = 14000, thin = 5,
                           pars = c("beta0", "beta1", "beta2", "sigma2_e", 
                                    "sigma2_u", "rho"),
                           cores = 4)
    
    resumen <- summary(fit_dagar)$summary
    
    matriz_beta0[i,]  <- resumen[1, col_select]
    matriz_beta1[i,]  <- resumen[2, col_select]
    matriz_beta2[i,]  <- resumen[3, col_select]
    matriz_sigmae[i,] <- resumen[4, col_select]
    matriz_sigmau[i,] <- resumen[5, col_select]
    matriz_rho[i,]    <- resumen[6, col_select]
  }

a <- rstan::extract(fit_dagar)

pairs(fit_dagar)

matriz_beta0 |> apply(2, mean)
matriz_beta1 |> apply(2, mean)
matriz_beta2 |> apply(2, mean)
matriz_sigmae |> apply(2, mean)
matriz_sigmau |> apply(2, mean)
matriz_rho |> apply(2, mean)

traceplot(fit_dagar)

summary(fit_dagar)

warnings()

resultados_dagar <- list(b0 =  matriz_beta0,
                         b1 =  matriz_beta1,
                         b2 =  matriz_beta2,
                         sigmae =  matriz_sigmae,
                         sigmau =  matriz_sigmau,
                         rho = matriz_rho
                         )

saveRDS(resultados_dagar, "resultados_dagar92.RDS")


resumen