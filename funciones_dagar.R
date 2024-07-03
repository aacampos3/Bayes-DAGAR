
# No todos lo grafos son DAG
# Esta función prueba valores de semilla hasta comprobar que en algún caso sea dag
choose_seed <- function(n, p, N){
  i <- 1
  interruptor <- TRUE
  while(i < N & interruptor){
    set.seed(i)
    cat("Iteracion ", i, "\r")
    
    grafo <- sample_gnp(n = n, p = p, directed = TRUE, loops = FALSE) #n=nodes, p=probability of an edge between two nodes
    
    if(is_dag(grafo)){
      interruptor <- FALSE
    }
    i <- i +1
    
  }

  return(i - 1)
}

# funciones que recibe un summary de stan y entrega los resultados en formato
# pipe
table_fit <- function(fit_modelo, rho, prob = FALSE){
  resultado <- resumen$summary |> as.data.frame()
  
  resultado <- resultado |> select(-c(se_mean, n_eff)) |> 
    slice(1:(n()-1)) |> 
    round(3)
  
  if(prob){
    resultado$valor_real <- c(4, 2, -1, 10, 1, rho, 0.02)
    parametro <- c("\\beta_0", "\\beta_1", "\\beta_2", 
                   "\\sigma2_e", "\\sigma2_u", "\\rho", "p")
  }
  
  else{
    resultado$valor_real <- c(4, 2, -1, 10, 1, rho)
    parametro <- c("\\beta_0", "\\beta_1", "\\beta_2", 
                   "\\sigma2_e", "\\sigma2_u", "\\rho")
  }
  
  rownames(resultado) <- NULL
  
  colnames(resultado) <- NULL
  
  resultado <- cbind(parametro, resultado)
  
  resultado |> knitr::kable(format = "pipe", #escape = FALSE, 
                    col.names = c("Parámetro", "Media", "Desv. est.", "2.5\\%", "25 \\%",
                                  "50 \\%", "75 \\%", "97.5 \\%", "Rhat", "Valor real"))
}

# funcion para obtener los vecinos
neighbors <- function(Minc){
  n <- nrow(Minc)
  unlist(
    lapply(2:n, function(i) which(Minc[i, (1:i)] == 1))
  )
}

# funcion para obtener la cantidad de vecinos por region y la suma acumulada de
# lo anterior
adj_index <- function(Minc){
  Minc.low = Minc
  Minc.low[upper.tri(Minc)] = 0
  list(N_nei = rowSums(Minc.low),
       adj_index = cumsum(rowSums(Minc.low)))
}

# 

resumen_resultados <- function(resultado){
  resumen_iteracion <- list(media = apply(resultado, 2, mean),
                            mediana = apply(resultado, 2, median),
                            sd = apply(resultado, 2, sd),
                            li = apply(resultado, 2, quantile, 0.025),
                            ls = apply(resultado, 2, quantile, 0.975))
  
  resumen_final <- map_vec(resumen_iteracion, mean)
  
  resumen_final <- round(resumen_final, 3)
  
  return(resumen_final)
}

