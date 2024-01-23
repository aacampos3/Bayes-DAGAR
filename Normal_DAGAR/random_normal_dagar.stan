functions{
  // funcion  ue entrega el numero de vecinos;
  int n_neighbors(matrix M) {

    int dim;
    dim = rows(M);
    int resultado;
    resultado = 0;
    
    for(i in 1:dim){
      for(j in 1:i){
        if(M[i,j] == 1){
          resultado = resultado +1;
        }
      }
    }
  return resultado;
  }
  
  // Entrega los vecinos de cada uno una de las regiones entregados en un vector
  // que concatena los vecinos.
  int[] neighbors(matrix M, int dim){
    
    int n_nei;
    //int dim;
    int pos;
  
    n_nei = n_neighbors(M);
    //dim = rows(M);
    
    int vector_nei[n_nei];
    pos = 1;
    
    for(i in 1:dim){
      for(j in 1:i){
        if(M[i,j] == 1){
        vector_nei[pos] = j;
        pos = pos + 1;
        }
      }
    }

    return vector_nei;
  }
  // Entrega dos valores el primero es la cantidad de vecinos por region
  vector N_vecinos(matrix M, int dim){
    //int dim;
    //dim = rows(M);
    
    matrix[dim, dim] M_aux;

    M_aux = M;
    
    for(i in 1:dim){
      for(j in i:dim){
        M_aux[i,j] = 0;
      }
    }
    
    vector[dim] vector_unos;
    vector[dim] total_nei;
    
    vector_unos = rep_vector(1, dim);
    total_nei = M_aux * vector_unos;
    
    return total_nei;
  }
  // Lo mismo que la funcion N_vecinos pero entrega un int array
  int[] Nint_neighbors(matrix M, int dim){
    //int dim;
    //dim = rows(M);

    int total_nei[dim];
    
    for (i in 1:dim) {
        total_nei[i] = 0;
        for (j in 1:i) {
            if (M[i, j] == 1) {
                total_nei[i] = total_nei[i] + 1;
            }
        }
    }
    return total_nei;
  }
  
  // Entrega la suma acumulada de Nint_neighbors
  int[] adjacency_index(int[] total_nei){
    int dim = num_elements(total_nei);
    int resultado[dim];

    for (i in 1:dim) {
      if (i == 1) {
        resultado[i] = total_nei[i];
      } else {
        resultado[i] = resultado[i - 1] + total_nei[i];
      }
    }
    return resultado;
  }
}
 
data {
  int <lower=1> N;
  // vector[N] N_nei; // Número de vecinos dirigidos para cada nodo
  // int <lower=1> N_edges;  // número total de aristas
  // int <lower=1, upper=N-1> nei[N_edges]; // Vector apilando los vecinos dirigidos de cada nodo
  // int<lower=0, upper=N_edges> adjacency_ends[N]; //Donde la adyacencia de cada nodo termina en el vector nei
  // 
  // int <lower=0, upper=1> conexion[N_edges];
  int M[N,N];

  real y[N];
  real z1[N];
  real z2[N];
}

transformed data{
  // la matriz original tiene formato int array de dos dimensiones debido a que 
  // las entradas son binarias pero para el calculo se necesita formato matricial
  matrix[N,N] M_m;
  M_m = to_matrix(M);
}

parameters {
  real <lower=0, upper=1> rho; // Spatial correlation parameter
  real <lower=0> sigma2_u; // Precision of the spatial effects
  real <lower=0> sigma2_e; // Error general
  //real <lower=0> tau; // Precision of the spatial effects
  vector[N] w; // Spatial Random Effect
  real beta0; // intercept
  real beta1;
  real beta2;
  real <lower=0, upper=1> prob;
}

model {
  vector[N] b; //
  vector[N] vec_var; //
  vector[N] t_rowsum; // only the rowsum of t is used
  vector[N] std_dev; // Rescaled std_dev by std_dev_w
  
  // definimos la versomilitud de la matriz de adyacencia, es decir,
  // agregamops incertidumbre a las conexiones vistas
  for(i in 1:N){
    for(j in 1:(i-1)){
      M[i, j] ~ binomial(1, prob);
    }
  }
  
  // funciones para poder implementar y calcular los pesos.
  int N_edges;
  vector[N] N_nei;
  
  N_nei = N_vecinos(M_m, N);
  N_edges = n_neighbors(M_m);
  
  int nei[N_edges];
  int adjacency_ends[N];
  
  nei = neighbors(M_m, N);
  int aux[N];
  aux = Nint_neighbors(M_m, N);

  adjacency_ends = adjacency_index(aux);

  // Construct w
  //vec_var = (1 - rho * rho) ./ (1 + (N_nei - rep_vector(1,N)) * rho * rho);
  //b = rho ./ (1 + (N_nei - rep_vector(1,N)) * rho * rho );
  
  // Construct w de las distribuciones condicionales de W
  vec_var = (1 - rho * rho) ./ (1 + (N_nei - rep_vector(1,N)) * rho * rho);
  b = rho ./ (1 + (N_nei - rep_vector(1,N)) * rho * rho );

  // Linear in number of edges
  t_rowsum[1] = 0;
  for(i in 2:N){
    t_rowsum[i] = sum(w[nei[(adjacency_ends[i-1] + 1):adjacency_ends[i]]]) * b[i];
  }

  std_dev = sqrt(sigma2_u) * sqrt(vec_var);
  
  // definios la verosimilitud de los efectos aleatorios y de las observaciones
  for(i in 1:N){
    w[i] ~ normal(t_rowsum[i], std_dev[i]);
    y[i] ~ normal(beta0 + beta1 * z1[i] + beta2 * z2[i] + w[i], sqrt(sigma2_e));
  }

  // definimos las prioris del modelo
  rho ~ beta(1,1);
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);
  sigma2_u ~ inv_gamma(0.01, 0.01);
  sigma2_e ~ inv_gamma(0.01, 0.01);
  prob ~ beta(1,1);
}

// linea en blanco;
