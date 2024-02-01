data {
  int <lower=1> N;
  vector[N] N_nei; // Número de vecinos dirigidos para cada nodo
  int <lower=1> N_edges;  // número total de aristas
  int <lower=1, upper=N-1> nei[N_edges]; // Vector apilando los vecinos dirigidos de cada nodo
  int<lower=0, upper=N_edges> adjacency_ends[N]; //Donde la adyacencia de cada nodo termina en el vector nei

  real y[N]; # datos 
  real z1[N]; # covariable
  real z2[N]; # covariable
}

parameters {
  real <lower=0, upper=1> rho; // Parametro de correlacion espacial
  real <lower=0> sigma2_u; //Varianza del efecto espacial
  real <lower=0> sigma2_e; // Varianza del modelo (marginal)
  //real <lower=0> tau;
  vector[N] w; // Efecto aleatorio espacial
  real beta0; // Intercepto
  real beta1; // Covariable 1
  real beta2; // Covariable 2
}

model {
  vector[N] b; // auxiliar (sirve para realizar calculos)
  vector[N] vec_var; // vector de varianzas del efecto aleatorio
  vector[N] t_rowsum; // suma de filas
  vector[N] std_dev; // varianza del vectoir aleatorio escalada


  // construccion de w
  vec_var = (1 - rho * rho) ./ (1 + (N_nei - rep_vector(1,N)) * rho * rho);
  b = rho ./ (1 + (N_nei - rep_vector(1,N)) * rho * rho );
  // suma de filas
  t_rowsum[1] = 0;
  for(i in 2:N){
    t_rowsum[i] = sum(w[nei[(adjacency_ends[i-1]+1):adjacency_ends[i]]]) * b[i];
  }

  std_dev = sqrt(sigma2_u) * sqrt(vec_var);
  
  # verosimilitud
  for(i in 1:N){
    w[i] ~ normal(t_rowsum[i], std_dev[i]);
    y[i] ~ normal(beta0 + beta1 * z1[i] + beta2 * z2[i] + w[i], sqrt(sigma2_e));
  }

  // distribuciones a priori
  rho ~ beta(1,1);
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);
  sigma2_u ~ inv_gamma(0.01, 0.01);
  sigma2_e ~ inv_gamma(0.01, 0.01);
}
