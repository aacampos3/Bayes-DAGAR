data {
  int <lower=1> N;
  vector[N] N_nei; // Número de vecinos dirigidos para cada nodo
  int <lower=1> N_edges;  // número total de aristas
  int <lower=1, upper=N-1> nei[N_edges]; // Vector apilando los vecinos dirigidos de cada nodo
  int<lower=0, upper=N_edges> adjacency_ends[N]; //Donde la adyacencia de cada nodo termina en el vector nei

  real y[N];
  real z1[N];
  real z2[N];
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
}

model {
  vector[N] b; //
  vector[N] vec_var; //
  vector[N] t_rowsum; // only the rowsum of t is used
  vector[N] std_dev; // Rescaled std_dev by std_dev_w


  // Construct w
  vec_var = (1 - rho * rho) ./ (1 + (N_nei - rep_vector(1,N)) * rho * rho);
  b = rho ./ (1 + (N_nei - rep_vector(1,N)) * rho * rho );
  // Linear in number of edges
  t_rowsum[1] = 0;
  for(i in 2:N){
    t_rowsum[i] = sum(w[nei[(adjacency_ends[i-1]+1):adjacency_ends[i]]]) * b[i];
  }

  std_dev = sqrt(sigma2_u) * sqrt(vec_var);

  for(i in 1:N){
    w[i] ~ normal(t_rowsum[i], std_dev[i]);
    y[i] ~ normal(beta0 + beta1 * z1[i] + beta2 * z2[i] + w[i], sqrt(sigma2_e));
  }

  //beta0 ~ normal(0, 100);
  rho ~ beta(1,1);
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);
  sigma2_u ~ inv_gamma(0.01, 0.01);
  sigma2_e ~ inv_gamma(0.01, 0.01);
}