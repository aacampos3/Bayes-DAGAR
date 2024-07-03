data {
  int N;
  vector[N] y;
  vector[N] z1;
  vector[N] z2;
  matrix<lower = 0, upper = 1>[N,N] W; //matriz de adyacencia
  matrix[N,N] dWI; //matriz diagnomal 
}

transformed data {
  vector[N] ceros;
  ceros = rep_vector(0, N);
}

parameters {
  real b0;
  real b1;
  real b2;
  real<lower = 0> sigma2_e;
  real<lower = 0> sigma2_u;
  real<lower=-1, upper = 1> rho;
  vector[N] u;
}

transformed parameters {
  real<lower = 0> tu;
  tu = (1/sigma2_u);
  real<lower = 0> sigma_e;
  sigma_e = sqrt(sigma2_e);
}

model {
  u ~ multi_normal_prec(ceros, tu * (dWI - rho * W));
  vector[N] mu;
  // se define la verosimilitud
  for(i in 1:N){
    mu[i] = b0 + b1 * z1[i] + b2 * z2[i] + u[i];
    
    // se define el modelo
    y[i] ~ normal(mu[i], sigma_e);
  }
  
  b0 ~ normal(0, 1000); // priori para b0
  b1 ~ normal(0, 1000); // priori para b1
  b2 ~ normal(0, 1000); // priori para b2
  sigma2_u ~ inv_gamma(0.01, 0.01);
  sigma2_e ~ inv_gamma(0.01, 0.01);
  rho ~ beta(1,1);
}
