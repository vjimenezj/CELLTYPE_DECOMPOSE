data {
  int<lower=1> N;	// total number of observations
  vector[N] Y; 		// response variable (column vector with all measurements performed in the experiment)
  int<lower=1> G;  	// number of genes
  matrix[N, G] X1;  	// indicator matrix for genes
  int<lower=1> GC; 	// number of genes * (Conditions - 1)
  matrix[N, GC] X2;	// indicator matrix for condition-specific gene effects (differential expression associated to each gene)
  int<lower=1> R;  	// number of cell types
  matrix[N, R] FT;  	// regulatory matrix (association of cell-types with their specific gene markers)
}

transformed data {
}

parameters {
  vector[G] gene;  	// gene mean expression in the reference condition (WT)
  vector[GC] b_gene;  	// gene-level condition specific effects
  vector[R] b_reg;  	// cell-type-level condition specific effects
  real<lower=0> sigma;  // common dispersion parameter for all genes
}

transformed parameters {
  real lprior = 0;  	// prior contributions to the log posterior
  lprior += normal_lpdf(sigma | 0, 0.5)
  - 1 * normal_lpdf(0 | 0, 0.5);
  lprior += normal_lpdf(gene | 10, 3);
  lprior += normal_lpdf(b_gene | 0, 1);
  lprior += normal_lpdf(b_reg | 0, 1);
}

model {
  vector[N] expr_estimation;
  expr_estimation = X1 * gene + X2 * b_gene + FT * b_reg;
  target += normal_lpdf(Y | expr_estimation, sigma);
  target += lprior;
}
