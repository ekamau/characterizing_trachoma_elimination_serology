//
data {
  int N; // length of dataframe = nrows
  int z[N]; // vector - seropositive
  int n[N]; // vector - tested
  int age[N];
  int m; // no. of clusters
  int<lower = 0, upper = m> clusid[N]; // cluster id
  
}

parameters {
  vector<lower=0, upper = 1>[m] clus_foi; // cluster-level lambda
  real<lower=0, upper = 1> lambda; // eu-level lambda
  
}

transformed parameters {
  vector[N] prob_seropos;
  
  for(i in 1:N){
    prob_seropos[i] = 1 - exp(-(clus_foi[clusid[i]] * age[i]));
  }
  
}

model {
  // priors
  lambda ~ exponential(1); // hyperprior
  clus_foi ~ exponential(1/lambda);
  
  // likelihood
  z ~ binomial(n, prob_seropos);
  
}

generated quantities {
  vector[N] z_sim;
  
  for (i in 1:N) 
    z_sim[i] = binomial_rng(n[i], prob_seropos[i]); 
    
}
