data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of samples (nrow)
    int<lower=1> Species[S];    // index of species, each of which will have a different value for p11 and p10
    int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
    int<lower=1> L[S];   // index of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> K[S];   // number of replicates per site (ncol)
    int<lower=0> N[S]; // number of detections among these replicates
    int z[S];   // integer flag to help estimate psi parameter
  }
  parameters{/////////////////////////////////////////////////////////////////////
  real<lower=0,upper=1> psi[Nloc];  //commonness parameter
  real<lower=0,upper=1> p11[Nspecies]; //true positive detection rate
  real<lower=0,upper=1> p10[Nspecies]; //false positive detection rate
  }
  transformed parameters{/////////////////////////////////////////////////////////////////////
  }
  model{/////////////////////////////////////////////////////////////////////
  real p[S];

  for (i in 1:S){
  z[i] ~ bernoulli(psi[L[i]]);
  p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
  N[i] ~ binomial(K[i], p[i]);
  }; 

  //priors
  psi ~ beta(2,2); 
  p11 ~ beta(2,2); 
  p10 ~ beta(1,10);
  }
  generated quantities{
  real<lower=0,upper=1> Occupancy_prob[S];    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015

  for (i in 1:S){
  Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  / ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
  + (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
  );
  }
  }

  
