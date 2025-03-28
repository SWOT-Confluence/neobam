functions {
  // Conversion from array to vector takes place by row.
  // Nested elements are appended in sequence.

  // Convert an array to a vector based on a binary matrix
  // indicating non-missing data
  vector ragged_vec(vector[] x, int[,] bin) {
    vector[num_elements(x)] out;
    int ind;

    ind = 1;
    for (i in 1:size(x)) {
      for (t in 1:num_elements(x[1])) {
        if (bin[i, t] == 1) {
          out[ind] = x[i, t];
          ind += 1;
        }
      }
    }

    // print(out);
    return(out[1:(ind - 1)]);
  }

  // Repeat elements of a "row" vector to match with 2-D array vectorization
  vector ragged_row(vector x, int[,] bin) {
    vector[num_elements(bin)] out;
    int ind;

    ind = 0;
    for (i in 1:size(bin)) {
      for (t in 1:num_elements(bin[1])) {
        if (bin[i, t] == 1) {
          ind += 1;
          out[ind] = x[t];
        }
      }
    }
    return(out[1:ind]);
  }

  // Repeat elements of a "column" vector to match with 2-D array vectorization
  vector ragged_col(vector x, int[,] bin) {
    vector[num_elements(bin)] out;
    int ind;

    ind = 0;
    for (i in 1:size(bin)) {
      for (t in 1:num_elements(bin[1])) {
        if (bin[i, t] == 1) {
          ind += 1;
          out[ind] = x[i];
        }
      }
    }
    return(out[1:ind]);
  }

  // indices of vectorized bin1 that are also in vectorized bin2
  int[] commoninds(int[,] bin1, int[,] bin2) {
    int out[num_elements(bin1)];
    int vecinds[size(bin1), num_elements(bin1[1])];
    int ctr;
    int ind;

    ctr = 0;
    for (i in 1:size(bin1)) {
      for (t in 1:num_elements(bin1[1])) {
        if (bin1[i, t] == 1) {
          ctr += 1;
          vecinds[i, t] = ctr;
        }
      }
    }

    ind = 0;
    for (i in 1:size(vecinds)) {
      for (t in 1:size(vecinds[1])) {
        if (bin2[i, t] == 1) {
          ind += 1;
          out[ind] = vecinds[i, t];
        }
      }
    }
    return(out[1:ind]);
  }
  
  vector vector_pow(vector vector1, vector vector2){
      vector[num_elements(vector1)] out;
     int N;
     N=num_elements(vector1);
     
     for (i in 1:N) {
        out[i] = pow(vector1[i],vector2[i]);
      }
      
      return(out);
  
  }
}

data {
  // Dimensions
  int<lower=0> nx; // number of reaches
  int<lower=0> nt; // number of times
  int<lower=0> ntot; // total number of non-missing widths

  // Missing data
  int<lower=0,upper=1> hasdat[nx, nt]; // matrix of 0 (missing), 1 (not missing)

  // *Actual* data
  vector[nt] Sobs[nx];
  //vector[nt] Wobs[nx];
  vector[nt] Hobs[nx];
  
  real<lower=0> Serr_sd;
  //real<lower=0> Werr_sd;
  real<lower=0> Herr_sd;
  
  real<lower=0> logSerr_sd;
  //real<lower=0> logWerr_sd;
  real<lower=0> logHerr_sd;
  
  //geophysical
  real stationvec[nx];

  // Hard bounds on parameters
  real lowerbound_logQ;
  real upperbound_logQ;

  real lowerbound_r;
  real upperbound_r;

  real lowerbound_logn;
  real upperbound_logn;
  
  real lowerbound_logDb;
  real upperbound_logDb;
  
  real lowerbound_logWb;
  real upperbound_logWb;
  
  real lowerbound_Zo1;
  real upperbound_Zo1;
  
  //real lowerbound_Zo;
  //real upperbound_Zo;


  // *Known* likelihood parameters
  vector<lower=0>[nt] sigma_man[nx];

  // Hyperparameters
  vector[nt] logQ_hat; // prior mean on logQ

  real logWb_hat[nx];
  real r_hat[nx];
  real logn_hat[nx];
  real logDb_hat[nx];
  real Zo1_hat;
  //real Zo_hat[nx];

  vector<lower=0>[nt] logQ_sd;
  real<lower=0> logWb_sd[nx];
  real<lower=0> r_sd[nx];
  real<lower=0> logn_sd[nx];
  real <lower=0> logDb_sd[nx];
  real <lower=0> Zo1_sd;
   //real<lower=0> Zo_sd[nx];
}

transformed data {
  // Transformed data are *vectors*, not arrays. This to allow ragged structure
  // there are three flavors of input data : height, width, and slope
  // since everything is precompiled, this means we need 6 variables as the
  // size of the variable is different depending on whether or not we have widths only
  // or height, width, and slope.
  // further, we need to declare the log transforms of these before taking the logs.
  // we could, I suppose, accept the log transform as an input directly
  vector[ntot] Sobsvec;
  //vector[ntot] Wobsvec;
  vector[ntot] Hobsvec;
  vector[ntot] logSobsvec;
  //vector[ntot] logWobsvec;
  vector[ntot] logHobsvec;
  vector[ntot] sigma_vec_man;

  // convert pseudo-ragged arrays to vectors
  Sobsvec = ragged_vec(Sobs, hasdat);
  //Wobsvec= ragged_vec(Wobs,hasdat);
  Hobsvec= ragged_vec(Hobs,hasdat);
  
   // convert pseudo-ragged arrays to vectors
  //logWobsvec = log(Wobsvec);
  logSobsvec = log(Sobsvec);
  logHobsvec = log(Hobsvec);


  
  sigma_vec_man = ragged_vec(sigma_man, hasdat);
  

}

parameters {
  // what is passed to the model function to do the sampling.
  // in essence, what are the important terms we need to run the MCMC?
  // DOES NOT INCLUDE OBSERVED DATA
  // pure declaration here, no manipulation
  
  //db, wb, r,n, Q, Zo are pure parameters

  vector<lower=lowerbound_logQ,upper=upperbound_logQ>[nt] logQ;
  vector<lower=lowerbound_r, upper=upperbound_r>[nx] r[1];
  vector<lower=lowerbound_logWb, upper=upperbound_logWb>[nx] logWb[1];
  vector<lower=lowerbound_logn, upper=upperbound_logn>[nx] logn[1];
  vector<lower=lowerbound_logDb, upper=upperbound_logDb>[nx] logDb[1];
  //vector<lower=lowerbound_Zo, upper=upperbound_Zo>[nx] Zo[1];
 
  
  //real<lower=lowerbound_r, upper=upperbound_r> r;
  //real<lower=lowerbound_logWb, upper=upperbound_logWb> logWb;
  //real<lower=lowerbound_logn, upper=upperbound_logn> logn;
  //real<lower=lowerbound_logDb, upper=upperbound_logDb> logDb;
  real<lower=lowerbound_Zo1, upper=upperbound_Zo1> Zo1;
 
  // we'll turn on measurement error later, so these need to be declared as theoretical sampling
  // targets
  
  vector<lower=0>[ntot] Sact[1];
  //vector<lower=0>[ntot] Wact[1];
  vector<lower=0>[ntot] Hact[1];
  
}

transformed parameters {
  // in this block we write the algebraic expressions that will form
  // the basis of the sampling later in the 'model' block.
  // we can pull in anything from data, transformed data, or parameters

  //low level language- declare variables
  //in essence, sample Q then repeat it across all nx. Q is an nt variable, so it is sampled at each nt
  
      vector[ntot] logQ_rep[1]; // location-repeated logQ
      logQ_rep[1]=ragged_row(logQ, hasdat);
  
  //make the formula easier to type by putting the LHS together here
  
      vector[nx] LHS[1]; 
      LHS[1] = log(Hobsvec -  Zo_transform[1]).*(  rep_vector(1.66,ntot) + (  rep_vector(1,ntot)./ragged_col(r[1],hasdat) )  );
  
  //make a vector that goes downhill
  //Zo1 is a constant, so we don't need to transform it
  //zo_transform is ntot, as are Sobsvec and stationvec2
  
       vector[ntot] Zo_transform[1];  
       Zo_transform[1] = Zo1 - Sobsvec.*stationvec2;
                
  //now make a RHS given the samples we have above
  
      vector[ntot] RHS[1];
      RHS[1] = (  logQ_rep[1]  ) +
                  (( rep_vector(1,ntot)./ragged_col(r[1],hasdat) ).*ragged_col(logDb[1],hasdat)  ) +
                  (  ragged_col(logn[1],hasdat)  ) -
                  (  ragged_col(logWb[1],hasdat)  ) -
                  (  rep_vector(1.66,ntot).*log(ragged_col(r[1],hasdat) ./ (ragged_col(r[1],hasdat) + 1) ))  ;
  
 
}
model {
  // Priors
  // these are the 'pure priors'. they are only informed indirectly through the sampling of the 
  // flow law later on. Constrained by input operations
  
  //this one then gets location repeated by the variable transform to logQ_rep
  //sample a vector of Q that doesn't follow the flow law, it is just a Q sample
  logQ ~ normal(logQ_hat, logQ_sd);
  r[1] ~ normal(r_hat, r_sd);
  logn[1] ~ normal(logn_hat, logn_sd);
  logWb[1] ~normal(logWb_hat,logWb_sd);
  logDb[1] ~normal(logDb_hat,logDb_sd);
  Zo1 ~normal(Zo1_hat,Zo1_sd);
  
  //the flow laws
  //observations, or combinations of observations, on the LHS
  // read this is as, what is the likelihood of the LHS given the value of the thing on the RHS
  // the RHS is calculated from all of the samples above
  // given whatever is on the LHS, ideally observations

 LHS[1] ~ normal(RHS[1],0.26)
  
    // Latent vars for measurement error
  Sact[1] ~ normal(Sobsvec[1],Serr_sd); // S meas err
  Hact[1] ~ normal(Hobsvec[1],Herr_sd);
  
  

}





