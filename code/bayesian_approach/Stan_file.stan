functions {
  vector predict_eq(matrix B1,       
               vector Ei,
               int N,
               vector y) {

    matrix[N,N] my_B;
    vector[N] predicted;
    real my_pred;

    for(i in 1:N){
      for(j in 1:N){
        my_B[i,j] = B1[i,j];
      }
    }

    for(i in 1:N){
      if(Ei[i]==0){
        for(j in 1:N){
          my_B[i,j]=0;
          my_B[j,i]=0;
        }
        my_B[i,i]=1;
      }
    }
    predicted = inverse(my_B)*y;

    for(i in 1:N){
      if(predicted[i]<=0){
        my_pred = predicted[i];
        predicted[i] = 10^(-20+my_pred);
      }
    }

    return predicted;

  }
}
data {
  int N;                      // number of species
  int M;                      // number of endpoints
  matrix[M,N] E;              // matrix of equilibrium    
  matrix[N,N] B_upper;
  matrix[N,N] B_lower;
  vector[N] y;
  real maxB;                // vector of -1  
}
parameters {
  matrix[N,N] B;    // the regression parameters
  vector<lower=0, upper=1.5>[N] sigmax;    // lognormal error, per species 
}
model {
  vector[N] my_x; // vector of predicted abundances

   for(i in 1:N){
      sigmax[i] ~ normal(0,0.25)T[0,1.5]; //normal(sigmax_prior[i],se_sigmax_prior[i]);
  }

  for(i in 1:N){
    for(j in 1:N){
      if(B_upper[i,j]==0){
        B[i,j] ~ normal(0,maxB)T[,0];
      }
      else if(B_lower[i,j]==0){
        B[i,j] ~ normal(0,maxB)T[0,];
      }        
      else{
        B[i,j] ~ normal(0,maxB);
      }
    }
  }


    for(i in 1:M){
      my_x = predict_eq(B, to_vector(E[i,]), N, y);
      for(j in 1:N){
        if(E[i,j]>0){
          target +=  lognormal_lpdf(E[i,j] | log(my_x[j]), sigmax[j]); 
        }
      }
    }

}

