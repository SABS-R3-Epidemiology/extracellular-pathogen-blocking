functions {
  array[] real CompoundModel(real t,
                             array[] real y,
                             array[] real theta,
                             array[] real x_r,
                             array[] int x_i) {
    array[5] real dydt;
  
    real beta = theta[1];
    real alpha = theta[2];
    real p = theta[3];
    real gamma = theta[4];
    real r = theta[5];
    real delta = theta[6];
    real phi = theta[7];
    real k = theta[8];
    
    real T = y[1];  // Susceptible cells
    real I = y[2];  // Cells infected with virus
    real V = y[3];  // Free virus
    real P = y[4];  // Cells which have acquired protection by the compound
    real C = y[5];  // Compound level
  
    // T
    dydt[1] = -beta * T * V - alpha * C * T;
    
    // I
    dydt[2] = beta * T * V - alpha * C * I;
  
    // V
    dydt[3] = p * I - gamma * C * V;
  
    // P
    dydt[4] = alpha * C * T + alpha * C * I;
    
    // C
    dydt[5] = r - delta * C * T - delta * C * I - k * C - phi * C * V;
    
    return dydt;
  }

  array[] real CompareModelBasic(real t,
                                 array[] real y,
                                 array[] real theta,
                                 array[] real x_r,
                                 array[] int x_i) {
    array[3] real dydt;
  
    real beta = theta[1];
    real p = theta[2];
    
    real T = y[1];  // Susceptible cells
    real I = y[2];  // Cells infected with virus
    real V = y[3];  // Free virus
  
    // T
    dydt[1] = -beta * T * V;
    
    // I
    dydt[2] = beta * T * V;
  
    // V
    dydt[3] = p * I;
  
    return dydt;
  }

  array[] real CompareModelExtended(real t,
                                    array[] real y,
                                    array[] real theta,
                                    array[] real x_r,
                                    array[] int x_i) {
    array[4] real dydt;
  
    real beta = theta[1];
    real p = theta[2];
    real L = theta[3];
    real e = theta[4];
    real s = theta[5];
    real A1 = theta[6];
    real A2 = theta[7];
    real l1 = theta[8];
    real l2 = theta[9];
    real tstar = theta[10];
    
    real T = y[1];  // Susceptible cells
    real I = y[2];  // Cells infected with virus
    real V = y[3];  // Free virus
    real E = y[4];  // Cells exposed to virus but not yet producing virus
  
    // T
    dydt[1] = -beta * T * V / (L + V);
  
    // I
    dydt[2] = e * E;
  
    // V
    real F;
    if (t <= tstar) {
      F = A1 * exp(l1 * t);
    } else {
      F = A2 * exp(-l2 * t);
    }
    dydt[3] = p * I / (1.0 + s * F);
    
    // E
    dydt[4] = beta * T * V / (L + V) - e * E;
    
    return dydt;
  }
}

    
data {
  int T;  // Number of time points in each time series
  
  array[T] real ts;  // Time points
  
  int N;  // Number of time series (i.e., 6 for 0%, 20%, ... , 100% Wolbachia)
  
  array[T, N] real v;  // Virus over cell area data (corresponding to I compartment)
  
  int model_type; // model_type = 1 : Protective compound
                  // model_type = 2 : Basic comparator model
                  // model_type = 3 : Extended comparator model
                  
  int use_priors; // use_priors = 0 : No priors (maximum likelihood)
                  // use_priors = 1 : With prior (MAP estimate)
}


transformed data {
  array[0] real x_r;
  array[0] int x_i;
  
  int num_theta_params;  // Number of unknown parameters in the theta parameter vector
  int num_updated_theta_params;  // Number of parameters in the vector which will be sent to the ODE RHS in
                                 // the compound model. The other models send theta to the
                                 // RHS directly, and have num_updated_theta_params equal to zero. See 
                                 // comment in transformed parameters.
  if (model_type == 1) {
    num_theta_params = 8;
    num_updated_theta_params = 8;
  } else if (model_type == 2) {
    num_theta_params = 2;
    num_updated_theta_params = 0;
  } else if (model_type == 3) {
    num_theta_params = 10;
    num_updated_theta_params = 0;
  }
  
  int num_y_dims;  // Number of dimensions of the state of the ODE
  if (model_type == 1) {
    num_y_dims = 5;
  } else if (model_type == 2) {
    num_y_dims = 3;
  } else if (model_type == 3) {
    num_y_dims = 4;
  }
}


parameters {
  real<lower=0.8, upper=1.0> area_ratio;
  real<lower=0.001, upper=1000> V_0;
  array[N] real<lower=0, upper=1000> sigma;
  array[num_theta_params] real<lower=0, upper=1000> theta;
  real<lower=0.1, upper=1> colonized_prop;
}


transformed parameters {
  array[T, num_y_dims, N] real y;
  array[num_y_dims, N] real<lower=0> y0_full;  // Initial condition
  
  for (i in 1:N){
    y0_full[1, i] = (1.0 - colonized_prop * 0.2 * (i - 1)) * area_ratio;
    y0_full[2, i] = 0.0;
    y0_full[3, i] = V_0;
    if (model_type == 1) {
      y0_full[4, i] = 0.0;
      y0_full[5, i] = 0.0;
    } else if (model_type == 3) {
      y0_full[4, i] = 0.0;
    }
  }
  
  real t0=-0.0001;
  
  // For compound model, make a theta_complete vector whose fifth entry will
  // be the compound production rate (for the given Wolbachia percentage, 
  // indexed here by i) computed as a function of the parameter theta[5].
  // For the other models, this step is unneeded and theta parameter vector
  // itself will be sent to the ODE.
  array[num_updated_theta_params] real<lower=0> theta_complete;
  
  for (i in 1:N){
    if (model_type == 1) {
      theta_complete[1] = theta[1];
      theta_complete[2] = theta[2];
      theta_complete[3] = theta[3];
      theta_complete[4] = theta[4];
      theta_complete[5] = theta[5] * 0.2 * (i - 1);
      theta_complete[6] = theta[6];
      theta_complete[7] = theta[7];
      theta_complete[8] = theta[8];
    }

    if (model_type == 1) {
      y[,,i] = integrate_ode_rk45(CompoundModel, y0_full[,i], t0, ts, theta_complete, x_r, x_i);
    } else if (model_type == 2) {
      y[,,i] = integrate_ode_rk45(CompareModelBasic, y0_full[,i], t0, ts, theta, x_r, x_i);
    } else if (model_type == 3) {
      y[,,i] = integrate_ode_rk45(CompareModelExtended, y0_full[,i], t0, ts, theta, x_r, x_i);
    }
  }
}


model {
  if (use_priors == 1) {
    theta ~ normal(0, 1) T[0,];
    area_ratio ~ uniform(0.8, 1.0);
    V_0 ~ normal(0, 0.1) T[0.001, 1000];
    sigma ~ normal(0, 0.01) T[0,];
    colonized_prop ~ uniform(0.1, 1.0);
  }

  for (i in 1:N){
    for (j in 1:T) {
      v[j, i] ~ normal(y[j, 2, i], sigma[i]);
    }
  }
}


generated quantities {
  array[T, num_y_dims, N] real simulated_y;
  simulated_y = y;
  real log_likelihood = 0;
  for (i in 1:N){
    for (j in 1:T) {
      log_likelihood += normal_lpdf(v[j, i] | y[j, 2, i], sigma[i]);
    }
  }
}