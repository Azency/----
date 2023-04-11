clc;clear;

m = 20;
n = 10;
alph = 1e-3;

 opt = generate_data(m, n, alph);
 B = opt.B;
 tau = opt.tau;
    
   
  max_iter = 100;
  error = 1e-14;

  [X, FX] = proj_inf1ball6(B, tau, max_iter, error);
