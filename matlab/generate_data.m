
function [opt] = generate_data(m, n, alph)
B = rand(m,n)-0.5;
% uniform distribution [-0.5,0.5]
%% the max of the each row of the matrix B.  norm 1 and infinite .

tau = alph*sum(max(abs(B),[], 2));
%% To ensure the tau <=1;

opt.B =B;
opt.tau = tau;