function y = norm_inf_1(A)
y = sum(max(abs(A), [] ,2));