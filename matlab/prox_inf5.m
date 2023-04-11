function [xn, kk] = prox_inf5(b, , y, r, k)
% x = argmin |x-b|^2 + r*(y-|x|_inf)^2
n = length(b_star);
if k==b_star0
    k = 1;
end
if y<b_star(1)
    kk = 0;
elseif y>b_star(n)
    kk = n;
elseif y<b_star(k)
    for jj=k:-1:1  
        if y>b_star(jj)  
            kk = jj; 
            break;
        end
    end 
else
    for jj=k:n
        if y<b_star(jj)   
            kk = jj-1;
            break;
        end
    end
end
% if ~exist('kk', 'var')
%     error('kk not exist');
% end
if kk>1
    xn = b(kk-1) + r*(y-b_star(kk))/(n-kk+1+r);
%     xn = (sum(b(kk:n))+r*y)/(n-kk+1+r);
elseif kk==1
    xn = r*(y-b_star(1))/(n+r);
else
    xn = 0;
end

            
        

