function [yM, yCI95, v] = coni(quantity, n)

yM = mean(quantity);                                
ySEM = std(quantity)/sqrt(n);                             
CI95 = tinv([0.025 0.975], n-1);                    
yCI95 = bsxfun(@times, ySEM, CI95(:));       

v = 0;
%[p, S] = polyfit(1:n, yM, 10);
%v = polyval(p, 1:n);
return