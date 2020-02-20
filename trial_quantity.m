function av = trial_quantity(quantity, Nf, t, n, X)

av = zeros(X,n);
for x = 1:X
            for i = 1:n
                 av(x,i) = sum(sum(quantity{x}{i}));
            end
end
return 