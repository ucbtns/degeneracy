function av = average_quantity(quantity, Nf, t, n, X)

for x = 1:X
   quantity_lower = quantity{x} ;
   av = zeros(Nf,t);
    for f = 1:Nf
        for j = 1:t
            collector =0;
            for i = 1:n
                   collector = collector + quantity_lower{i}(f,j);
            end
            av(f,j) = av(f,j) + collector/n;
        end
    end
    av2{x} =  av;
end

av=0;
for x = 1:X
    av = av + av2{x};  
end

av = av/(X);
av = sum(sum(av));
return 