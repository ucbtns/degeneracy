function [o, av_o,  s_o] = calculate_average(data, X, n, t)

if t == 4
    n_o = data(:, 10:13);
elseif t == 3
    n_o = data(:,8:10);
end

for i = 1:length(n_o)
    o(i, 1) = sum(unique(n_o(i,:)) ~= 0);
    if ( n_o(i,1) == n_o(i,3) ||  n_o(i,1) == n_o(i,2))
        o(i,2) = 1;
    end
end

% Individual Total:
av_o = zeros(X, 1);
%num = 490;
num = 1;

for i = 1:X
    av_o(i,1) =  sum(o(num:(num+(n-1)),2))/(n);
    %av_o(i,1) =  sum(o(num:(num+(9)),2))/(10);
    num = num + n; 
end

% Total:
s_o = mean(av_o);
