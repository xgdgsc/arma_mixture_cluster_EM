function bic = BIC(posterior, prior, v)
%FUNCTION BIC
% Decide the number of clusters when the # of clusters has not been specified by the user.
% The method uses Baysican information criterion(BIC) to decide which
% cluster_num is best.
%
%INPUT:
% posterior: the posterior probablities generated by EM algorithm.
% prior: the prior probablities.
% v: the number of parameters in each component model, in our AR model, 
%    v = p + q + 2;
%OUTPUT:
% bic: the bic value.

% M: cluster#, N: time series #.
[M,N] = size(posterior); 

sum_N = 0;
for i = 1:N
    sum_M = 0;
    for k = 1:M
        sum_M = sum_M + posterior(k,i).*prior(k,1);
    end
    sum_N = sum_N + log(sum_M);
end

bic = sum_N - (1/2)*(M*v+M-1)*log(N);