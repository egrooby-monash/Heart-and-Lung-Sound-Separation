function [W, H] = nmf(V, K, MAXITER)

F = size(V,1);
T = size(V,2);

rand('seed',0)
W = 1+rand(F, K);
% W = W./repmat(sum(W),F,1);
H = 1+rand(K, T);

ONES = ones(F,T);

for i=1:MAXITER 
    % update activations
    H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES);
    
    % update dictionaries
    W = W .* ((V./(W*H+eps))*H') ./(ONES*H');

end


% normalize W to sum to 1
sumW = sum(W);
W = W*diag(1./sumW);
H = diag(sumW)*H;