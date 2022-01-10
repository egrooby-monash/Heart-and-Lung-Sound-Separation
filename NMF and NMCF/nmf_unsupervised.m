function [W, H,cost] = nmf_unsupervised(V, K, MAXITER)
%% Paper Information
% Single-Channel Source Separation Tutorial Mini-Series
% https://ccrma.stanford.edu/~njb/teaching/sstutorial/ 
%% Purpose
% Single source sound separation wiht non-negative matrix factorisation
%% Inputs
% V= positive matrix
% K= number of sources in mixture
% MAXITER= maximum number of iterations
%% Output
% Separated W and H matrices
% cost= cost function at each iteration

F = size(V,1);
T = size(V,2);

% Initialise
rand('seed',0)
W = 1+rand(F, K);
% W = W./repmat(sum(W),F,1);
H = 1+rand(K, T);

ONES = ones(F,T);
cost=zeros(MAXITER,1);
for i=1:MAXITER 
    % update activations
    H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES);
    
    % update dictionaries
    W = W .* ((V./(W*H+eps))*H') ./(ONES*H');
    
    cost(i)=nmf_cost(1,V,W*H,H,0); 
end

% normalize W to sum to 1
sumW = sum(W);
W = W*diag(1./sumW);
H = diag(sumW)*H;

end