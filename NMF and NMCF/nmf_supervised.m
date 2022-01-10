function [W, H,cost] = nmf_supervised(V, K, W, MAXITER, fixedInds)
%% Paper Information
% Single-Channel Source Separation Tutorial Mini-Series
% https://ccrma.stanford.edu/~njb/teaching/sstutorial/ 
%% Purpose
% Supervised sound separation with non-negative matrix factorisation
%% Inputs
% V= positive matrix
% K= number of sources in mixture
% MAXITER= maximum number of iterations
% fixedInds= indices that are fixed and cannot change
%% Output
% Separated W and H matrices
% cost= cost function at each iteration

% Initialize base matrix W and activation matrix H
F = size(V,1); 
T = size(V,2);

% Initialise
rand('seed',0)
if isempty(W)
    W = 1+rand(F, sum(K));
end
H = 1+rand(sum(K), T);

inds = setdiff(1:sum(K),fixedInds);
ONES = ones(F,T);
cost=zeros(MAXITER,1);
for i=1:MAXITER 
    % update activations
    H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES);
    
    % update dictionaries
    W(:,inds) = W(:,inds) .* ((V./(W*H+eps))*H(inds,:)') ./(ONES*H(inds,:)');
    cost(i)=nmf_cost(1,V,W*H,H,0); 
end

% normalize W to sum to 1
sumW = sum(W);
W = W*diag(1./sumW);
H = diag(sumW)*H;

end


