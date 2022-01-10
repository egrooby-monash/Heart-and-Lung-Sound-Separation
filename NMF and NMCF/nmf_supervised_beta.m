function [W, H,cost]=nmf_supervised_beta(V,K,W,MAXITER,fixedInds,beta_loss)
%% Paper Information
% Sparse NMF--half-baked or well done?
% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.698.2772&rep=rep1&type=pdf 
%% Purpose
% Supervised sound separation with non-negative matrix factorisation with
% general cost function
%% Inputs
% V= positive matrix
% K= number of sources in mixture
% MAXITER= maximum number of iterations
% fixedInds= indices that are fixed and cannot change
% beta_loss= beta loss for cost function
%% Output
% Separated W and H matrices
% cost= cost function at each iteration

% Initialize base matrix W and activation matrix H
F = size(V,1);
T = size(V,2);

if isempty(W)
    W = 1+rand(F, sum(K));
end
H = 1+rand(sum(K), T);

inds = setdiff(1:sum(K),fixedInds);

WH=W*H+eps;
cost=zeros(MAXITER,1);
for i=1:MAXITER
    % H updates
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1);
    
    H=H.*(num./dom);
    
    WH=W*H+eps;
    
    % W updates
    num=(WH.^(beta_loss-2).*V)*H(inds,:)';
    dom=WH.^(beta_loss-1)*H(inds,:)';
    
    W(:,inds)=W(:,inds).*(num./dom);
    
    WH=W*H+eps;
    cost(i)=nmf_cost(beta_loss,V,WH,H,0); 
end
% normalize W to sum to 1
sumW = sum(W);
W = W*diag(1./sumW);
H = diag(sumW)*H;
end