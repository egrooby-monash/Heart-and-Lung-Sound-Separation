function [W,H,cost]=nmf_supervised_sparse(V,K,W,MAXITER,fixedInds,beta_loss,sparsity)
%% Paper Information
% Sparse NMF--half-baked or well done?
% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.698.2772&rep=rep1&type=pdf 
%% Purpose
% Supervised sound separation with non-negative matrix factorisation with
% general cost function and sparsity
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

[W]=matrix_normalization(W);
WH=W*H+eps;
cost=zeros(MAXITER,1);
for i=1:MAXITER
    %H updates
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1)+sparsity;
    
    H=H.*(num./dom);
    
    WH=W*H+eps;
    
    % W updates
    temp=(WH.^(beta_loss-1)*H(inds,:)').*W(:,inds);
    temp=ones(F)*temp;
    temp=temp.*W(:,inds);
    negdp=(WH.^(beta_loss-2).*V)*H(inds,:)'+temp;
    
    temp=((WH.^(beta_loss-2).*V)*H(inds,:)').*W(:,inds);
    temp=ones(F)*temp;
    temp=temp.*W(:,inds);
    posdp=WH.^(beta_loss-1)*H(inds,:)'+temp;
    
    W(:,inds)=W(:,inds).*(negdp./posdp);
    
    [W]=matrix_normalization(W);
    WH=W*H+eps;
    cost(i)=nmf_cost(beta_loss,V,WH,H,sparsity); 
end

end
