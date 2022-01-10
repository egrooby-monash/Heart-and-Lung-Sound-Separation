function [W,H_all,cost]=nmf_supervised_sparse_multi(V_all,K,W,MAXITER,fixedInds,beta_loss,sparsity)
%% Paper Information
% Sparse NMF--half-baked or well done?
% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.698.2772&rep=rep1&type=pdf 
%% Purpose
% Supervised sound separation with non-negative matrix factorisation with
% general cost function and sparsity for multiple recordings
%% Inputs
% V= positive matrix
% K= number of sources in mixture
% MAXITER= maximum number of iterations
% fixedInds= indices that are fixed and cannot change
% beta_loss= beta loss for cost function
% sparsity= activation L1 penalty
%% Output
% Separated W and H matrices
% cost= cost function at each iteration

% Initialize base matrix W and activation matrix H
F = size(V_all{1},1);
T = size(V_all{1},2);

if isempty(W)
    W = 1+rand(F, sum(K));
end

num=length(V_all);
H_all = cell(1,num);
for i=1:length(V_all)
    H_all{i} = 1+rand(sum(K), T);
end

inds = setdiff(1:sum(K),fixedInds);

[W]=matrix_normalization(W);

cost=zeros(MAXITER,1);
for i=1:MAXITER
    for idx=1:length(V_all)
    H=H_all{idx}+eps;
    V=V_all{idx}+eps;    
    WH=W*H+eps;
    % H updates
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1)+sparsity;
    H_all{idx}=H.*(num./dom);
    end
    
    negdp=zeros(size(W(:,inds)));
    posdp=zeros(size(W(:,inds)));
    for idx=1:length(V_all)
        H=H_all{idx}+eps;
        V=V_all{idx}+eps;
        WH=W*H+eps;
        
        % W updates
        temp=(WH.^(beta_loss-1)*H(inds,:)').*W(:,inds);
        temp=ones(F)*temp;
        temp=temp.*W(:,inds);
        negdp=(WH.^(beta_loss-2).*V)*H(inds,:)'+temp+negdp;

        temp=((WH.^(beta_loss-2).*V)*H(inds,:)').*W(:,inds);
        temp=ones(F)*temp;
        temp=temp.*W(:,inds);
        posdp=WH.^(beta_loss-1)*H(inds,:)'+temp+posdp;
    end
    W(:,inds)=W(:,inds).*(negdp./posdp);
    [W]=matrix_normalization(W);
    cost(i)=nmf_cost_multi(beta_loss,V_all,H_all,W,sparsity);
end

end
