function cost_overall=nmf_cost(beta_loss,V,WH,H,sparsity)
%% Paper Information
% Sparse NMF--half-baked or well done?
% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.698.2772&rep=rep1&type=pdf 
%% Purpose
% Calculate cost function
%% Inputs
% beta_loss= beta loss in cost functoin
% V= actual data
% WH= estimate of data
% H= activation matrix
% sparstiy= L1 penalty on activation matrix 

if beta_loss==0
    % Itakura-Saito (IS) distance
    cost=V./WH-log(V./WH)-1; 
elseif beta_loss==1
    % Kullback-Leibler (KL) divergence
    cost=V.*(log(V)-log(WH))+(WH-V); 
else
    % if beta=2 then Euclidean distance
    cost=(1/(beta_loss*(beta_loss-1)))*(V.^beta_loss-WH.^beta_loss-beta_loss*WH.^(beta_loss-1).*(V-WH)); 
end
    
cost_overall=sum(cost,'all')+sparsity*sum(abs(H),'all');

end


