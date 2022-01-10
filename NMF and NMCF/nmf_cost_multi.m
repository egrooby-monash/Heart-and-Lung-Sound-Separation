function cost_overall=nmf_cost_multi(beta_loss,V_all,H_all,W,sparsity)
%% Paper Information
% Sparse NMF--half-baked or well done?
% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.698.2772&rep=rep1&type=pdf 
%% Purpose
% Calculate cost function with multiple recordings
%% Inputs
% beta_loss= beta loss in cost functoin
% V= actual data
% WH= estimate of data
% H= activation matrix
% sparstiy= L1 penalty on activation matrix 
cost_overall=0; 
for idx=1:length(V_all)
    V=V_all{idx}+eps; 
    H=H_all{idx}; 
    WH=W*H+eps; 
    cost_overall=nmf_cost(beta_loss,V,WH,H,sparsity)+cost_overall; 
end
end


