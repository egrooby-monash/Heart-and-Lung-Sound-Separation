function [H_1_n, H_2_n, H_3_n, H_4_n, H_5_n, H_6_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n, WH_5_n, WH_6_n] = ...
    h_updates2(V_n, sparsity, beta_loss, W_1, W_2, W_3, W_4, W_5, W_6, H_1_n, H_2_n, H_3_n, H_4_n, H_5_n, H_6_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n, WH_5_n, WH_6_n)
%% Purpose
% Perform activation H updates for all 4 groups (Heart, Noise Unsupervsied,
% Noise, Lung)
%% Inputs
% V_n= desired matrix with all recording examples
% sparsity= activation L1 penalty
% beta_loss= beta loss for cost function 
% W_1, W_2, W_3, W_4, W_5, W_6= base matrices for heart, noise
% unsupervised, noise (cry),lung, stethoscope movement and respiratory support respectively
% H_1_n, H_2_n, H_3_n, H_4_n, H_5_n, H_6_n= activation matrices for eart, noise
% unsupervised, noise and lung respectively for particular group of
% recordings
% WH_1_n, WH_2_n, WH_3_n, WH_4_n, WH_5_n, WH_6_n= estimates of the recording examples V_n
%% Outputs
% Updated H_1_n, H_2_n, H_3_n, H_4_n, H_5_n, H_6_n and WH_1_n, WH_2_n, WH_3_n, WH_4_n, WH_5_n, WH_6_n
for idx=1:length(V_n)
    WH = WH_1_n{idx} + WH_2_n{idx} + WH_3_n{idx} + WH_4_n{idx}+ WH_5_n{idx}+ WH_6_n{idx} +eps;
    V=V_n{idx}+eps;
    
    % H_1_n update
    W=W_1+eps;
    H=H_1_n{idx}+eps;
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1)+sparsity;
    
    H_1_n{idx}=H.*(num./dom);
    
    WH_1_n{idx} = W_1 * H_1_n{idx};
    
    % H_2_n update
    W=W_2+eps;
    H=H_2_n{idx}+eps;
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1)+sparsity;
    
    H_2_n{idx}=H.*(num./dom);
    
    WH_2_n{idx} = W_2 * H_2_n{idx};
    
    % H_3_n update
    W=W_3+eps;
    H=H_3_n{idx}+eps;
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1)+sparsity;
    
    H_3_n{idx}=H.*(num./dom);
    
    WH_3_n{idx} = W_3 * H_3_n{idx};
    
    % H_4_n update
    W=W_4+eps;
    H=H_4_n{idx}+eps;
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1)+sparsity;
    
    H_4_n{idx}=H.*(num./dom);
    
    WH_4_n{idx} = W_4 * H_4_n{idx};  
    
    % H_5_n update
    W=W_5+eps;
    H=H_5_n{idx}+eps;
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1)+sparsity;
    
    H_5_n{idx}=H.*(num./dom);
    
    WH_5_n{idx} = W_5 * H_5_n{idx};  
    
    % H_6_n update
    W=W_6+eps;
    H=H_6_n{idx}+eps;
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1)+sparsity;
    
    H_6_n{idx}=H.*(num./dom);
    
    WH_6_n{idx} = W_6 * H_6_n{idx};  
end
end