function [negdp, posdp]=...
    w_updates(V_n, weight, F, negdp, posdp, beta_loss, W_1, H_1_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n)
%% Purpose
% Perform base W matrix updates
%% Inputs
%% Inputs
% V_n= desired matrix with all recording examples
% negdp and posdp= current update numerator and denominator 
% w_n = weight
% beta_loss= beta loss for cost function 
% W_1, W_2, W_3, W_4= base matrices for heart, noise unsupervised, noise
% and lung respectively
% WH_1_n, WH_2_n, WH_3_n, WH_4_n= estimates of the recording examples V_n
%% Outputs
% negdp and posdp= Uupdate numerator and denominator 

for idx=1:length(V_n)
    WH = WH_1_n{idx} + WH_2_n{idx} + WH_3_n{idx} + WH_4_n{idx}+eps;
    V=V_n{idx}+eps;
    H=H_1_n{idx}+eps;
    W=W_1+eps;
    temp=(WH.^(beta_loss-1)*H').*W;
    temp=ones(F)*temp;
    temp=temp.*W;
    negdp=weight*((WH.^(beta_loss-2).*V)*H'+temp) +negdp;
    
    temp=((WH.^(beta_loss-2).*V)*H').*W;
    temp=ones(F)*temp;
    temp=temp.*W;
    posdp=weight*(WH.^(beta_loss-1)*H'+temp) + posdp;
    
end
end




