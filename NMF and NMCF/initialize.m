function [H_1_n, H_2_n, H_3_n, H_4_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n] = ...
    initialize(V, r, T, W_1, W_2, W_3, W_4, H_1_n, H_2_n, H_3_n, H_4_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n)
%% Purpose
% Initialise H and WH matrices
%% Inputs
% V= desired matrix with all recording examples
% r= contain number of bases for each components
% T= column dimension of activation matrices
% Empty matrices W_1, W_2, W_3, W_4, H_1_n, H_2_n, H_3_n, H_4_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n
%% Outputs
% Initialised matrices W_1, W_2, W_3, W_4, H_1_n, H_2_n, H_3_n, H_4_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n

for i=1:length(V)
    H_1=rand(r(1),T)+1;
    H_1_n{i}=H_1;
    WH_1_n{i}=W_1*H_1;
    
    H_2=rand(r(2),T)+1;
    H_2_n{i}=H_2;
    WH_2_n{i}=W_2*H_2;
    
    H_3=rand(r(3),T)+1;
    H_3_n{i}=H_3;
    WH_3_n{i}=W_3*H_3;
    
    H_4=rand(r(4),T)+1;
    H_4_n{i}=H_4;
    WH_4_n{i}=W_4*H_4;
end
end