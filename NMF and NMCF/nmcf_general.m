function [W_1, W_2, W_3, W_4, H_1_m, H_2_m, H_3_m, H_4_m,cost] = ...
    nmcf_general(V_h, V_l, V_n, V_m, r, options_nmf, MAXITER)
%% Purpose
% Perform NMCF with the following setup
% Heart Supervised
% Lung Supervsied
% Noise Supervised
% Noise Unsupervised
%% Inputs
% V_h and w_h are heart segments and weights
% V_l and w_l are lung segments and weights
% V_n and w_n are noise segments and weights
% V_m is mixture segments
% r is inner dimension set with length 4
% max_iter is maximum number of iterations
% beta_loss= beta loss for cost function 
% sparsity= activation L1 penalty
%% Outputs
% W_1 heart, W_2 unsupervised, W_3 noise, W_4 lung base matrices
% And their corresponding activations for the mixture recordings 

if nargin<1
    V_h=[];
    V_l=[];
    V_n=[];
    V_m=[];
    r=[];
    options_nmf.sparsity=0;
    options_nmf.beta_loss=1;
    MAXITER=200;
end



% Segment dimensions
F=size(V_h{1},1);
T=size(V_h{1},2);

% Initialise base matrix W
% W_1 heart, W_2 unsupervised, W_3 noise, W_4 lung
W_1=1+rand(F,r(1));
W_2=1+rand(F,r(2));
W_3=1+rand(F,r(3));
W_4=1+rand(F,r(4));

W_1 = matrix_normalization(W_1);
W_2 = matrix_normalization(W_2);
W_3 = matrix_normalization(W_3);
W_4 = matrix_normalization(W_4);

% Initialise activation matrix H
% convention: 1,2,3,4 refers to separation and m,h,n refers to of what dataset

% mixture group
num=length(V_m);
H_1_m = cell(1,num);
H_2_m = cell(1,num);
H_3_m = cell(1,num);
H_4_m = cell(1,num);
WH_1_m = cell(1,num);
WH_2_m = cell(1,num);
WH_3_m = cell(1,num);
WH_4_m = cell(1,num);

% heart group
num=length(V_h);
H_1_h = cell(1,num);
H_2_h = cell(1,num);
H_3_h = cell(1,num);
H_4_h = cell(1,num);
WH_1_h = cell(1,num);
WH_2_h = cell(1,num);
WH_3_h = cell(1,num);
WH_4_h = cell(1,num);

% noise group
num=length(V_n);
H_1_n = cell(1,num);
H_2_n = cell(1,num);
H_3_n = cell(1,num);
H_4_n = cell(1,num);
WH_1_n = cell(1,num);
WH_2_n = cell(1,num);
WH_3_n = cell(1,num);
WH_4_n = cell(1,num);

% lung group
num=length(V_l);
H_1_l = cell(1,num);
H_2_l = cell(1,num);
H_3_l = cell(1,num);
H_4_l = cell(1,num);
WH_1_l = cell(1,num);
WH_2_l = cell(1,num);
WH_3_l = cell(1,num);
WH_4_l = cell(1,num);

%initialization of variables
[H_1_n, H_2_n, H_3_n, H_4_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n] = ...
    initialize(V_n, r, T, W_1, W_2, W_3, W_4, H_1_n, H_2_n, H_3_n, H_4_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n);

[H_1_m, H_2_m, H_3_m, H_4_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m] = ...
    initialize(V_m, r, T, W_1, W_2, W_3, W_4, H_1_m, H_2_m, H_3_m, H_4_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m);

[H_1_h, H_2_h, H_3_h, H_4_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h] = ...
    initialize(V_h, r, T, W_1, W_2, W_3, W_4, H_1_h, H_2_h, H_3_h, H_4_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h);

[H_1_l, H_2_l, H_3_l, H_4_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l] = ...
    initialize(V_l, r, T, W_1, W_2, W_3, W_4, H_1_l, H_2_l, H_3_l, H_4_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l);

% Update
cost=zeros(MAXITER,1);
for i=1:MAXITER
    % H_n updates
    [H_1_n, H_2_n, H_3_n, H_4_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n] = ...
        h_updates(V_n, options_nmf.sparsity, options_nmf.beta_loss, W_1, W_2, W_3, W_4, H_1_n, H_2_n, H_3_n, H_4_n, WH_1_n, WH_2_n, WH_3_n,WH_4_n);
    
    % H_m updates
    [H_1_m, H_2_m, H_3_m, H_4_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m] = ...
        h_updates(V_m, options_nmf.sparsity, options_nmf.beta_loss, W_1, W_2, W_3, W_4, H_1_m, H_2_m, H_3_m, H_4_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m);
    
    % H_h updates
    [H_1_h, H_2_h, H_3_h, H_4_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h] = ...
        h_updates(V_h, options_nmf.sparsity, options_nmf.beta_loss, W_1, W_2, W_3, W_4, H_1_h, H_2_h, H_3_h, H_4_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h);
    
    % H_l updates
    [H_1_l, H_2_l, H_3_l, H_4_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l] = ...
        h_updates(V_l, options_nmf.sparsity, options_nmf.beta_loss, W_1, W_2, W_3, W_4, H_1_l, H_2_l, H_3_l, H_4_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l);
    
    % W_1 update
    negdp=zeros(F,r(1));
    posdp=zeros(F,r(1));
%     [negdp, posdp] = ...
%         w_updates(V_n, w_n, 0, weight_i, F, negdp, posdp, beta_loss, W_1, H_1_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n);
    [negdp, posdp] = ...
        w_updates(V_m, 1, F, negdp, posdp, options_nmf.beta_loss, W_1, H_1_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m);
    [negdp, posdp] = ...
        w_updates(V_h, 1/length(V_h), F, negdp, posdp, options_nmf.beta_loss, W_1, H_1_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h);
%     [negdp, posdp] = ...
%         w_updates(V_l, w_l, 0, weight_i, F, negdp, posdp, beta_loss, W_1, H_1_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l);
    
    W_1=W_1.*(negdp./posdp);
    W_1=matrix_normalization(W_1);
    for idx= 1:length(H_1_n)
        WH_1_n{idx} = W_1 * H_1_n{idx}+eps;
    end
    for idx= 1:length(H_1_m)
        WH_1_m{idx} = W_1 * H_1_m{idx}+eps;
    end
    for idx= 1:length(H_1_h)
        WH_1_h{idx} = W_1 * H_1_h{idx}+eps;
    end
    for idx= 1:length(H_1_l)
        WH_1_l{idx} = W_1 * H_1_l{idx}+eps;
    end
    
    % W_2 update
    negdp=zeros(F,r(2));
    posdp=zeros(F,r(2));
      
%     [negdp, posdp] = ...
%         w_updates(V_n, w_n, 0, weight_i, F, negdp, posdp, beta_loss, W_2, H_2_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n);
    [negdp, posdp] = ...
        w_updates(V_m, 1, F, negdp, posdp, options_nmf.beta_loss, W_2, H_2_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m);
%     [negdp, posdp] = ...
%         w_updates(V_h, w_h, 0, weight_i, F, negdp, posdp, beta_loss, W_2, H_2_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h);
%     [negdp, posdp] = ...
%         w_updates(V_l, w_l, 0, weight_i, F, negdp, posdp, beta_loss, W_2, H_2_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l);
    
    W_2=W_2.*(negdp./posdp);
    W_2=matrix_normalization(W_2);
    
    for idx= 1:length(H_2_n)
        WH_2_n{idx} = W_2 * H_2_n{idx}+eps;
    end
    for idx= 1:length(H_2_m)
        WH_2_m{idx} = W_2 * H_2_m{idx}+eps;
    end
    for idx= 1:length(H_2_h)
        WH_2_h{idx} = W_2 * H_2_h{idx}+eps;
    end
    for idx= 1:length(H_2_l)
        WH_2_l{idx} = W_2 * H_2_l{idx}+eps;
    end
    
    %W_3 update
    negdp=zeros(F,r(3));
    posdp=zeros(F,r(3));
    [negdp, posdp] = ...
        w_updates(V_n, 1/length(V_n), F, negdp, posdp, options_nmf.beta_loss, W_3, H_3_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n);
    [negdp, posdp] = ...
        w_updates(V_m, 1, F, negdp, posdp, options_nmf.beta_loss, W_3, H_3_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m);
%     [negdp, posdp] = ...
%         w_updates(V_h, w_h, 0, weight_i, F, negdp, posdp, beta_loss, W_3, H_3_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h);
%     [negdp, posdp] = ...
%         w_updates(V_l, w_l, 0, weight_i, F, negdp, posdp, beta_loss, W_3, H_3_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l);
    
    W_3=W_3.*(negdp./posdp);
    W_3=matrix_normalization(W_3);
    
    for idx= 1:length(H_3_n)
        WH_3_n{idx} = W_3 * H_3_n{idx}+eps;
    end
    for idx= 1:length(H_3_m)
        WH_3_m{idx} = W_3 * H_3_m{idx}+eps;
    end
    for idx= 1:length(H_3_h)
        WH_3_h{idx} = W_3 * H_3_h{idx}+eps;
    end
    for idx= 1:length(H_3_l)
        WH_3_l{idx} = W_3 * H_3_l{idx}+eps;
    end
    
    % W_4 update
    negdp=zeros(F,r(4));
    posdp=zeros(F,r(4));
    
%     [negdp, posdp] = ...
%         w_updates(V_n, w_n, 0, weight_i, F, negdp, posdp, beta_loss, W_4, H_4_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n);
    [negdp, posdp] = ...
        w_updates(V_m, 1, F, negdp, posdp, options_nmf.beta_loss, W_4, H_4_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m);
%     [negdp, posdp] = ...
%         w_updates(V_h, w_h, 0, weight_i, F, negdp, posdp, beta_loss, W_4, H_4_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h);
    [negdp, posdp] = ...
        w_updates(V_l, 1/length(V_l), F, negdp, posdp, options_nmf.beta_loss, W_4, H_4_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l);
    
    W_4=W_4.*(negdp./posdp);
    W_4=matrix_normalization(W_4);
    
    for idx= 1:length(H_4_n)
        WH_4_n{idx} = W_4 * H_4_n{idx}+eps;
    end
    for idx= 1:length(H_4_m)
        WH_4_m{idx} = W_4 * H_4_m{idx}+eps;
    end
    for idx= 1:length(H_4_h)
        WH_4_h{idx} = W_4 * H_4_h{idx}+eps;
    end
    for idx= 1:length(H_4_l)
        WH_4_l{idx} = W_4 * H_4_l{idx}+eps;
    end
    
    
    cost_overall=0; 
    for idx=1:length(V_l)    
        V=V_l{idx}+eps;
        H=H_4_l{idx}+eps;
        %W=W_4+eps; 
        WH=WH_1_l{idx} + WH_2_l{idx} + WH_3_l{idx} + WH_4_l{idx}+eps;%W*H+eps;
        cost_overall=(1/length(V_l))*nmf_cost(options_nmf.beta_loss,V,WH,H,options_nmf.sparsity)+cost_overall;
    end
    
    for idx=1:length(V_h)
        V=V_h{idx}+eps;
        H=H_1_h{idx};
        %W=W_1;
        WH=WH_1_h{idx} + WH_2_h{idx} + WH_3_h{idx} + WH_4_h{idx}+eps;%WH=W*H+eps;
        cost_overall=(1/length(V_h))*nmf_cost(options_nmf.beta_loss,V,WH,H,options_nmf.sparsity)+cost_overall;
    end
    
    for idx=1:length(V_n)
        V=V_n{idx}+eps;
        H=H_3_n{idx}; 
        %W=W_3;
        WH=WH_1_n{idx} + WH_2_n{idx} + WH_3_n{idx} + WH_4_n{idx}+eps;%WH=W*H+eps;
        cost_overall=(1/length(V_n))*nmf_cost(options_nmf.beta_loss,V,WH,H,options_nmf.sparsity)+cost_overall;
    end
    
    for idx=1:length(V_m)
        V=V_m{idx}+eps;
        H=[H_1_m{idx};H_2_m{idx};H_3_m{idx};H_4_m{idx}];
        W=[W_1,W_2,W_3,W_4];
        WH=W*H+eps;
        cost_overall=nmf_cost(options_nmf.beta_loss,V,WH,H,options_nmf.sparsity)+cost_overall;
    end
    
    cost(i)=cost_overall;
end
end











