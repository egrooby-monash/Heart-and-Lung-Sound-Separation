function [W_1, W_2, W_3, W_4, W_5, W_6, H_1_m, H_2_m, H_3_m, H_4_m, H_5_m, H_6_m, cost] = ...
    nmcf_general2(V_h, V_l, V_n, V_s, V_r, V_m, r, options_nmf, MAXITER)
%% Purpose
% Perform NMCF with the following setup
% Heart Supervised
% Lung Supervsied
% Cry Noise, Stethoscope Movement Noise and Respiratory Support Noise Supervised
% Noise Unsupervised
%% Inputs
% V_h is are heart segments
% V_l is lung segments
% V_n is noise segments (cry)
% V_s is stethoscope movement segments
% V_r is respiratory support segments
% V_m is mixture segments
% r is inner dimension set with length 6
% max_iter is maximum number of iterations
% beta_loss= beta loss for cost function
% sparsity= activation L1 penalty
%% Outputs
% W_1 heart, W_2 unsupervised, W_3 cry, W_4 lung W_5 stethoscope movement W_6 respiratory
% support base matrices
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
    options_nmf.W1=1;
    options_nmf.W3=1;
    options_nmf.W4=1;
    options_nmf.W5=1;
    options_nmf.W6=1;
end



% Segment dimensions
F=size(V_h{1},1);
T=size(V_h{1},2);

% Initialise base matrix W
% W_1 heart, W_2 unsupervised, W_3 noise, W_4 lung, W_5 stethoscope
% movement, W_6 respiratory support
W_1=1+rand(F,r(1));
W_2=1+rand(F,r(2));
W_3=1+rand(F,r(3));
W_4=1+rand(F,r(4));
W_5=1+rand(F,r(5));
W_6=1+rand(F,r(6));

W_1 = matrix_normalization(W_1);
W_2 = matrix_normalization(W_2);
W_3 = matrix_normalization(W_3);
W_4 = matrix_normalization(W_4);
W_5 = matrix_normalization(W_5);
W_6 = matrix_normalization(W_6);

% Initialise activation matrix H
% convention: 1,2,3,4,5,6 refers to separation and m,h,n refers to of what dataset

% mixture group
num=length(V_m);
H_1_m = cell(1,num);
H_2_m = cell(1,num);
H_3_m = cell(1,num);
H_4_m = cell(1,num);
H_5_m = cell(1,num);
H_6_m = cell(1,num);
WH_1_m = cell(1,num);
WH_2_m = cell(1,num);
WH_3_m = cell(1,num);
WH_4_m = cell(1,num);
WH_5_m = cell(1,num);
WH_6_m = cell(1,num);

% heart group
num=length(V_h);
H_1_h = cell(1,num);
H_2_h = cell(1,num);
H_3_h = cell(1,num);
H_4_h = cell(1,num);
H_5_h = cell(1,num);
H_6_h = cell(1,num);
WH_1_h = cell(1,num);
WH_2_h = cell(1,num);
WH_3_h = cell(1,num);
WH_4_h = cell(1,num);
WH_5_h = cell(1,num);
WH_6_h = cell(1,num);

% noise group
num=length(V_n);
H_1_n = cell(1,num);
H_2_n = cell(1,num);
H_3_n = cell(1,num);
H_4_n = cell(1,num);
H_5_n = cell(1,num);
H_6_n = cell(1,num);
WH_1_n = cell(1,num);
WH_2_n = cell(1,num);
WH_3_n = cell(1,num);
WH_4_n = cell(1,num);
WH_5_n = cell(1,num);
WH_6_n = cell(1,num);

% lung group
num=length(V_l);
H_1_l = cell(1,num);
H_2_l = cell(1,num);
H_3_l = cell(1,num);
H_4_l = cell(1,num);
H_5_l = cell(1,num);
H_6_l = cell(1,num);
WH_1_l = cell(1,num);
WH_2_l = cell(1,num);
WH_3_l = cell(1,num);
WH_4_l = cell(1,num);
WH_5_l = cell(1,num);
WH_6_l = cell(1,num);

% stmv group
num=length(V_s);
H_1_s = cell(1,num);
H_2_s = cell(1,num);
H_3_s = cell(1,num);
H_4_s = cell(1,num);
H_5_s = cell(1,num);
H_6_s = cell(1,num);
WH_1_s = cell(1,num);
WH_2_s = cell(1,num);
WH_3_s = cell(1,num);
WH_4_s = cell(1,num);
WH_5_s = cell(1,num);
WH_6_s = cell(1,num);

% respiratory group
num=length(V_r);
H_1_r = cell(1,num);
H_2_r = cell(1,num);
H_3_r = cell(1,num);
H_4_r = cell(1,num);
H_5_r = cell(1,num);
H_6_r = cell(1,num);
WH_1_r = cell(1,num);
WH_2_r = cell(1,num);
WH_3_r = cell(1,num);
WH_4_r = cell(1,num);
WH_5_r = cell(1,num);
WH_6_r = cell(1,num);

%initialization of variables
[H_1_n, H_2_n, H_3_n, H_4_n, H_5_n, H_6_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n, WH_5_n, WH_6_n] = ...
    initialize2(V_n, r, T, W_1, W_2, W_3, W_4, W_5, W_6, H_1_n, H_2_n, H_3_n, H_4_n, H_5_n, H_6_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n, WH_5_n, WH_6_n);

[H_1_m, H_2_m, H_3_m, H_4_m, H_5_m, H_6_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m, WH_5_m, WH_6_m] = ...
    initialize2(V_m, r, T, W_1, W_2, W_3, W_4, W_5, W_6, H_1_m, H_2_m, H_3_m, H_4_m, H_5_m, H_6_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m, WH_5_m, WH_6_m);

[H_1_h, H_2_h, H_3_h, H_4_h, H_5_h, H_6_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h, WH_5_h, WH_6_h] = ...
    initialize2(V_h, r, T, W_1, W_2, W_3, W_4, W_5, W_6, H_1_h, H_2_h, H_3_h, H_4_h, H_5_h, H_6_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h, WH_5_h, WH_6_h);

[H_1_l, H_2_l, H_3_l, H_4_l, H_5_l, H_6_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l, WH_5_l, WH_6_l] = ...
    initialize2(V_l, r, T, W_1, W_2, W_3, W_4, W_5, W_6, H_1_l, H_2_l, H_3_l, H_4_l, H_5_l, H_6_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l, WH_5_l, WH_6_l);

[H_1_s, H_2_s, H_3_s, H_4_s, H_5_s, H_6_s, WH_1_s, WH_2_s, WH_3_s, WH_4_s, WH_5_s, WH_6_s] = ...
    initialize2(V_s, r, T, W_1, W_2, W_3, W_4, W_5, W_6, H_1_s, H_2_s, H_3_s, H_4_s, H_5_s, H_6_s, WH_1_s, WH_2_s, WH_3_s, WH_4_s, WH_5_s, WH_6_s);

[H_1_r, H_2_r, H_3_r, H_4_r, H_5_r, H_6_r, WH_1_r, WH_2_r, WH_3_r, WH_4_r, WH_5_r, WH_6_r] = ...
    initialize2(V_r, r, T, W_1, W_2, W_3, W_4, W_5, W_6, H_1_r, H_2_r, H_3_r, H_4_r, H_5_r, H_6_r, WH_1_r, WH_2_r, WH_3_r, WH_4_r, WH_5_r, WH_6_r);

% Update
cost=zeros(MAXITER,1);
for i=1:MAXITER
    % H_n updates
    [H_1_n, H_2_n, H_3_n, H_4_n, H_5_n, H_6_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n, WH_5_n, WH_6_n] = ...
        h_updates2(V_n, options_nmf.sparsity, options_nmf.beta_loss, W_1, W_2, W_3, W_4, W_5, W_6, H_1_n, H_2_n, H_3_n, H_4_n, H_5_n, H_6_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n, WH_5_n, WH_6_n);
    
    % H_m updates
    [H_1_m, H_2_m, H_3_m, H_4_m, H_5_m, H_6_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m, WH_5_m, WH_6_m] = ...
        h_updates2(V_m, options_nmf.sparsity, options_nmf.beta_loss, W_1, W_2, W_3, W_4, W_5, W_6, H_1_m, H_2_m, H_3_m, H_4_m, H_5_m, H_6_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m, WH_5_m, WH_6_m);
    
    % H_h updates
    [H_1_h, H_2_h, H_3_h, H_4_h, H_5_h, H_6_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h, WH_5_h, WH_6_h] = ...
        h_updates2(V_h, options_nmf.sparsity, options_nmf.beta_loss, W_1, W_2, W_3, W_4, W_5, W_6, H_1_h, H_2_h, H_3_h, H_4_h, H_5_h, H_6_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h, WH_5_h, WH_6_h);
    
    % H_l updates
    [H_1_l, H_2_l, H_3_l, H_4_l, H_5_l, H_6_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l, WH_5_l, WH_6_l] = ...
        h_updates2(V_l, options_nmf.sparsity, options_nmf.beta_loss, W_1, W_2, W_3, W_4, W_5, W_6, H_1_l, H_2_l, H_3_l, H_4_l, H_5_l, H_6_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l, WH_5_l, WH_6_l);
    
    % H_s updates
    [H_1_s, H_2_s, H_3_s, H_4_s, H_5_s, H_6_s, WH_1_s, WH_2_s, WH_3_s, WH_4_s, WH_5_s, WH_6_s] = ...
        h_updates2(V_s, options_nmf.sparsity, options_nmf.beta_loss, W_1, W_2, W_3, W_4, W_5, W_6, H_1_s, H_2_s, H_3_s, H_4_s, H_5_s, H_6_s, WH_1_s, WH_2_s, WH_3_s, WH_4_s, WH_5_s, WH_6_s);
    
    % H_r updates
    [H_1_r, H_2_r, H_3_r, H_4_r, H_5_r, H_6_r, WH_1_r, WH_2_r, WH_3_r, WH_4_r, WH_5_r, WH_6_r] = ...
        h_updates2(V_r, options_nmf.sparsity, options_nmf.beta_loss, W_1, W_2, W_3, W_4, W_5, W_6, H_1_r, H_2_r, H_3_r, H_4_r, H_5_r, H_6_r, WH_1_r, WH_2_r, WH_3_r, WH_4_r, WH_5_r, WH_6_r);
    
    % W_1 update
    negdp=zeros(F,r(1));
    posdp=zeros(F,r(1));
    if options_nmf.W1~=0
        [negdp, posdp] = ...
            w_updates2(V_m, options_nmf.W1, F, negdp, posdp, options_nmf.beta_loss, W_1, H_1_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m, WH_5_m, WH_6_m);
    end
    [negdp, posdp] = ...
        w_updates2(V_h, 1/length(V_h), F, negdp, posdp, options_nmf.beta_loss, W_1, H_1_h, WH_1_h, WH_2_h, WH_3_h, WH_4_h, WH_5_h, WH_6_h);
    
    W_1=W_1.*(negdp./posdp);
    W_1=matrix_normalization(W_1);
    for idx= 1:length(H_1_n)
        WH_1_n{idx} = W_1 * H_1_n{idx};
    end
    for idx= 1:length(H_1_m)
        WH_1_m{idx} = W_1 * H_1_m{idx};
    end
    for idx= 1:length(H_1_h)
        WH_1_h{idx} = W_1 * H_1_h{idx};
    end
    for idx= 1:length(H_1_l)
        WH_1_l{idx} = W_1 * H_1_l{idx};
    end
    for idx= 1:length(H_1_s)
        WH_1_s{idx} = W_1 * H_1_s{idx};
    end
    for idx= 1:length(H_1_r)
        WH_1_r{idx} = W_1 * H_1_r{idx};
    end
    
    % W_2 update
    negdp=zeros(F,r(2));
    posdp=zeros(F,r(2));
    [negdp, posdp] = ...
        w_updates2(V_m, 1, F, negdp, posdp, options_nmf.beta_loss, W_2, H_2_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m, WH_5_m, WH_6_m);
    
    W_2=W_2.*(negdp./posdp);
    W_2=matrix_normalization(W_2);
    
    for idx= 1:length(H_2_n)
        WH_2_n{idx} = W_2 * H_2_n{idx};
    end
    for idx= 1:length(H_2_m)
        WH_2_m{idx} = W_2 * H_2_m{idx};
    end
    for idx= 1:length(H_2_h)
        WH_2_h{idx} = W_2 * H_2_h{idx};
    end
    for idx= 1:length(H_2_l)
        WH_2_l{idx} = W_2 * H_2_l{idx};
    end
    for idx= 1:length(H_2_s)
        WH_2_s{idx} = W_2 * H_2_s{idx};
    end
    for idx= 1:length(H_2_r)
        WH_2_r{idx} = W_2 * H_2_r{idx};
    end
    
    
    %W_3 update
    negdp=zeros(F,r(3));
    posdp=zeros(F,r(3));
    [negdp, posdp] = ...
        w_updates2(V_n, 1/length(V_n), F, negdp, posdp, options_nmf.beta_loss, W_3, H_3_n, WH_1_n, WH_2_n, WH_3_n, WH_4_n, WH_5_n, WH_6_n);
    if options_nmf.W3~=0
        [negdp, posdp] = ...
            w_updates2(V_m, options_nmf.W3, F, negdp, posdp, options_nmf.beta_loss, W_3, H_3_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m, WH_5_m, WH_6_m);
    end
    W_3=W_3.*(negdp./posdp);
    W_3=matrix_normalization(W_3);
    
    for idx= 1:length(H_3_n)
        WH_3_n{idx} = W_3 * H_3_n{idx};
    end
    for idx= 1:length(H_3_m)
        WH_3_m{idx} = W_3 * H_3_m{idx};
    end
    for idx= 1:length(H_3_h)
        WH_3_h{idx} = W_3 * H_3_h{idx};
    end
    for idx= 1:length(H_3_l)
        WH_3_l{idx} = W_3 * H_3_l{idx};
    end
    for idx= 1:length(H_3_s)
        WH_3_s{idx} = W_3 * H_3_s{idx};
    end
    for idx= 1:length(H_3_r)
        WH_3_r{idx} = W_3 * H_3_r{idx};
    end
    
    % W_4 update
    negdp=zeros(F,r(4));
    posdp=zeros(F,r(4));
    if options_nmf.W4~=0
        [negdp, posdp] = ...
            w_updates2(V_m, options_nmf.W4, F, negdp, posdp, options_nmf.beta_loss, W_4, H_4_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m, WH_5_m, WH_6_m);
    end
    [negdp, posdp] = ...
        w_updates2(V_l, 1/length(V_l), F, negdp, posdp, options_nmf.beta_loss, W_4, H_4_l, WH_1_l, WH_2_l, WH_3_l, WH_4_l, WH_5_l, WH_6_l);
    
    W_4=W_4.*(negdp./posdp);
    W_4=matrix_normalization(W_4);
    
    for idx= 1:length(H_4_n)
        WH_4_n{idx} = W_4 * H_4_n{idx};
    end
    for idx= 1:length(H_4_m)
        WH_4_m{idx} = W_4 * H_4_m{idx};
    end
    for idx= 1:length(H_4_h)
        WH_4_h{idx} = W_4 * H_4_h{idx};
    end
    for idx= 1:length(H_4_l)
        WH_4_l{idx} = W_4 * H_4_l{idx};
    end
    for idx= 1:length(H_4_s)
        WH_4_s{idx} = W_4 * H_4_s{idx};
    end
    for idx= 1:length(H_4_r)
        WH_4_r{idx} = W_4 * H_4_r{idx};
    end
    
    % W_5 update
    negdp=zeros(F,r(5));
    posdp=zeros(F,r(5));
    if options_nmf.W5~=0
        [negdp, posdp] = ...
            w_updates2(V_m, options_nmf.W5, F, negdp, posdp, options_nmf.beta_loss, W_5, H_5_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m, WH_5_m, WH_6_m);
    end
    [negdp, posdp] = ...
        w_updates2(V_s, 1/length(V_s), F, negdp, posdp, options_nmf.beta_loss, W_5, H_5_s, WH_1_s, WH_2_s, WH_3_s, WH_4_s, WH_5_s, WH_6_s);
    
    W_5=W_5.*(negdp./posdp);
    W_5=matrix_normalization(W_5);
    
    for idx= 1:length(H_5_n)
        WH_5_n{idx} = W_5 * H_5_n{idx};
    end
    for idx= 1:length(H_5_m)
        WH_5_m{idx} = W_5 * H_5_m{idx};
    end
    for idx= 1:length(H_5_h)
        WH_5_h{idx} = W_5 * H_5_h{idx};
    end
    for idx= 1:length(H_5_l)
        WH_5_l{idx} = W_5 * H_5_l{idx};
    end
    for idx= 1:length(H_5_s)
        WH_5_s{idx} = W_5 * H_5_s{idx};
    end
    for idx= 1:length(H_5_r)
        WH_5_r{idx} = W_5 * H_5_r{idx};
    end
    
    % W_6 update
    negdp=zeros(F,r(6));
    posdp=zeros(F,r(6));
    if options_nmf.W6~=0
        [negdp, posdp] = ...
            w_updates2(V_m, options_nmf.W6, F, negdp, posdp, options_nmf.beta_loss, W_6, H_6_m, WH_1_m, WH_2_m, WH_3_m, WH_4_m, WH_5_m, WH_6_m);
    end
    [negdp, posdp] = ...
        w_updates2(V_r, 1/length(V_r), F, negdp, posdp, options_nmf.beta_loss, W_6, H_6_r, WH_1_r, WH_2_r, WH_3_r, WH_4_r, WH_5_r, WH_6_r);
    
    W_6=W_6.*(negdp./posdp);
    W_6=matrix_normalization(W_6);
    
    for idx= 1:length(H_6_n)
        WH_6_n{idx} = W_6 * H_6_n{idx};
    end
    for idx= 1:length(H_6_m)
        WH_6_m{idx} = W_6 * H_6_m{idx};
    end
    for idx= 1:length(H_6_h)
        WH_6_h{idx} = W_6 * H_6_h{idx};
    end
    for idx= 1:length(H_6_l)
        WH_6_l{idx} = W_6 * H_6_l{idx};
    end
    for idx= 1:length(H_6_s)
        WH_6_s{idx} = W_6 * H_6_s{idx};
    end
    for idx= 1:length(H_6_r)
        WH_6_r{idx} = W_6 * H_6_r{idx};
    end
    
    cost=0;
%     cost_overall=0;
%     for idx=1:length(V_l)
%         V=V_l{idx}+eps;
%         H=H_4_l{idx}+eps;
%         WH=WH_1_l{idx} + WH_2_l{idx} + WH_3_l{idx} + WH_4_l{idx}+ WH_5_l{idx}+ WH_6_l{idx}+eps;
%         cost_overall=(options_nmf.W4/length(V_l))*nmf_cost(options_nmf.beta_loss,V,WH,H,options_nmf.sparsity)+cost_overall;
%     end
%     
%     for idx=1:length(V_h)
%         V=V_h{idx}+eps;
%         H=H_1_h{idx};
%         WH=WH_1_h{idx} + WH_2_h{idx} + WH_3_h{idx} + WH_4_h{idx}+ WH_5_h{idx}+ WH_6_h{idx}+eps;
%         cost_overall=(options_nmf.W1/length(V_h))*nmf_cost(options_nmf.beta_loss,V,WH,H,options_nmf.sparsity)+cost_overall;
%     end
%     
%     for idx=1:length(V_n)
%         V=V_n{idx}+eps;
%         H=H_3_n{idx};
%         WH=WH_1_n{idx} + WH_2_n{idx} + WH_3_n{idx} + WH_4_n{idx}+ WH_5_n{idx}+ WH_6_n{idx}+eps;
%         cost_overall=(options_nmf.W3/length(V_n))*nmf_cost(options_nmf.beta_loss,V,WH,H,options_nmf.sparsity)+cost_overall;
%     end
%     
%     for idx=1:length(V_s)
%         V=V_s{idx}+eps;
%         H=H_5_s{idx};
%         WH=WH_1_s{idx} + WH_2_s{idx} + WH_3_s{idx} + WH_4_s{idx}+ WH_5_s{idx}+ WH_6_s{idx}+eps;
%         cost_overall=(options_nmf.W5/length(V_s))*nmf_cost(options_nmf.beta_loss,V,WH,H,options_nmf.sparsity)+cost_overall;
%     end
%     
%     for idx=1:length(V_r)
%         V=V_r{idx}+eps;
%         H=H_6_r{idx};
%         WH=WH_1_r{idx} + WH_2_r{idx} + WH_3_r{idx} + WH_4_r{idx}+ WH_5_r{idx}+ WH_6_r{idx}+eps;
%         cost_overall=(options_nmf.W3/length(V_n))*nmf_cost(options_nmf.beta_loss,V,WH,H,options_nmf.sparsity)+cost_overall;
%     end
%     
%     for idx=1:length(V_m)
%         V=V_m{idx}+eps;
%         H=[H_1_m{idx};H_2_m{idx};H_3_m{idx};H_4_m{idx};H_5_m{idx};H_6_m{idx}];
%         W=[W_1,W_2,W_3,W_4,W_5,W_6];
%         WH=W*H+eps;
%         cost_overall=nmf_cost(options_nmf.beta_loss,V,WH,H,options_nmf.sparsity)+cost_overall;
%     end
%     
%     cost(i)=cost_overall;
end
end











