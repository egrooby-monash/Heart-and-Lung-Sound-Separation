function [V_h, V_l, V_n,V_s, V_r, W_1, W_3, W_4, W_5, W_6]=load_example2(respiratory_support, TF,options_tf, max_examples,heart_path,lung_path,cry_path,stmv_path,bubble_path,cpap_path,options_nmf,MAXITER,K)
%% Purpose
% Pretrain nmf/nmcf method with reference sounds
%% Inputs
% respiratory_support= type of respiratory support mixture recording was
% TF= time frequency represetnation
% options_tf= parameters for time frequency representation
% max_examples= max_examples to read for each folder
% heart/lung/cry/stmv/bubble/cpap_path= folder location for these reference
% files
% options_nmf= parameters for nmf/nmcf method
% MAXITER= maximum number of iterations
% K= number of bases
if nargin<1
    respiratory_support="none"; 
    TF="STFT"; 
    max_examples=20; 
end 

% Heart Examples
[V_h,~]=read_files(heart_path,TF,options_tf,[],[],max_examples);
    
% Lung Examples  
[V_l,~]=read_files(lung_path,TF,options_tf,[],[],max_examples);

% Noise Examples 
% cry
[V_n, ~] = read_files(cry_path,TF,options_tf,[],[],max_examples);
% stethoscope
[V_s, ~] = read_files(stmv_path,TF,options_tf, [], [],max_examples);
% V_n=[];
% V_n_pha=[]; 
%respiratory support
switch respiratory_support
    case "none"
    case "CPAP"
        [V_r, ~] = read_files(cpap_path,TF,options_tf, [], [],max_examples);
    case "Bubble"
        [V_r, ~] = read_files(bubble_path,TF,options_tf, [], [],max_examples);
end

% nmf pretraining
% W_1 heart, W_2 unsupervised, W_3 noise, W_4 lung
[W_1,~,~]=nmf_supervised_sparse_multi(V_h,K(1),[],MAXITER,[],options_nmf.beta_loss,options_nmf.sparsity);
[W_3,~,~]=nmf_supervised_sparse_multi(V_n,K(3),[],MAXITER,[],options_nmf.beta_loss,options_nmf.sparsity);
[W_4,~,~]=nmf_supervised_sparse_multi(V_l,K(4),[],MAXITER,[],options_nmf.beta_loss,options_nmf.sparsity);
[W_5,~,~]=nmf_supervised_sparse_multi(V_s,K(5),[],MAXITER,[],options_nmf.beta_loss,options_nmf.sparsity);
[W_6,~,~]=nmf_supervised_sparse_multi(V_r,K(6),[],MAXITER,[],options_nmf.beta_loss,options_nmf.sparsity);
end
