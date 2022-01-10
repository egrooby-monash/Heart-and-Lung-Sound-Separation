% mixed= 10 second noisy neonatal chest sound recording

%% Existing Methods
% Adaptive Line Enhancement
mu=0.000001;
L=256;
Delta = 10;
iterations=50;
adapt_method='Normalized LMS';
[heart,lung]=adaptive_line_enhancement(mixed,mu,L,Delta,iterations,adapt_method);

% RLS filter
[heart,lung]=rls_filter(mixed,4000);

% WTST FIlter
[heart,lung]=WTST_NST_filter(mixed,4000,3.5);

% Frequency interpolation
[heart,lung]=frequency_interpolation(mixed,4000,'springer','all','linear');

% Modulation filtering
[heart,lung]=modulation_filtering(mixed,4000,[4 20]);

% Singular Spectrum Analysis
[heart,lung]=singular_spectrum_analysis(mixed,4000);

% Standard Filtering 
[heart,lung]=recursive_filter(mixed,4000) ;

% Adaptive Fourier Transform 
[heart,lung]=adf_filtering(mixed,4000,'springer',5);

% EMD 
[heart,lung]=emd_separation(mixed,4000,'eemd','log energy');

% Wavelet SSA
heart=wavelet_ssa(mixed,4000);

% NMF Cluster 1
xhat=nmf_cluster1(mixed,4000);
heart=xhat(:,1);
lung=xhat(:,2);

% NMF Cluster 2
% heart_path= folder containing reference clean heart sounds
[Wh]=nmf_cluster2_training(heart_path);
xhat= nmf_cluster2(mixed,4000,Wh);
heart=resample(xhat(:,1),4000,8000);
lung=resample(xhat(:,2),4000,8000);


%% NMCF
respiratory_support="Bubble";
options_tf.FFTSIZE = 1024;
options_tf.HOPSIZE = 256;
options_tf.WINDOWSIZE =512;
max_examples=10;
% heart_path= folder with clean heart sound examples
% lung_path= folder with clean lung sound examples
% cry_path= folder with cry noises 
% stmv_path= folder with stethoscope movement noises
% bubble_path= folder with bubble cpap respiratory support noises
% cpap_path= folder with ventilator cpap respiratory support noises

TF='STFT';

options_nmf.W1=0;
options_nmf.W3=1;
options_nmf.W4=0;
options_nmf.W5=0.25;
options_nmf.W6=0.25;

options_nmf.beta_loss=1;
options_nmf.sparsity=0.1;
MAXITER = 100;
K=[20 10 20 20 20];

[V_h, V_l, V_n,V_s, V_r, W_1, W_3, W_4, W_5, W_6]=...
    load_example2(respiratory_support, TF,options_tf, max_examples,heart_path,lung_path,cry_path,stmv_path,bubble_path,cpap_path,options_nmf,MAXITER,K);

nmf_method="nmcf";
supervision="HS LS NS NU";
reconstruction="Filtering";

xhat=...
    nmcf_overall2(mixed,4000,nmf_method, supervision, TF,options_tf, reconstruction,K, options_nmf, MAXITER, W_1, [], W_3, W_4, W_5, W_6, V_h, V_l, V_n, V_s, V_r);
heart=xhat(:,1);
lung=xhat(:,4);
heart=heart/max(abs(heart));
lung=lung/max(abs(lung));

%% NMF 
respiratory_support="Bubble";
options_tf.FFTSIZE = 1024;
options_tf.HOPSIZE = 256;
options_tf.WINDOWSIZE =512;
max_examples=10;
% heart_path= folder with clean heart sound examples
% lung_path= folder with clean lung sound examples
% cry_path= folder with cry noises 
% stmv_path= folder with stethoscope movement noises
% bubble_path= folder with bubble cpap respiratory support noises
% cpap_path= folder with ventilator cpap respiratory support noises

TF='STFT';

options_nmf.W1=0;
options_nmf.W3=1;
options_nmf.W4=0;
options_nmf.W5=0.25;
options_nmf.W6=0.25;

options_nmf.beta_loss=1;
options_nmf.sparsity=0.1;
MAXITER = 100;
K=[20 10 20 20 20];

[V_h, V_l, V_n,V_s, V_r, W_1, W_3, W_4, W_5, W_6]=...
    load_example2(respiratory_support, TF,options_tf, max_examples,heart_path,lung_path,cry_path,stmv_path,bubble_path,cpap_path,options_nmf,MAXITER,K);

nmf_method='nmf';
supervision="HS LS NS NU";
reconstruction="Filtering";

xhat=...
    nmcf_overall2(mixed,4000,nmf_method, supervision, TF,options_tf, reconstruction,K, options_nmf, MAXITER, W_1, [], W_3, W_4, W_5, W_6, V_h, V_l, V_n, V_s, V_r);
heart=xhat(:,1);
lung=xhat(:,4);
heart=heart/max(abs(heart));
lung=lung/max(abs(lung)); 

 