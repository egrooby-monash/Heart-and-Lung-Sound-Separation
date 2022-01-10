function xhat=...
    nmcf_overall2(xm,fs,nmf_method, supervision, TF,options_tf, reconstruction,K, options_nmf, MAXITER, W_1, W_2, W_3, W_4, W_5, W_6, V_h, V_l, V_n, V_s, V_r)
%% Purpose
% Perform NMCF with the following setup
% Heart Supervised
% Lung Supervsied
% Cry Noise, Stethoscope Movement Noise and Respiratory Support Noise Supervised
% Noise Unsupervised
%% Inputs
% xm= mixture recording
% fs= sampling frequency
% nmf_method= nmf or nmcf method
% TF= time-frequency representation
% options_tf= parameter settings related to time-frequency representation
% reconstruction= method of reconstructing separated signals in time domain
% K is inner dimension set with length 6
% options_nmf= parameter settings related to NMF e.g. .sparsity (activation L1 penalty) and
% beta_loss (beta loss for cost function)
% max_iter is maximum number of iterations
% W_1, W_2, W_3, W_4, W_5, W_6= base matrices for heart, noise
% unsupervised, noise (cry),lung, stethoscope movement and respiratory support respectively
% V_h is heart segments 
% V_l is lung segments 
% V_n is noise segments 
% V_s is stethoscope movement segments
% V_r is respiratory support segments
%% Outputs
% xhat= column matrix of separated sounds which contains heart, noise unsupervised, noise
% and lung respectively                   
if nargin<1
end

% Mixture Examples
g=[];
fshifts=[];
switch TF
    case "STFT"
        Xm = myspectrogram_modified(xm,options_tf.FFTSIZE,fs,hann(options_tf.WINDOWSIZE),-options_tf.HOPSIZE);
        Vm = abs(Xm(1:(options_tf.FFTSIZE/2+1),:));
        F = size(Vm,1);
        T = size(Vm,2);
    case "Cochleagram"
        %% Paper Information
        % Unsupervised SIngle-Channel Separation of Non-Stationary Signals using
        % Gammatone Filterbank and Itakura Saito Non-negative Matrix
        % Two-Dimensional Factorization
        % https://www.researchgate.net/profile/Wai-Lok-Woo/publication/236007062_Unsupervised_Single-Channel_Separation_of_Nonstationary_Signals_Using_Gammatone_Filterbank_and_Itakura-Saito_Nonnegative_Matrix_Two-Dimensional_Factorizations/links/5efe1aed45851550508533a6/Unsupervised-Single-Channel-Separation-of-Nonstationary-Signals-Using-Gammatone-Filterbank-and-Itakura-Saito-Nonnegative-Matrix-Two-Dimensional-Factorizations.pdf
        % https://www.mathworks.com/matlabcentral/fileexchange/48622-cochleagram-and-is-nmf2d-for-blind-source-separation
        
        % gammatone filterbank (cochleagram) which is more separable than the
        % spectrogram or log-frequency spectrogram (constant Q-transform)
        % Construct the cochleagram use Gammatone filterbank
        gf = gammatoneFast_modified(xm,options_tf.num_bin,options_tf.fRange,fs,options_tf.filterOrder,options_tf.gL);
        Xm = cochleagram_modified(gf,options_tf.cgL);
        Vm=abs(Xm);
        F = size(Vm,1);
        T = size(Vm,2);
    case "Q-transform"
        [Xm,~,g,fshifts] = cqt(xm,'SamplingFrequency' ,fs,'BinsPerOctave',options_tf.bins);
        Vm = abs(Xm);
        F = size(Vm,1);
        T = size(Vm,2);
    case "Q-transform2"
        freq=py.librosa.cqt(py.numpy.array(xm),pyargs(...
            'sr',int32(fs),...
            'n_bins',int32(options_tf.total_bin),...
            'bins_per_octave',int32(options_tf.bins),...
            'sparsity',options_tf.sparsity,...
            'hop_length',int32(options_tf.HOPSIZE)));
        Xm=complex(double(py.numpy.array(freq.real)),double(py.numpy.array(freq.imag)));
        Vm = abs(Xm);
        F = size(Vm,1);
        T = size(Vm,2);
        
    case "STFT2"
        xm(40001:40192)=0;
        [Xm,f] = stft(xm,fs,'Window',hann(options_tf.WINDOWSIZE,'periodic'),...
            'OverlapLength',options_tf.WINDOWSIZE-options_tf.HOPSIZE,'FFTLength',...
            options_tf.FFTSIZE,'FrequencyRange','onesided');
        Vm = abs(Xm);
        F = size(Vm,1);
        T = size(Vm,2);
end
V_m{1}=Vm;
V_m_pha{1}=angle(Xm);
%[V_m,V_m_pha]=read_files(mixture_path,TF,options_tf,[],[],max_examples);

% W_1 heart, W_2 unsupervised, W_3 noise, W_4 lung
switch supervision
    case "HS LU"
        % W_1 heart, W_2 unsupervised, W_3 NA, W_4 NA
        K(3) = 0;
        K(4) = 0;
        W_3=[];
        W_4=[];
    case "HS LS NU"
        % W_1 heart, W_2 unsupervised, W_3 NA, W_4 lung
        K(3) = 0;
        W_3=[];
        W_4=[];
    case "HS LU NS"
        % W_1 heart, W_2 unsupervised, W_3 noise, W_4 NA
        K(4) = 0; 
        W_4=[];
    case "HS LS NS"
        % W_1 heart, W_2 NA, W_3 noise, W_4 lung
        K(2) = 0; 
    case "HS LS NS NU"
        % W_1 heart, W_2 unsupervised, W_3 noise, W_4 lung
end

switch nmf_method
    case "nmcf"
        [W_1, W_2, W_3, W_4, W_5, W_6, H_1_m, H_2_m, H_3_m, H_4_m, H_5_m, H_6_m] = ...
    nmcf_general2(V_h, V_l, V_n, V_s, V_r, V_m, K, options_nmf, MAXITER);
    case "nmf"
        fixedInds=1:(K(1)+K(3)+K(4)+K(5)+K(6)); 
        W_2=rand(size(W_1,1), K(2)); 
        [W,H_all,~]=nmf_supervised_sparse_multi(V_m,K,[W_1 W_3 W_4 W_5 W_6 W_2],MAXITER,fixedInds,options_nmf.beta_loss,options_nmf.sparsity);
        c = [0 cumsum([K(1) K(3) K(4) K(5) K(6) K(2)])];
        W_2=W(:,c(6)+1:c(7)); 
        H_1_m{1}=H_all{1}(c(1)+1:c(2),:); 
        H_3_m{1}=H_all{1}(c(2)+1:c(3),:); 
        H_4_m{1}=H_all{1}(c(3)+1:c(4),:); 
        H_5_m{1}=H_all{1}(c(4)+1:c(5),:); 
        H_6_m{1}=H_all{1}(c(5)+1:c(6),:); 
        H_2_m{1}=H_all{1}(c(6)+1:c(7),:); 
        
end

xhat=compute_output2(reconstruction,TF,options_tf, V_m, V_m_pha, W_1, W_2, W_3, W_4, W_5, W_6, H_1_m, H_2_m, H_3_m, H_4_m, H_5_m, H_6_m,xm,g,fshifts,fs);

end
                 
    
