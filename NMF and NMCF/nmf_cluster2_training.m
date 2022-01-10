function [Wh]=nmf_cluster2_training(path,TF,options_tf,max_examples,nmf_method,options_nmf,MAXITER,Kg)
%% Paper Information
% A non-negative matrix factorization approach based on spectro-temporal clustering to extract heart sounds
% https://www-sciencedirect-com.ezproxy.lib.monash.edu.au/science/article/pii/S0003682X16304923 
%% Purpose
% Obtain reference basis vector from heart sound database
%% Inputs
% path= reference folder location
% TF= time frequency 
% options_tf= parameters for time frequency representation
% max_examples= maximum number of files to read in a folder
% nmf_method= to have sparse nmf or not
% options_nmf= parameter settings related to NMF e.g. .sparsity (activation L1 penalty) and
% beta_loss (beta loss for cost function)
% MAXITER= maximum iterations 
% Kg= number of bases the reference heart sounds are represented as 
%% Outputs
% Wh= reference heart frequency bases

subcount=1; 
if nargin<2
    max_examples=20; 
    TF="STFT";
    switch TF
        case "STFT"
            options_tf.FFTSIZE = 1024;
            options_tf.HOPSIZE = 256;
            options_tf.WINDOWSIZE = 512;
        case "Cochleagram"
            % fRange= frequency range
            options_tf.fRange = [50, fs/2];
            % num_bin= number of bins/channels for gammatone filterbank
            options_tf.num_bin = 128;
            % filterOrder= filter order of Gammatone filterbank
            options_tf.filterOrder = 4;
            % gammatone filter length or 128 ms for 16 kHz sampling rate
            options_tf.gL = round(0.128*fs);
            % default window length in sample points which is 20 ms for 16 KHz sampling frequency
            options_tf.cgL=0.02*fs;
        case 'Q-transform'
            % Number of bins per octave
            options_tf.bins=24;
    end
    nmf_method="nmf_general";
    switch nmf_method
        case "nmf_kl"
        case "nmf_general"
            options_nmf.beta_loss=1;
        case "nmf_sparse"
            options_nmf.beta_loss=1;
            options_nmf.sparsity=0.15;
    end
    MAXITER = 500;
    Kg=4;
end

files = dir(path);

% iterate over all reference heart sounds
Wh=[];
for i=1:length(files)
    % skip empty files
    if strcmp(files(i).name,'.')
        continue
    elseif strcmp(files(i).name,'..')
        continue
    elseif contains(files(i).name,'.png')
        continue
    elseif strcmp(files(i).name,'.DS_Store')
        continue
    end
    [xm,fs]=audioread(sprintf('%s/%s',path,files(i).name));
    xm=resample(xm,8000,fs); 
    fs=8000;
    %% Time Frequency Representation
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
    end
    switch nmf_method
        case "nmf_kl"
            [W, ~]   = nmf_supervised(Vm, Kg, [], MAXITER, []);
        case "nmf_general"
            [W, ~]=nmf_supervised_beta(Vm,Kg,[],MAXITER,[],options_nmf.beta_loss);
        case "nmf_sparse"
            [W,~]=nmf_supervised_sparse(Vm,Kg,[],MAXITER,[],options_nmf.beta_loss,options_nmf.sparsity);
    end
    if isempty(Wh)
        Wh=W;
    else
        Wh=[Wh W];
    end
    % break if max_examples is exceeded
    subcount=subcount+1;
    if subcount>max_examples
        return
    end
end
end
