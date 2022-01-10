function xhat= nmf_overall(xm,fs,K,reconstruction,TF,nmf_method,MAXITER,options_tf,options_nmf,supervised,x1,x2)
%% Paper Information
% Single-Channel Source Separation Tutorial Mini-Series
% https://ccrma.stanford.edu/~njb/teaching/sstutorial/
%% Purpose
% decompose mixture recording xm
%% Inputs
% xm= mixture signel to be separated
% fs= sampling frequency
% k= number of components to separate it into in the format [a b] which
% means separated into 2 components with 1st components having a bases and
% second component having b bases
% FFTSIZE, HOPSIZE, WINDOWSIZE= parameters of stft
% supervised is a 1x2 logic array which states if the situation is
% unsupervised, semi-supervised or fully supervised
% x1 and x2 are the reference signals
%% Output
% xhat which is signal length x k matrix. Each column represents a
% separated signal
if nargin<4
    reconstruction="Filtering";
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
            options_nmf.sparsity=0;
    end
    MAXITER = 100;
end
if nargin<9
    supervised=[false false];
end

%% Time Frequency Representation
switch TF
    case "STFT"
        Xm = myspectrogram_modified(xm,options_tf.FFTSIZE,fs,hann(options_tf.WINDOWSIZE),-options_tf.HOPSIZE);
        Vm = abs(Xm(1:(options_tf.FFTSIZE/2+1),:));
        if supervised(1)
            X1 = myspectrogram_modified(x1,options_tf.FFTSIZE,fs,hann(options_tf.WINDOWSIZE),-options_tf.HOPSIZE);
            V1 = abs(X1(1:(options_tf.FFTSIZE/2+1),:));
        end
        if supervised(2)
            X2 = myspectrogram_modified(x2,options_tf.FFTSIZE,fs,hann(options_tf.WINDOWSIZE),-options_tf.HOPSIZE);
            V2 = abs(X2(1:(options_tf.FFTSIZE/2+1),:));
        end
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
        if supervised(1)
            gf = gammatoneFast_modified(x1,options_tf.num_bin,options_tf.fRange,fs,options_tf.filterOrder,options_tf.gL);
            X1 = cochleagram_modified(gf,options_tf.cgL);
            V1=abs(X1);
        end
        if supervised(2)
            gf = gammatoneFast_modified(x2,options_tf.num_bin,options_tf.fRange,fs,options_tf.filterOrder,options_tf.gL);
            X2 = cochleagram_modified(gf,options_tf.cgL);
            V2=abs(X2);
        end
        F = size(Vm,1);
        T = size(Vm,2);
    case "Q-transform"
        [Xm,~,g,fshifts] = cqt(xm,'SamplingFrequency' ,fs,'BinsPerOctave',options_tf.bins);
        Vm = abs(Xm);
        if supervised(1)
            [X1,~,~,~] = cqt(x1,'SamplingFrequency' ,fs,'BinsPerOctave',options_tf.bins);
            V1 = abs(X1);
        end
        if supervised(2)
            [X2,~,~,~] = cqt(x2,'SamplingFrequency' ,fs,'BinsPerOctave',options_tf.bins);
            V2 = abs(X2);
        end
        F = size(Vm,1);
        T = size(Vm,2);
end

%% NMF
if supervised(1)
    switch nmf_method
        case "nmf_kl"
            [W1, ~] = nmf_supervised(V1, K(1), [], MAXITER,[]);
        case "nmf_general"
            [W1, ~]=nmf_supervised_beta(V1,K(1),[],MAXITER,[],options_nmf.beta_loss);
        case "nmf_sparse"
            [W1,~]=nmf_supervised_sparse(V1,K(1),[],MAXITER,[],options_nmf.beta_loss,options_nmf.sparsity);
    end
else
    W1 = 1+rand(F, K(1));
end

if supervised(2)
    switch nmf_method
        case "nmf_kl"
            [W2, ~] = nmf_supervised(V2, K(2), [], MAXITER,[]);
        case "nmf_general"
            [W2, ~]=nmf_supervised_beta(V2,K(2),[],MAXITER,[],options_nmf.beta_loss);
        case "nmf_sparse"
            [W2,~]=nmf_supervised_sparse(V2,K(2),[],MAXITER,[],options_nmf.beta_loss,options_nmf.sparsity);
    end
else
    W2 = 1+rand(F, K(2));
end


if supervised(1) && supervised(2)
    fixedInds = 1:sum(K);
elseif supervised(1)
    fixedInds = 1:K(1);
elseif supervised(2)
    fixedInds = (K(1)+1):sum(K);
else
    fixedInds = [];
end

switch nmf_method
    case "nmf_kl"
        [W, H]   = nmf_supervised(Vm, K, [W1 W2], MAXITER, fixedInds);
    case "nmf_general"
        [W, H]=nmf_supervised_beta(Vm,K,[W1 W2],MAXITER,fixedInds,options_nmf.beta_loss);
    case "nmf_sparse"
        [W,H]=nmf_supervised_sparse(Vm,K,[W1 W2],MAXITER,fixedInds,options_nmf.beta_loss,options_nmf.sparsity);
end

%% ISTFT / RECONSTRUCTION
switch reconstruction
    case "Synthesis"
        % get the mixture phase
        phi = angle(Xm);
        for i=1:K(1)
            switch TF
                case "STFT"
                    XmagHat = W(:,i)*H(i,:);
                    
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    
                    xhat(:,i) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                case "Q-transform"
                    XmagHat_full = W(:,i)*H(i,:);
                    
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    
                    xhat(:,i) = real(icqt(XHat,g,fshifts))';
            end
        end
    case "No Mask"
        % get the mixture phase
        phi = angle(Xm);
        c = [0 cumsum(K)];
        for i=1:length(K)
            switch TF
                case "STFT"
                    % filter
                    XmagHat = W(:,c(i)+1:c(i+1))*H(c(i)+1:c(i+1),:);
                    
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    
                    xhat(:,i) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                case "Q-transform"
                    % filter
                    XmagHat_full = W(:,c(i)+1:c(i+1))*H(c(i)+1:c(i+1),:);
                    
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    
                    xhat(:,i) = real(icqt(XHat,g,fshifts))';
            end
        end
    case "Filtering"
        % get the mixture phase
        phi = angle(Xm);
        c = [0 cumsum(K)];
        for i=1:length(K)
            
            % create masking filter
            Mask =  W(:,c(i)+1:c(i+1))*H(c(i)+1:c(i+1),:)./(W*H);
            
            switch TF
                case "STFT"
                    % filter
                    XmagHat = Vm.*Mask;
                    
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    
                    xhat(:,i) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                case "Q-transform"
                    % filter
                    XmagHat_full = Vm.*Mask;
                    
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    
                    xhat(:,i) = real(icqt(XHat,g,fshifts))';
                case "Cochleagram"
                    xhat(:,i) = synthesisFast(xm,Mask,options_tf.fRange,options_tf.cgL,fs);
            end
            
        end
    case "Best Mask"
        % get the mixture phase
        phi = angle(Xm);
        c = [0 cumsum(K)];
        for i=1:length(K)
            Rec(:,:,i)=W(:,c(i)+1:c(i+1))*H(c(i)+1:c(i+1),:);
        end
        for i=1:length(K)
            % create masking filter
            Mask = logical(Rec(:,:,i)==max(Rec,[],3));
            
            switch TF
                case "STFT"
                    % filter
                    XmagHat = Vm.*Mask;
                    
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    
                    xhat(:,i) = real(invmyspectrogram(XHat,options_tf.HOPSIZE))';
                case "Q-transform"
                    % filter
                    XmagHat_full = Vm.*Mask;
                    
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    
                    xhat(:,i) = real(icqt(XHat,g,fshifts))';
                case "Cochleagram"
                    xhat(:,i) = synthesisFast(xm,Mask,options_tf.fRange,options_tf.cgL,fs)';
            end
        end
end

end



