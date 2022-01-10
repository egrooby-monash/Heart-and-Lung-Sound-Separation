function xhat= nmf_cluster2(xm,fs,Wh,K,Kc,Kl,reconstruction,TF,nmf_method,MAXITER,options_tf,options_nmf,heart_reference)
%% Paper Information
% A non-negative matrix factorization approach based on spectro-temporal clustering to extract heart sounds
% https://www-sciencedirect-com.ezproxy.lib.monash.edu.au/science/article/pii/S0003682X16304923 
%% Purpose
% Sound separation
%% Inputs
% xm= mixture signel to be separated
% fs= sampling frequency
% K= total number of bases to decompose mixture recording into
% Kc and Kl are number of bases that represent heart and lung sounds
% respectively
% reconstruction= method of reconstructing separated signals in time domain
% MAXITER= maximum iterations 
% nmf_method= to have sparse nmf or not
% options_nmf= parameter settings related to NMF e.g. .sparsity (activation L1 penalty) and
% beta_loss (beta loss for cost function)
% TF= time-frequency represetnation
% options_tf= parameter settings related to time-frequency representation- FFTSIZE, HOPSIZE, WINDOWSIZE= parameters of stft
% heart_reference= reference heart sound bases
%% Output
% xhat which contains heart and lung sounds in first and second columns
% respectively 
if nargin<4
    xm=resample(xm,8000,fs);
    fs=8000; 
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
    nmf_method="nmf_sparse";
    switch nmf_method
        case "nmf_kl"
            
        case "nmf_general"
            options_nmf.beta_loss=1;
        case "nmf_sparse"
            options_nmf.beta_loss=1;
            options_nmf.sparsity=0.15;
    end
    MAXITER = 500;
    Kc=55;
    Kl=64;
    K=Kc+Kl;
    heart_reference= 'springer';
end

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

%% NMF
switch nmf_method
    case "nmf_kl"
        [W, H]   = nmf_supervised(Vm, K, [], MAXITER, []);
    case "nmf_general"
        [W, H]=nmf_supervised_beta(Vm,K,[],MAXITER,[],options_nmf.beta_loss);
    case "nmf_sparse"
        [W,H]=nmf_supervised_sparse(Vm,K,[],MAXITER,[],options_nmf.beta_loss,options_nmf.sparsity);
end

% Spectral Correlation
SC=W'*Wh./(sqrt(sum(W.^2,1))'*sqrt(sum(Wh.^2,1)));
SC=max(SC,[],2)';
Uf=0.4;
SCh=SC>=Uf;

% Roll off
phi = angle(Xm);
for i=1:K
    XmagHat = W(:,i)*H(i,:);
    
    % create upper half of frequency before istft
    XmagHat = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
    
    % Multiply with phase
    XHat = XmagHat.*exp(1i*phi);
    xhat1(:,i) = real(invmyspectrogram(XHat,options_tf.HOPSIZE))';
end

ROh=zeros(1,size(xhat1,2));
RO_frac=zeros(1,size(xhat1,2));
for i=1:size(xhat1,2)
    [s,f]=stft(xhat1(:,i),fs);
    power= sum(abs(s).^2,2);
    f0=260;
    RO=sum(power(f>=0 & f<=f0));
    TP=sum(power(f>=0 & f<=fs/2));
    RO_frac(i)=RO/TP;
    if RO>=0.85*TP
        ROh(i)=1;
    end
end

% Temporal correlation
filter_length=400;
filter_band=[50 250];
b = fir1(filter_length,filter_band/(fs/2));
bandpassed= filtfilt(b,1,xm);
heartSound=heart_localisation(bandpassed,fs,heart_reference); 

bandpassed(heartSound==0)=0; 
s2 =myspectrogram_modified(bandpassed,options_tf.FFTSIZE,fs,hann(options_tf.WINDOWSIZE),-options_tf.HOPSIZE);
power=sum(abs(s2));
P=power>mean(power);


H_new=H>mean(H,2);

TC= (1/(T-1))*(H_new-mean(H_new,2)./std(H_new,[],2))*(P-mean(P)/std(P))';

Ut=0;
TCh=TC>=Ut;
TCh=TCh';


temp=normalize(TC')+normalize(SC)+normalize(RO_frac); 
[~,I]=sort(temp); 

W1=W(:,I(65:end));
W2=W(:,I(1:64));
H1=H(I(65:end),:);
H2=H(I(1:64),:);
W=[W1,W2];
H=[H1;H2];
K=[Kc Kl]; 

% Clustering
%W1=W(:,TCh|ROh|SCh);
%W2=W(:,~(TCh|ROh|SCh));
%H1=H(TCh|ROh|SCh,:);
%H2=H(~(TCh|ROh|SCh),:);
%W=[W1,W2];
%H=[H1;H2];

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




