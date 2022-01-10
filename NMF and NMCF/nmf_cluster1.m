function xhat= nmf_cluster1(xm,fs,K,reconstruction,TF,nmf_method,MAXITER,options_tf,options_nmf)
%% Paper Information
% On the blind recovery of cardiac and respiratory sounds
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6879427  
%% Purpose
% Sound separation 
%% Inputs
% xm= mixture signel to be separated
% fs= sampling frequency
% k= number of components to separate it into in the format [a b] which
% means separated into 2 components with 1st components having a bases and
% second component having b bases
% Reconstruction= method of reconstructing separated signals in time domain
% MAXITER= maximum iterations 
% nmf_method= to have sparse nmf or not
% options_nmf= parameter settings related to NMF e.g. .sparsity (activation L1 penalty) and
% beta_loss (beta loss for cost function)
% TF= time-frequency represetnation
% options_tf= parameter settings related to time-frequency representation- FFTSIZE, HOPSIZE, WINDOWSIZE= parameters of stft
%% Output
% xhat which contains heart and lung in first and second columns
% respectively

if nargin<3
    K=20;
    reconstruction="Best Mask";
    TF="STFT";
    xm=resample(xm,4000,fs);
    fs=4000;
    switch TF
        case "STFT"
            options_tf.FFTSIZE = 1024;
            options_tf.WINDOWSIZE = round(23.2/(1/fs)/1000);
            options_tf.HOPSIZE = round(options_tf.WINDOWSIZE*0.5);
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
            options_nmf.beta_loss=2;
        case "nmf_sparse"
            options_nmf.beta_loss=2;
            options_nmf.sparsity=0;
    end
    MAXITER = 130;
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
        [W, ~]   = nmf_supervised(Vm, K, [], MAXITER, []);
    case "nmf_general"
        [W, ~]=nmf_supervised_beta(Vm,K,[],MAXITER,[],options_nmf.beta_loss);
    case "nmf_sparse"
        [W,~]=nmf_supervised_sparse(Vm,K,[],MAXITER,[],options_nmf.beta_loss,options_nmf.sparsity);
end


%% Select initial reference basis spectrum option 2
Wt=W;
iter=K/2-1;

f=(0:1/(options_tf.FFTSIZE/2):1)*(fs/2);

cl=sum(W(f>=250&f<=1000,:),1);
ch=sum(W(f<=250,:),1);
%cl=sum(W(f>=300,:),1);
%ch=sum(W(f<=100,:),1);

[maxl,locl]=max(cl);
[maxh,loch]=max(ch);

if locl==loch
    if maxh>=maxl
        cl(locl)=-100;
        [~,locl]=max(cl);
    else
        ch(loch)=-100;
        [~,loch]=max(ch);
    end
end

Wl=Wt(:,locl);
Wh=Wt(:,loch);
Wt(:,[locl,loch])=[];

%% Clustering
for i=1:iter
    % 2. Assign closest basis vector to each reference basis spectrum
    % - define as correlation values
    % - spectral correlation as Wk vs Wi ref
    % - correlation between basis vectors of the mixed and reference
 
    cl=Wt'*Wl./(sqrt(sum(Wt.^2,1))'*sqrt(sum(Wl.^2,1)));
    ch=Wt'*Wh./(sqrt(sum(Wt.^2,1))'*sqrt(sum(Wh.^2,1)));
    
    [maxl,locl]=max(cl);
    [maxh,loch]=max(ch);
    if locl==loch
        if maxh>=maxl
            cl(locl)=-100;
            [~,locl]=max(cl);
        else
            ch(loch)=-100;
            [~,loch]=max(ch);
        end
    end
    % 3. Compute new RBS by adding the newly asigned basis to the previous RBS
    
    Wl=i*Wl/(i+1)+Wt(:,locl)/(i+1);
    Wh=i*Wh/(i+1)+Wt(:,loch)/(i+1);
    Wt(:,[locl,loch])=[];
end


%% Return final reference basis spectrum as complete BS of sources and perform NMF again
Wnew=[Wh,Wl];
K=[1,1];
fixedInds = 1:sum(K);
switch nmf_method
    case "nmf_kl"
        [W, H]   = nmf_supervised(Vm, K, Wnew, MAXITER, fixedInds);
    case "nmf_general"
        [W, H]=nmf_supervised_beta(Vm,K,Wnew,MAXITER,fixedInds,options_nmf.beta_loss);
    case "nmf_sparse"
        [W,H]=nmf_supervised_sparse(Vm,K,Wnew,MAXITER,fixedInds,options_nmf.beta_loss,options_nmf.sparsity);
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



