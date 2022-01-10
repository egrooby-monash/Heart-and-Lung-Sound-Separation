function [heartSound,lungSound]=WTST_NST_filter(X,fs1,Fadj,e,window_length,maxiter,type)
%% Paper Information
% A wavelet-based reduction of heart sound noise from lung sounds
% https://www-sciencedirect-com.ezproxy.lib.monash.edu.au/science/article/pii/S1386505698001373#BIB10 
% 
% Separation of Discontinuous Adventitious Sounds from Vesicular Sounds Using a Wavelet-Based Filter
% https://ieeexplore-ieee-org.ezproxy.lib.monash.edu.au/stamp/stamp.jsp?tp=&arnumber=649999 
% 
% AN OVERVIEW OF HEART-NOISE REDUCTION OF LUNG SOUND USING WAVELET TRANSFORM BASED FILTER
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1279719
%% Purpose
% Main purpose is to clean up the lung sound
%% Inputs
% X= mixture
% fs1= sampling frequency
% Fadj= hard thresholding factor 
% e= early stopping criterion
% window_length= length of window 
% maxiter= maximum iterations
% type= discrete or continuous WT filter 
%% Outputs
% Separated heartSound and lungSound

if nargin<3
    % Adjusting multiplicative factor for hard thresholding 
    % As stated in original paper Fadj=2.7
    % Though 3rd paper stated For low flow rate it was varied between 2.5 
    % and 2.7 and for medium flow rate it was adjusted between 2.7 to 4
    Fadj=2.7; 
end    
if nargin< 4
    % stopping criterion
    e=0.00001;
    window_length=2048; 
    maxiter=20; 
    %'continuous'; %'discrete'
    type='discrete';
end

% Resample to 2500Hz
X = resample(X,2500,fs1); 
fs=2500; 
% Normalise 
X=normalize(X);

N=length(X);
heartSound=zeros(size(X));
lungSound=zeros(size(X));
for i=1:ceil(N/window_length)
    disp('Second')
    disp(i)
    % window signal
    if i==ceil(N/window_length)
        f=X((i-1)*window_length+1:end);
    else
        f=X((i-1)*window_length+1:i*window_length);
    end 
    
    N2=length(f);
    % WT scales
    M=floor(log2(N2));
    %m=1:M;
    Ro=0; 
    %X=X(1:2^M);
    WT={};
    WTc={};
    WTr={};
    C={};
    R={};
        
    %Iteration part- usually 8 is enough
    for k=1:maxiter 
        % wavelet multi-resolution decomposition
        switch type
            case 'continuous'
                [WT{k},F]=cwt(f,fs,'VoicesPerOctave',4);
                elements=size(WT{k},1);
            case 'discrete'
                fin=N2-mod(N2,2^M);
                f=f(1:fin);
                % Decompose = WaveletTransform
                wname='db8';
                % discrete meyer wavelet frequency filters
                SWC=swt(f,M,wname);
                WT{k}=SWC;
                elements=M+1;
        end

        WTc{k}=WT{k};
        WTr{k}=WT{k};
        % loop over all wavelet levels
        for j=1:elements
            % hard thresholding
            sigma=std(WT{k}(j,:));
            THR=sigma*Fadj; 
            WTc{k}(j,abs(WT{k}(j,:))<THR)=0;
            WTr{k}(j,abs(WT{k}(j,:))>=THR)=0;
        end
        disp('iteration')
        
        % wavelet multi-resolution reconstruction
        switch type
            case 'continuous'
                C{k}=icwt(WTc{k},F,[min(F) max(F)])';
                R{k} = icwt(WTr{k},F,[min(F) max(F)])';
            case 'discrete'
                C{k} = iswt(WTc{k},'db8')';
                R{k} = iswt(WTr{k},'db8')';
        end
        
        disp(k)
        
        % check early stopping criteria
        if k==1
            STC=abs(mean(Ro^2)-mean(R{k}.^2));
        else
            STC=abs(mean(R{k-1}.^2)-mean(R{k}.^2));
        end
        if STC>=e
            f=R{k};
        else
            disp('complete')
            break
        end
    end
    
    %non-stationary part= heart sounds
    DAS=sum([C{:}],2); 
    %stationary part= lung sounds
    PVS=R{end}; 
     if i==ceil(N/window_length)
        heartSound((i-1)*window_length+1:(i-1)*window_length+length(DAS))=DAS;
        lungSound((i-1)*window_length+1:(i-1)*window_length+length(PVS))=PVS;
     else
        heartSound((i-1)*window_length+1:i*window_length)=DAS;
        lungSound((i-1)*window_length+1:i*window_length)=PVS;
    end 
end

% resample to fs1
heartSound = resample(heartSound,fs1,fs); 
lungSound = resample(lungSound,fs1,fs); 
heartSound=heartSound/max(abs(heartSound));
lungSound=lungSound/max(abs(lungSound));
end