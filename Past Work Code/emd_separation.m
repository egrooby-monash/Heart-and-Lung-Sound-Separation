function [heartSound,lungSound]=emd_separation(x,fs,emd_method,heart_reference)
%% Paper Information
% Reduction of heart sound interference from lung sound signals using empirical mode decomposition technique
% https://www.tandfonline.com/doi/abs/10.3109/03091902.2011.595529
%% Purpose
% Separation of heart and lung sounds
%% Inputs
% x= mixture recording
% fs= sampling frequency
% emd_method= method of empirical mode decomposition
%  - 'MATLAB emd'= matlab's verion of empircal mode decomposition
%  - 'emd'= empirical mode decomposition
%  - 'eemd'= ensemble empirical mode decomposition
%  - 'ceemd'= complete ensemble empirical mode decomposition
% heart_reference= method to obtaining S1 and S2 peak locations
%% Outputs
% Separation heartSounds and lungSounds

if nargin<3
    emd_method='ceemd'; 
    heart_reference='log energy';
end

% normalise
x=x/max(abs(x));

% decompose using empirical mode decomposition
% max iterations
MaxIter = 100;
% standard deviation of noise
Nstd = 0.2;
% number of realisations
NR = 5;
% if equals 1, then the SNR increases for every stage, as in [1].
% If equals 2, then the SNR is the same for all stages, as in [2].
SNRFlag=1; 
switch emd_method
    case 'MATLAB emd'
        [imf,~] = emd(x);
    case 'emd'
        [imf, ~, ~]=emd_alt(x,'MAXITERATIONS',MaxIter);
         imf=imf'; 
    case 'eemd'
        [imf, ~]=eemd(x,Nstd,NR,MaxIter); 
        imf=imf'; 
    case 'ceemd'
        [imf,~]=ceemdan_improved(x,Nstd,NR,MaxIter,SNRFlag); 
        imf=imf'; 
end
 
% identify heart included segments (HI) and heart free segments in each IMF
HI=zeros(size(imf));
HF=zeros(size(imf));

for i=1:size(imf,2)
    % 50-250Hz bandpass filter
    filter_length=400;
    filter_band=[50 250];
    b = fir1(filter_length,filter_band/(fs/2));
    bandpassed= filtfilt(b,1,imf(:,i));
    
    % S1 and S2 peak detection
    heart_locations=heart_localisation(bandpassed,fs,heart_reference);
    
    % extend site of S1 and S2 peaks by 50ms back and 100ms forward
    change_states=diff(heart_locations);
    start_heart=find(change_states==1);
    new_start_heart=max(start_heart-round(50/1000*fs),1);
    
    for j=1:length(start_heart)
        heart_locations(new_start_heart(j):start_heart(j))=1;
    end
    
    end_heart=find(change_states==-1);
    new_end_heart=min(end_heart+round(100/1000*fs),length(x));
    
    for j=1:length(end_heart)
        heart_locations(end_heart(j):new_end_heart(j))=1;
    end
    
    HI(:,i)= imf(:,i);
    HI(heart_locations==0,i)=0;
    HF(:,i)= imf(:,i);
    HF(heart_locations==1,i)=0;
end
H=sum(HI,2);
L=sum(HF,2);

% High pass H included segments to obtain lung information
[b,a] = butter(10,150/(fs/2),'high');
H_high = filtfilt(b,a,H);
[b,a] = butter(10,150/(fs/2));
heartSound = filtfilt(b,a,H);

lungSound=L+H_high;
end