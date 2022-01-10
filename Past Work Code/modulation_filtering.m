function [heartSound, lungSound]=modulation_filtering(x,fs,filter_band)        
%% Paper Information
% Modulation filtering for heart and lung sound separation from breath sound recordings
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4649547 
% Study of lung and heart diagnosis application based on frequency separation of breath sound
% https://www.sciencedirect.com/science/article/pii/S221201731300457X
%% Purpose
% To separate heart and lung sounds
%% Inputs
% x= mixture
% fs= sampling frequency
% filter_band= modulation frequency pass and stop band for filters
%% Outputs
% Separate heart and lung sounds

if nargin<3
    filter_band=[1 20];
end

% STFT
[s,~] =stft(x,fs,'Window',hann(fs*0.02,'periodic'),'OverlapLength',fs*0.01);

% Modulation bandpass and bandstop filters
fs_freq=100;
filter_length=150;
b_bandpass = fir1(filter_length,filter_band/(fs_freq/2));
b_bandstop = fir1(filter_length,filter_band/(fs_freq/2),'stop');

% Filtering magnitude of STFT
XmagHat  = abs(s);
XmagHat_heart = zeros(size(s)); 
XmagHat_lung = zeros(size(s)); 
for i=1:size(s,1)
    XmagHat_heart(i,:)= filtfilt(b_bandpass,1,XmagHat(i,:));
    XmagHat_lung(i,:)= filtfilt(b_bandstop,1,XmagHat(i,:));
end

% Half-wave recitifier
XmagHat_heart(XmagHat_heart<0)=0; 
XmagHat_lung(XmagHat_lung<0)=0; 

% Reconstructing STFT for heart and lung
phi = angle(s);
XHat_heart = XmagHat_heart.*exp(1i*phi);
XHat_lung = XmagHat_lung.*exp(1i*phi);

heartSound=real(istft(XHat_heart,'Window',hann(fs*0.02,'periodic'),'OverlapLength',fs*0.01));
lungSound=real(istft(XHat_lung,'Window',hann(fs*0.02,'periodic'),'OverlapLength',fs*0.01));

% High-pass filter to remove some artefacts
filter_length=150;
b = fir1(filter_length,50/(fs/2),'high');
heartSound=filtfilt(b,1,heartSound);
lungSound=filtfilt(b,1,lungSound);
end





