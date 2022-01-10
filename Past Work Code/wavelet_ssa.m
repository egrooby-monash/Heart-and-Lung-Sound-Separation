function heartSound=wavelet_ssa(x,fs1)
%% Paper Information
% A Noise Reduction Technique Based on Nonlinear Kernel Function for Heart Sound Analysis
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7851051 
%% Purpose
% Obtain clean heart sounds
%% Inputs
% x= mixture
% fs= sampling frequency
%% Output
% Clean heart sound

% resample to 4000Hz
x = resample(x,4000,fs1); 
fs=4000; 

% Wavelet decomposition and reconstruction of a2 coefficients to get
% 0-500Hz frequency band

% no. of levels for SWT
N=6; 
% wavelet name= Symlet order 5
wname='coif5';
% Discrete stationary wavelet transform 1-D
% correcting the length of signal
fin=size(x,1)-mod(size(x,1),2^N);
x_shortened=x(1:fin);
[swa,swd] =swt(x_shortened,N,wname);
x2 = iswt(swa(3,:),swd(3,:),wname); 

% Singular spectrum analysis 
[heartSound,~]=singular_spectrum_analysis(x2,fs); 

% Adaptive threshold to obtain just S1 and S2 peaks 
hilbert_envelope = Hilbert_Envelope(heartSound, fs);
heart_locations=hilbert_envelope>std(hilbert_envelope); 
heartSound(~heart_locations)=0; 
end
