function [signal, heart_location]=heart_localisation_wavelet(x)
%% Paper Information
% Heart Sound Cancellation Based on Multiscale Products and Linear
% Prediction 
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1404075
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4067109
%% Purpose
% Determine location of S1 and S2 heart sounds
%% Inputs
% x= mixture recording
%% Outputs
% signal= wavelet multiresolution product
% heart_location= 1s where S1 and S2 sounds occur and zeros elsewhere


% no. of levels for SWT
N=3;
% wavelet name= Symlet order 5
wname='sym5';
%Discrete stationary wavelet transform 1-D
% correcting the length of signal
fin=size(x,1)-mod(size(x,1),2^N);
x_shortened=x(1:fin);
[swa,~] =swt(x_shortened,N,wname);

signal=swa(1,:).*swa(2,:).*swa(3,:); 

mag_signal=abs(signal); 
heart_location=zeros(size(signal));
heart_location(mag_signal>mean(mag_signal)+std(mag_signal))=1; 
end