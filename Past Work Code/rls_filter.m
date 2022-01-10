function [heartSound,lungSound]=rls_filter(x,fs,filter_length,filter_band,heart_reference,method,weight_L,forgetting_factor,iterations)
%% Paper Information
%Recursive Least Squares Adaptive Noise Cancellation Filtering for Heart Sound Reduction in Lung Sounds Recordings
%https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1280403 
%% Purpose
% To remove heart sound interference from lung sound recordings
%% Inputs
% x= mixutre
% fs= sampling frequency
% filter_length= length of FIR filter to obtain reference heart sound
% filter_band= frequency band to obtain S1 and S2 heart peaks
% heart_reference= method to obtain heart sound reference
% method= method of adaptive filter e.g. RLS, LMS etc. 
% weight_L= adaptive filter length
% forgetting_factor= forgetting factor of RLS type adaptive filters
% iterations= iterations of the adaptive filter to stablise
%% Outputs
% Separated heartSound and lungSound

if nargin<3
    filter_length=400; 
    filter_band=[20 300];
    heart_reference='springer';
    
    method='Conventional RLS'; 
    weight_L = 2;
    forgetting_factor=1;
    iterations=50; 
end

% filter to get heart sound
b = fir1(filter_length,filter_band/(fs/2));
bandpassed= filtfilt(b,1,x);
% find S1 and S2 heart peaks
heart_location=heart_localisation(bandpassed,fs,heart_reference); 
bandpassed(heart_location==0)=0;

% run adaptive filter to obtain heart and lung sounds
w = zeros(weight_L,1);
for K=1:iterations
    RLS = dsp.RLSFilter('Length',weight_L, 'Method',method,...
    'InitialCoefficients',w,'ForgettingFactor',forgetting_factor);
    [heartSound, lungSound] = RLS(bandpassed, x);
    w=RLS.Coefficients;
end

end
