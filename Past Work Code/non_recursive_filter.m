function [heartSound,lungSound]=non_recursive_filter(x,fs,filter_length,heart_type,heart_cutoff,lung_type,lung_cutoff)
%% Paper Information
%Heart sound and lung sound separation algorithms: a review
% https://www.tandfonline.com/doi/pdf/10.1080/03091902.2016.1209589?needAccess=true

%high pass linear filter is developed for a cut-off frequency of 50-150Hz
%250th order FIR filter
%Butterworth
% Cutoff between 50-150Hz
%% Purpose
% Separate heart and lung with FIR filter
%https://www.mathworks.com/help/signal/ug/iir-filter-design.html
%% Inputs
% x= mixture
% fs= sampling frequency
% heart/lung_cutoff is frequency range of filter
% heart/lung_type is high/low/bandpass type of filter
%% Outputs
% Separated heartSound and lungSound

%% Dafault Values
if nargin<3
    filter_length=250;
    heart_type='bandpass';
    heart_cutoff=[50 250];
    
    lung_type='bandpass';
    lung_cutoff=[200 1000];
    %heart_type='low';
    %heart_cutoff=100;
    %
    %lung_type='high';
    %lung_cutoff=100;
end

b = fir1(filter_length,lung_cutoff/(fs/2),lung_type);
lungSound= filtfilt(b,1,x);
b = fir1(filter_length,heart_cutoff/(fs/2),heart_type);
heartSound= filtfilt(b,1,x);

heartSound=heartSound/max(abs(heartSound));
lungSound=lungSound/max(abs(lungSound));
end

