function [e,heart_location]=heart_localisation_log_energy(x,fs)
%% Paper Information
% Lung-Heart Sound Separation Using Noise Assisted Multivariate Empirical Mode Decomposition
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6704645 
%% Purpose
% Determine location of S1 and S2 heart sounds
%% Inputs
% x= mixture recording
% fs= sampling frequency
%% Outputs
% e= Shannon entropy
% heart_location= 1s where S1 and S2 sounds occur and zeros elsewhere

filter_order=10; 
cutoff=100; 
[b,a] = butter(filter_order,cutoff/(fs/2));
x = filtfilt(b,a,x);

% parameters for sliding window
window_seconds = 12.5*10^-3;                
window_length = window_seconds*fs;
hop_length = 0.5*window_length; 
num_samples= floor((length(x)-window_length)/hop_length+1);
e=zeros(num_samples,1); 
for i=1:num_samples
    % Going through sliding window with original signal
    y_seg = x((i-1)*hop_length+1:(i-1)*hop_length+window_length);
    e(i) = log10((1/window_length)*sum(y_seg.^2)); 
end

heart_location=zeros(size(e));
heart_location(e>mean(e)+std(e))=1; 
heart_location=round(resample(heart_location,fs,160));

end
