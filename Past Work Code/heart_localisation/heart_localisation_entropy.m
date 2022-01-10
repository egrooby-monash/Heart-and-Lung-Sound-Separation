function [e,heart_location]=heart_localisation_entropy(x,fs)
%% Paper Information
% A robust method for heart sound localization using lung entropy
% https://ieeexplore.ieee.org/abstract/document/1597500/ 
%% Purpose
% Determine location of S1 and S2 heart sounds
%% Inputs
% x= mixture recording
% fs= sampling frequency
%% Outputs
% e= Shannon entropy
% heart_location= 1s where S1 and S2 sounds occur and zeros elsewhere

% parameters for sliding window
window_seconds = 20*10^-3;                
window_length = window_seconds*fs;
hop_length = 0.5*window_length; 
num_samples= floor((length(x)-window_length)/hop_length+1);
e=zeros(num_samples,1); 
for i=1:num_samples
    % Going through sliding window with original signal
    y_seg = x((i-1)*hop_length+1:(i-1)*hop_length+window_length);
    e(i) = wentropy(y_seg,'shannon');
end

heart_location=zeros(size(e));
heart_location(e>mean(e)+std(e))=1; 
heart_location=round(resample(heart_location,fs,100));

end

