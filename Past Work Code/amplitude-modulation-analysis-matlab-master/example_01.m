%% Example 01
% This example shows the use of the GUI to explore Amplitude Modulation
% for ECG data and EEG data
%
% The 'explore_strfft_ama_gui()' computes the Modulation Spectrogram.  
%   It uses the Short Time Real Fourier Fast Transform (STRFFT) to compute 
%   the Spectrogram, after rFFT is used to obtain the Modulation Spectrogram
%
% The 'explore_wavelet_ama_gui()' computes the Modulation Spectrogram using 
%   It uses the Wavelet transform with Complex Morlet wavelet to compute 
%   the Spectrogram, after rFFT is used to obtain the Modulation Spectrogram
%
% Usage for explore_*_ama_gui()
%
% Once the GUI is executed, it accepts the following commands
%
%   Key          Action
%   Up Arrow     Previous channel                  (-1 channel) 
%   Down Arrow   Next channel                      (+1 channel)
%   Left Arrow   Go back to the previous segment   (-1 segment shift)  
%   Right Arrow  Advance to the next segment       (+1 segment shift)
%   'W'          Previous channel                  (-5 channels) 
%   'S'          Next channel                      (+5 channel)
%   'A'          Go back to the previous segment   (-5 segment shift) 
%   'D'          Advance to the next segment       (+5 segment shift)
% 
%   'U'          Menu to update:
%                   parameters for Modulation Spectrogram 
%                   ranges for conventional and modulation frequency axes
%                   ranges for power in Spectrogram and Modulation Spectrogram 
%   ESC          Close the GUI
%

%% ECG data (1 channel) using STFFT-based Modulation Spectrogram
load('./example_data/ecg_data.mat');

% STFFT Modulation Spectrogram
explore_stfft_ama_gui(x, fs, 'ECG', 'jet');

%% ECG data (1 channel) using wavelet-based Modulation Spectrogram
load('./example_data/ecg_data.mat');

% Wavelet Modulation Spectrogram
explore_wavelet_ama_gui(x, fs, 'ECG');

%% EEG data (7 channels) using STFFT-based Modulation Spectrogram
load('./example_data/eeg_data.mat');

% STFFT Modulation Spectrogram
explore_stfft_ama_gui(x, fs, ch_names);

%% EEG data (7 channels) using wavelet-based Modulation Spectrogram
load('./example_data/eeg_data.mat');

% Wavelet Modulation Spectrogram
explore_wavelet_ama_gui(x, fs, ch_names);
