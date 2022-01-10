%% Single channel speech enhancement
%% Parameters
data.denoise_type = 'spec_sub_mag'; % 'spec_sub_power', 'spec_sub_mag', 'spec_sub_over', 'wiener', 'mmse_stsa'
data.plot = 0;

%% Generate noisy signal
speech_file = 'wav/input/female2.wav';
noise_file = 'wav/input/noise_babble.wav';
[x, fs] = audioread(speech_file);
[n, fs_noise] = audioread(noise_file);
n = resample(n, fs, fs_noise);

% Adjust input SNR using noise_gain
noise_gain = 0.7;
n = noise_gain*n(1:length(x));

y = x + n;

%% Make window
win_t = 0.03;                   % window size in seconds
win_s = round(fs*win_t);        % window size in samples
if (mod(win_s, 2)==0)           % make odd
    win_s = win_s - 1;
end
win = hann(win_s);
% normalize it so the power is equal to its length
win = win*sqrt(length(win)/sum(win.^2));

%% STFT of signal
hop_size = (win_s-1)/2; % hop size (half of window size)

num_frames = floor(length(y)/hop_size);
nfft = 8*win_s;   % over sample to prevent time aliasing of filters

data.est_Sx = zeros(nfft, num_frames); % estimate of clean speech spectrum
data.est_Pn = zeros(nfft, num_frames); % estimate of noise power spectrum
data.est_Mn = zeros(nfft, num_frames); % estimate of noise magnitude spectrum

for i = 1: num_frames
    data.iteration = i;
    %% FFT
    s = (i-1)*hop_size + 1;                 % start index
    e = min(s+win_s-1, length(y));          % end index    
    if(mod(e-s+1,2)==0) e=e-1; end          % make length odd
    
    l = e-s+1; % length of windowed signal
    
    yzp = zpzpwin(y(s:e), win(1:l), nfft);    % zero-pad, zero-phase windowing
    xzp = zpzpwin(x(s:e), win(1:l), nfft);   
    nzp = zpzpwin(n(s:e), win(1:l), nfft);
    
    data.Sy = fft(yzp); % FFT
    data.Sx = fft(xzp); % Sx, and Sn for analysis purposes only!!
    data.Sn = fft(nzp);
    
    %% Estimate power spectrum of noise
    data = noiseEstimationSNR(data);
    %% Do magic!!!
    [data, est_Sx] = enhanceSpeechFilled(data);
    data.est_Sx(:, i) = est_Sx; % for speed, set data element outside the function    
end

%% Overlap add
xhat = real(invmyspectrogram2(data.est_Sx, hop_size, win_s));
dc = sum(win(1:hop_size:end)); 
xhat = xhat/dc;

xhat = xhat(1:length(y));

%% If you want to hear it in MATLAB
% sound(y, fs);
% sound(xhat, fs);