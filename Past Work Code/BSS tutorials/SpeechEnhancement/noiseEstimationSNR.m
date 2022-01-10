function [data] = noiseEstimationSNR(data)
% SNR based noise estimation method. Update noise estimate if SNR in given
% frequency is low, and don't update if SNR is high
% Input: 
%   data.Sy: Fourier transform coefficients of the noisy spectrum
% Output:
%   data.est_Mn: Estimated magnitude spectrum of noise
%   data.est_Pn: Estimated power spectrum of noise
%   data.last_Mns: Most recent 10 samples of data.est_Mn
% Reference: 
%   Lin, L., Holmes, W., and Ambikairajah, E. Adaptive noise estimation
%   algorithm for speech enhancement.

beta = 0.5;
i = data.iteration;

if (i <= 10) % Just use noisey power spectra for first few  iterations
    data.est_Mn = abs(data.Sy);
    data.est_Pn = data.est_Mn.^2;
else
    gammak = abs(data.Sy).^2./mean(data.last_Mns.^2, 2);
    alpha = 1./(1+exp(-beta*(gammak-1.5)));

    %% Set new estimate
    data.est_Mn = alpha.*data.last_Mns(:, end) + (1-alpha).*abs(data.Sy);
    data.est_Pn = alpha.*data.last_Mns(:, end).^2 + (1-alpha).*abs(data.Sy).^2;
end

%% Keep track of last 10 estimates as an estimate for noise power spectra
if (i == 1)
    data.last_Mns = data.est_Mn;
elseif (i <= 10)
    data.last_Mns = [data.last_Mns, data.est_Mn];
else
    data.last_Mns = [data.last_Mns(:, end-9:end), data.est_Mn];
end
