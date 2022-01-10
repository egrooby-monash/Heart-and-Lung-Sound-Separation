function [data, est_Sx] = enhanceSpeech(data)
% Input: 
%  - data.Sy: Noisy speech spectrum. Y(w).
%  - data.est_Mn: Estimate of noise magnitude spectrum. E[|D(w)|]. 
%  - data.est_Pn: Estimate of noise power spectrum. E[|D(w)|^2]. 
%  - data.iteration: Frame index
% Output
%  - est_Sx: Estimate of clean speech spectrum. X(w).

switch (data.denoise_type)
    case 'spec_sub_mag'
        %% Magnitude spectral subtraction
    
    case 'spec_sub_power'
        %% Power spectral subtraction    

    case 'spec_sub_over'
        %% Oversubtraction method

    case 'wiener' 
        %% Wiener filter

    case 'mmse_stsa' 
        %% Ephraim & Malah

end

end