function [data, est_Sx] = enhanceSpeechFilled(data)
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
        estMagX = max(abs(data.Sy) - data.est_Mn, 1e-10);
        estAngleX = angle(data.Sy); 
        
        est_Sx = estMagX.*exp(complex(0,1)*estAngleX);        
    case 'spec_sub_power'
        %% Power spectral subtraction    
        estPowX = max(abs(data.Sy).^2 - data.est_Pn, 1e-10);
        estAngleX = angle(data.Sy);
        
        est_Sx = sqrt(estPowX).*exp(complex(0,1)*estAngleX);
    case 'spec_sub_over'
        %% Oversubtraction method
        % 1) fixed.
        alpha = 15;
        beta = 0.01;
        % 2) SNR varied.
%         snr = 10*log10(sum(abs(data.Sy).^2)./sum(data.est_Pn));
%         if (snr >= -5 && snr <= 20)
%             alpha = -3/20*snr + 4; 
%         elseif (snr > 20)
%             alpha = 1;
%         else
%             alpha = 4.75;
%         end
%         beta = 0.002;
        
        estPowX = max(abs(data.Sy).^2 - alpha*data.est_Pn, beta*data.est_Pn);
        estAngleX = angle(data.Sy);
        
        est_Sx = sqrt(estPowX).*exp(complex(0,1)*estAngleX);
    case 'wiener'
        gamma = abs(data.Sy).^2./data.est_Pn;
        
        % using inst_snr
        inst_snr = max(gamma - 1, 1e-10);
        
        % using decision directed approach
        if (data.iteration == 1)
            ksi = inst_snr;
        else
            alpha = 0.98;
            ksi = alpha*abs(data.prevEstSx).^2./data.prevEstPn + (1-alpha)*inst_snr;
        end
        
        H = ksi./(1+ksi);
        
        est_Sx = H.*data.Sy;
        
        % save previous estimates
        data.prevEstSx = est_Sx;
        data.prevEstPn = data.est_Pn;
    case 'mmse_stsa'
        gamma = abs(data.Sy).^2./data.est_Pn;
        gamma = min(gamma, 100);
        
        % using inst_snr
        inst_snr = max(gamma - 1, 1e-10);
        
        % using decision directed approach
        if (data.iteration == 1)
            ksi = inst_snr;
        else
            alpha = 0.98;
            ksi = alpha*abs(data.prevEstSx).^2./data.prevEstPn + (1-alpha)*inst_snr;
        end
        
        nu = ksi./(1+ksi).*gamma;
        H = sqrt(pi)/2*sqrt(nu)./gamma.*exp(-nu/2).*((1+nu).*besseli(0,nu/2) + nu.*besseli(1,nu/2));
        
        est_Sx = H.*data.Sy;
        
        % save previous estimates
        data.prevEstSx = est_Sx;
        data.prevEstPn = data.est_Pn;
end

end