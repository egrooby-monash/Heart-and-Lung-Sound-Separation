function [heartSound,lungSound]=frequency_interpolation(x,fs,heart_reference,remove_method,interpolation_method)
%% Paper Information
% Heart sound cancellation from lung sound recordings using adaptive threshold and 2D interpolation in time-frequency domain
% https://ieeexplore.ieee.org/abstract/document/1280444 
% Heart sound cancellation from lung sound recordings using time-frequency filtering
% https://link.springer.com/article/10.1007/s11517-006-0030-8 
%% Purpose
% To remove heart sound interference from lung sound recordings
%% Inputs
% x= mixture
% fs= sampling frequency
% heart_reference= method for finding location of S1 and S2 heart sounds
% remove_method= how to remove heart sounds from time-frequency domain
% interpolation_method= method of 2D interpolation and extrapolation 

%% Outputs

if nargin<3
    heart_reference='springer';
    remove_method='all';
    % interpolation method: 'linear' , 'nearest', 'natural'
    interpolation_method='linear'; 
end 

switch heart_reference
    case 'springer'
        % bandpass filter 50-250Hz
        filter_length=400; 
        filter_band=[50 250];
        b = fir1(filter_length,filter_band/(fs/2));
        bandpassed= filtfilt(b,1,x);
        % locate S1 and S2 heart sounds
        options.env='hilbert';
        options.autocorr='filtered';
        options.init_hr='autocorr_peak';
        options.systolic='yes';
        options.seg='springer';
        load('Springer_B_matrix.mat', 'Springer_B_matrix');
        load('Springer_pi_vector.mat', 'Springer_pi_vector');
        load('Springer_total_obs_distribution.mat', 'Springer_total_obs_distribution');
        options.seg_fs=50;
        options.seg_pi_vector=Springer_pi_vector;
        options.seg_b_matrix=Springer_B_matrix;
        options.seg_total_obs_dist=Springer_total_obs_distribution;
        max_HR=220;
        min_HR=70; 
        heart= get_hr_segmentation(bandpassed, fs, max_HR,min_HR,options);
        assigned_states_pks=heart.seg_states{1}; 
        assigned_states_pks = round(resample(assigned_states_pks,fs,1000)); 
        % bandpass signal with just S1 and S2 peaks
        bandpassed(assigned_states_pks==2)=0;
        bandpassed(assigned_states_pks==4)=0;
        heartSound=bandpassed; 
        % convert bandpass filter into same dimensions in time-frequency
        % domain
        [s2,~] =stft(bandpassed,fs,'Window',hann(fs*0.1,'periodic'),'OverlapLength',fs*0.05);
        power=sum(abs(s2)); 
        remove_sections=power>mean(power); 
end

% time-frequency domain
[s,f] =stft(x,fs,'Window',hann(fs*0.1,'periodic'),'OverlapLength',fs*0.05);
remove_sections=remove_sections(1:size(s,2));

% how to remove heart sounds
switch remove_method
    case 'section'
        % Remove 20-300Hz
        section=(f>=20 & f<=300)| (f<=-20 & f>=-300); 
        s(section,remove_sections)=NaN; 
    case 'all'
        % Remove entire section
        s(:,remove_sections)=NaN; 
end
    
% 2D interpolation
% returns interpolated values of a function of two variables at specific
% query points using linear interpolation. The results always pass through
% the original sampling of the function. X and Y contain the coordinates of
% the sample points. V contains the corresponding function values at each
% sample point. Xq and Yq contain the coordinates of the query points.
V=s(~isnan(s));
Xq=repmat(1:size(s,2),size(s,1),1); 
Yq=repmat((1:size(s,1))',1,size(s,2));
X=Xq(~isnan(s));
Y=Yq(~isnan(s));


% meshgrid format of 2D interpolation
% V=s(:,~remove_sections);
% Xq=repmat(1:size(s,2),size(s,1),1); 
% Yq=repmat((1:size(s,1))',1,size(s,2));
% X=Xq(:,~remove_sections);
% Y=Yq(:,~remove_sections);
%Vq = interp2(X,Y,V,Xq,Yq,'spline');
%'nearest' 'linear' 'spline'

% scatter version of 2D interpolation
% Vq = griddata(X(:),Y(:),V(:),Xq(:),Yq(:),'cubic');
% Vq=reshape(Vq,size(s));
% Vq=Vq(:,1:find(~remove_sections,1,'last')); 

% scatter version of 2D interpolation that also does extrapolation
F = scatteredInterpolant(X(:),Y(:),V(:),interpolation_method,interpolation_method);
Vq = F(Xq(:),Yq(:));
Vq=reshape(Vq,size(s));

% time domain
lungSound=real(istft(Vq,'Window',hann(fs*0.1,'periodic'),'OverlapLength',fs*0.05));
end



