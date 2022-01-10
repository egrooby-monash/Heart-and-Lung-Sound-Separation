function heart_location=heart_localisation(x,fs,heart_reference)
%% Purpose
% Identify location of S1 and S2 peaks
%% Inputs
% x= mixture recording
% fs= sampling frequency
% heart_reference= method to identify S1 and S2 peaks
%% Output
% Vector of zeros with 1s at site of S1 and S2 heart sounds

switch heart_reference
    case 'springer'
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
        heart= get_hr_segmentation(x, fs, max_HR,min_HR,options);
        assigned_states_pks=heart.seg_states{1}; 
        assigned_states_pks = round(resample(assigned_states_pks,fs,1000)); 
        heart_location=zeros(size(assigned_states_pks)); 
        heart_location(assigned_states_pks==1)=1; 
        heart_location(assigned_states_pks==3)=1; 
    case 'entropy'
        [~,heart_location]=heart_localisation_entropy(x,fs);
    case 'log energy'
        [~,heart_location]=heart_localisation_log_energy(x,fs);
    case 'wavelet'
        [~, heart_location]=heart_localisation_wavelet(x);
end







