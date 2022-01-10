function [heartSound,lungSound]=adf_filtering(x,fs,heart_reference,level)
%% Paper Information
% Adaptive Fourier decomposition approach for lung-heart sound separation
% https://ieeexplore.ieee.org/abstract/document/7158631
%% Purpose
% Separation of heart and lung sounds
%% Inputs
% x=mixture recording
% fs=sampling frequency 
% heart_reference= method for obtaining S1 and S2 peaks
% level= level of decomposition to reconstruction heart sounds
%% Outputs
% Separated heart and lung sounds 

if nargin<3
    heart_reference='springer';
    level=15; 
end

% 50-250Hz bandpass filter
filter_length=400;
filter_band=[50 250];
b = fir1(filter_length,filter_band/(fs/2));
bandpassed= filtfilt(b,1,x);

% S1 and S2 peak detection
heart_locations=heart_localisation(bandpassed,fs,heart_reference);

change_states=diff(heart_locations);
start_heart=find(change_states==1);

end_heart=find(change_states==-1);

if length(start_heart)>length(end_heart)
    end_heart=[end_heart;length(heart_locations)];
elseif length(end_heart)>length(start_heart)
    start_heart=[1;start_heart]; 
end
relevant_sections=end_heart-start_heart>10; 
start_heart(~relevant_sections)=[]; 
end_heart(~relevant_sections)=[];
heartSound=zeros(size(x)); 
for i=1:length(start_heart)
    segment=x(start_heart(i):end_heart(i)); 
    % Single Channel Fast AFD, circle searching dictionary, same phase (0~2\pi)
    % init AFD computation module
    afdcal_core=AFDCal();
    afdcal_core.setInputSignal(segment');
    % set decomposition method: Single Channel Fast AFD
    afdcal_core.setDecompMethod(2);
    % set AFD method: core AFD
    afdcal_core.setAFDMethod(1);
    % generate searching dictionary
    afdcal_core.genDic(0.02,0.95);
    % generate evaluators
    afdcal_core.genEva();
    % initilize decomposition
    afdcal_core.init_decomp()
    % decomposition 
    for n=1:level
        afdcal_core.nextDecomp()
    end
    reSig = afdcal_core.cal_reSig(level);
    heartSound(start_heart(i):end_heart(i))=reSig'; 
end

lungSound=x-heartSound;
end

