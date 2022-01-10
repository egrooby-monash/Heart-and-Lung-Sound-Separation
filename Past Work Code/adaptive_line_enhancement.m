function [heart,lung,difference]=adaptive_line_enhancement(x,mu,L,Delta,iterations,method)
%% Paper Information
% Separation of heart sound signal from lung sound signal by adaptive line enhancement
% https://ieeexplore.ieee.org/document/7099001
%% Purpose
% To separate heart sound from lung sound/white noise
% Adaptive line enhancement works by relying on the fact one signal is
% periodic and holds a strong correlation with the delayed version of the
% signal, while the other does not. 
%% Inputs
% x= mixture
% mu= learning rate of adaptive filter
% L= FIR filter length
% Delta= delay of reference desired signal
% iterations= iterations of the adaptive filter to stablise
% method= method of adaptive filter e.g. RLS, LMS etc. 
%% Output
% Separated heart and lung sounds
% difference= change in outputs in each iteration

if nargin<2
    % For heart and lung separation these were the parameters
    mu=0.0001;
    L=256;
    Delta = 375;
    % made up these ones
    iterations=50;
    method='Normalized LMS'; 
end

% generate adaptive filter
filt = dsp.LMSFilter('Method',method,'Length',L,'StepSize',mu);

% delayed version of input signal
delay = dsp.Delay(Delta);
d = delay(x);

% run through adaptive filter for numerous interataions
heart_prev=zeros(size(x));
lung_prev=zeros(size(x)); 
difference=zeros(iterations,1); 
for i=1:iterations
    desired=x;
    input=d;
    [output,error]=filt(input,desired);
    heart=output;
    lung=error;
    difference(i)=sum(abs(heart-heart_prev))+sum(abs(lung-lung_prev));
    heart_prev=heart;
    lung_prev=lung;
end
