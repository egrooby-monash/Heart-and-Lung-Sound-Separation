%% INIT
clear all; close all; clc;

% Read in audio file
[x fs] = wavread('Mary');


%% STFT
FFTSIZE = 1024;
HOPSIZE = 256;
WINDOWSIZE = 512;

X = myspectrogram(x,FFTSIZE,fs,hann(WINDOWSIZE),-HOPSIZE);
V = abs(X(1:(FFTSIZE/2+1),:));
F = size(V,1);
T = size(V,2);

imagesc(db(V))
set(gca,'YDir','normal')
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('Spectrogram of Mary Had a Little Lamb')
ylabel('Frequency')
xlabel('Time')

%sound(x,fs) 

%% NMF

K = 3; % number of basis vectors
MAXITER = 100; % total number of iterations to run


%% TODO: NMF

[W, H] = nmf(V, K, MAXITER);


%% ISTFT / RECONSTRUCTION METHOD 1 (SYNTHESIS)

% get the mixture phase
phi = angle(X);

for i=1:K
    
    XmagHat = W(:,i)*H(i,:);
    
    % create upper half of frequency before istft
    XmagHat = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
   
    % Multiply with phase
    XHat = XmagHat.*exp(1i*phi);
    
    xhat1(:,i) = real(invmyspectrogram(XHat,HOPSIZE))';
    
end

% sound(xhat1(:,1),fs) 
% sound(xhat1(:,2),fs) 
% sound(xhat1(:,3),fs) 

%% ISTFT / RECONSTRUCTION METHOD 2 (FILTERING)


% get the mixture phase
phi = angle(X);

for i=1:K
    
    % create masking filter
    Mask =  (W(:,i)*H(i,:)./(W*H));
    
    % filter
    XmagHat = V.*Mask; 
    
    % create upper half of frequency before istft
    XmagHat = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
   
    % Multiply with phase
    XHat = XmagHat.*exp(1i*phi);
    
    % create upper half of frequency before istft
    xhat2(:,i) = real(invmyspectrogram(XmagHat.*exp(1i*phi),HOPSIZE))';
    
end

% sound(xhat1(:,1),fs) 
% sound(xhat1(:,2),fs) 
% sound(xhat1(:,3),fs) 



%% PLOT BASIS Vectors
freq = linspace(0, fs/2, FFTSIZE/2+1);
time = linspace(0, length(x)/fs, T); 
% color = ['r','g','b',''];


figure;
for i=1:K
    plot((i-1)*max(max(W))+(1-W(:,i)),freq,'LineWidth', 3)
    hold on 
end
title('Basis Vectors')
ylabel('Frequency (Hz)')
xlabel('Basis')
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);


%% PLOT ACTIVATIONS

figure;
for i=1:K
    plot(time, (i-1)*max(max(H))+(H(i,:)),'LineWidth', 3)
    hold on
end
ylabel('Activations')
xlabel('Time (seconds)')

set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);

%% PLOT RECONSTRUCTION VS. ORIGINAL

figure;
subplot(2,1,1)
imagesc(db(V))
set(gca,'YDir','normal')
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('Original')
ylabel('Frequency')
xlabel('Time')

subplot(2,1,2)
imagesc(db(W*H))
set(gca,'YDir','normal')
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('Reconstruction')
ylabel('Frequency')
xlabel('Time')


%% PLOT LAYERS

figure;
CLIP = 100;
m = max(max(db(W*H)));
for i=1:K
    
    subplot(K,1,i)
    R = (max( db(W(:,i)*H(i,:)) - m, -CLIP) + CLIP)/CLIP;
    imagesc(R)
    set(gca,'YDir','normal')
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
    ylabel('Frequency')
    xlabel('Time')
    title(['Layer from Basis Vector ' int2str(i) ])

end


%% PLOT MASKS

figure;
CLIP = 100;
m = max(max(db(W*H)));
for i=1:K
    
    subplot(K,1,i)
    R = (max( db(W(:,i)*H(i,:)./(W*H)) - m, -CLIP) + CLIP)/CLIP;
    imagesc(R)
    set(gca,'YDir','normal')
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
    ylabel('Frequency')
    xlabel('Time')
    title(['Mask for Basis Vector ' int2str(i) ])

end

