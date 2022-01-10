%% INIT
clear all; close all; clc;

% Read in audio file
x1  = wavread('bass');
x2  = wavread('drums');
[xm fs] = wavread('drums+bass');

K = [25 25]; % number of basis vectors
supervised = [1 1]; % binary vector specifying what is supervised
MAXITER = 500; % total number of iterations to run

 
%% STFT
FFTSIZE = 1024;
HOPSIZE = 256;
WINDOWSIZE = 512;

X1 = myspectrogram(x1,FFTSIZE,fs,hann(WINDOWSIZE),-HOPSIZE);
V1 = abs(X1(1:(FFTSIZE/2+1),:));
X2 = myspectrogram(x2,FFTSIZE,fs,hann(WINDOWSIZE),-HOPSIZE);
V2 = abs(X2(1:(FFTSIZE/2+1),:));
Xm = myspectrogram(xm,FFTSIZE,fs,hann(WINDOWSIZE),-HOPSIZE);
Vm = abs(Xm(1:(FFTSIZE/2+1),:)); maxV = max(max(db(Vm)));

F = size(Vm,1);
T = size(Vm,2);

figure;
subplot(3,1,1)
imagesc(db(V1))
set(gca,'YDir','normal')
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('Drum + Bass')
ylabel('Frequency')
xlabel('Time')


subplot(3,1,2)
imagesc(db(V2))
set(gca,'YDir','normal')
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('Drum + Bass')
ylabel('Frequency')
xlabel('Time')


subplot(3,1,3)
imagesc(db(Vm))
set(gca,'YDir','normal')
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('Drum + Bass')
ylabel('Frequency')
xlabel('Time')


%% TODO: NMF


 
%% ISTFT / RECONSTRUCTION METHOD 2 (FILTERING)


% get the mixture phase
phi = angle(Xm);
c = [1 cumsum(K)];
for i=1:length(K)
    
    % create masking filter
    Mask =  W(:,c(i):c(i+1))*H(c(i):c(i+1),:)./(W*H);
    
    % filter
    XmagHat = Vm.*Mask; 
    
    % create upper half of frequency before istft
    XmagHat = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
   
    % Multiply with phase
    XHat = XmagHat.*exp(1i*phi);
    
    % create upper half of frequency before istft
    xhat(:,i) = real(invmyspectrogram(XmagHat.*exp(1i*phi),HOPSIZE))';
    
%     sound(xhat(:,i),fs) 
    
end
 

%% PLOT BASIS Vectors
freq = linspace(0, fs/2, FFTSIZE/2+1);
time = linspace(0, length(xm)/fs, T); 

figure;
for i=1:sum(K)
    plot((i-1)*max(max(W))+(1-W(:,i)),freq,'LineWidth', 3)
    hold on 
end
title('Basis Vectors')
ylabel('Frequency (Hz)')
xlabel('Basis')
ylim([0 1000])
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);


%% PLOT ACTIVATIONS

figure;
for i=1:sum(K)
    plot(time, (i-1)*max(max(H))+(H(i,:)),'LineWidth', 3)
    hold on
end
ylabel('Activations')
xlabel('Time (seconds)')

set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);

%% PLOT SOURCE ACTIVATIONS

figure;
c = [1 cumsum(K)];
for i=1:length(K)
    
    plot(time, (i-1)*max(max(sum(H)))+sum(H(c(i):c(i+1),:)),'LineWidth', 3)
    hold on
end
ylabel('Source Activations')
xlabel('Time (seconds)')
set(gca,'YDir','normal')
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
set(gca, 'YTick', []);
set(gca, 'XTick', []);

set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);



%% PLOT RECONSTRUCTION VS. ORIGINAL

figure;
subplot(2,1,1)
imagesc(db(Vm))
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
c = [1 cumsum(K)];
for i=1:length(K)
    
    subplot(length(K),1,i)
    R = (max( db(W(:,c(i):c(i+1))*H(c(i):c(i+1),:)) - m, -CLIP) + CLIP)/CLIP;
    imagesc(R)
    set(gca,'YDir','normal')
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
    ylabel('Frequency')
    xlabel('Time')
    title(['Layer from Source ' int2str(i) ])

end


%% PLOT MASKS

figure;
CLIP = 100;
m = max(max(db(W*H)));
c = [1 cumsum(K)];
for i=1:length(K)
    
    subplot(length(K),1,i)
    R = (max( db(W(:,c(i):c(i+1))*H(c(i):c(i+1),:)./(W*H)) - m, -CLIP) + CLIP)/CLIP;
    imagesc(R)
    set(gca,'YDir','normal')
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    set(gca, 'YTick', []);
    set(gca, 'XTick', []);
    ylabel('Frequency')
    xlabel('Time')
    title(['Mask for Source ' int2str(i) ])

end

