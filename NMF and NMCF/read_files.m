function [V,V_pha]=read_files(path,TF,options_tf,V,V_pha,max_examples)
%% Purpose
% Read all files and calculation magnitude and phase
%% Inputs
% path= folder location of files to calculation time frequency
% respresentaiton
% TF= time frequency representation method
% options_tf= parameters for time frequency representaiton
% V and V_pha are magnitude and phase of time frequency representation
% max_examples= max number of files to read
%% Output
% V and V_pha are magnitude and phase of time frequency representation

count=length(V);
subcount=1;
if nargin<2
    TF="STFT";
    switch TF
        case "STFT"
            options_tf.FFTSIZE = 1024;
            options_tf.HOPSIZE = 256;
            options_tf.WINDOWSIZE = 512;
        case "Cochleagram"
            % fRange= frequency range
            options_tf.fRange = [50, fs/2];
            % num_bin= number of bins/channels for gammatone filterbank
            options_tf.num_bin = 128;
            % filterOrder= filter order of Gammatone filterbank
            options_tf.filterOrder = 4;
            % gammatone filter length or 128 ms for 16 kHz sampling rate
            options_tf.gL = round(0.128*fs);
            % default window length in sample points which is 20 ms for 16 KHz sampling frequency
            options_tf.cgL=0.02*fs;
        case 'Q-transform'
            % Number of bins per octave
            options_tf.bins=24;
        case 'Q-transform2'
            options_tf.bins=64;
            options_tf.total_bin=380;
            options_tf.sparsity=0.01;
            options_tf.HOPSIZE=256;
        case 'STFT2'
            options_tf.FFTSIZE = 1024;
            options_tf.HOPSIZE = 256;
            options_tf.WINDOWSIZE = 512;
    end
end

% loading all files
files = dir(path);
for i=1:length(files)
    % skip irrelevant file names
    if strcmp(files(i).name,'.')
        continue
    elseif strcmp(files(i).name,'..')
        continue
    elseif contains(files(i).name,'.png')
        continue
    elseif strcmp(files(i).name,'.DS_Store')
        continue
    end
    % read file
    [xm,fs]=audioread(sprintf('%s/%s',path,files(i).name));
    % dealing with stethoscope movement
    if length(xm)<10*4000
        temp=zeros(10*4000,1);
        insert_point=round(rand(1)*(10*4000-length(xm)));
        end_insert_point=insert_point+length(xm);
        temp(insert_point+1:end_insert_point)=xm;
        xm=temp;
    end
    %% Time Frequency Representation
    switch TF
        case "STFT"
            Xm = myspectrogram_modified(xm,options_tf.FFTSIZE,fs,hann(options_tf.WINDOWSIZE),-options_tf.HOPSIZE);
            Vm = abs(Xm(1:(options_tf.FFTSIZE/2+1),:));
            F = size(Vm,1);
            T = size(Vm,2);
        case "Cochleagram"
            %% Paper Information
            % Unsupervised SIngle-Channel Separation of Non-Stationary Signals using
            % Gammatone Filterbank and Itakura Saito Non-negative Matrix
            % Two-Dimensional Factorization
            % https://www.researchgate.net/profile/Wai-Lok-Woo/publication/236007062_Unsupervised_Single-Channel_Separation_of_Nonstationary_Signals_Using_Gammatone_Filterbank_and_Itakura-Saito_Nonnegative_Matrix_Two-Dimensional_Factorizations/links/5efe1aed45851550508533a6/Unsupervised-Single-Channel-Separation-of-Nonstationary-Signals-Using-Gammatone-Filterbank-and-Itakura-Saito-Nonnegative-Matrix-Two-Dimensional-Factorizations.pdf
            % https://www.mathworks.com/matlabcentral/fileexchange/48622-cochleagram-and-is-nmf2d-for-blind-source-separation
            
            % gammatone filterbank (cochleagram) which is more separable than the
            % spectrogram or log-frequency spectrogram (constant Q-transform)
            % Construct the cochleagram use Gammatone filterbank
            gf = gammatoneFast_modified(xm,options_tf.num_bin,options_tf.fRange,fs,options_tf.filterOrder,options_tf.gL);
            Xm = cochleagram_modified(gf,options_tf.cgL);
            Vm=abs(Xm);
            F = size(Vm,1);
            T = size(Vm,2);
        case "Q-transform"
            [Xm,~,g,fshifts] = cqt(xm,'SamplingFrequency' ,fs,'BinsPerOctave',options_tf.bins);
            Vm = abs(Xm);
            F = size(Vm,1);
            T = size(Vm,2);
        case "Q-transform2"
            freq=py.librosa.cqt(py.numpy.array(xm),pyargs(...
                'sr',int32(fs),...
                'n_bins',int32(options_tf.total_bin),...
                'bins_per_octave',int32(options_tf.bins),...
                'sparsity',options_tf.sparsity,...
                'hop_length',int32(options_tf.HOPSIZE)));
            Xm=complex(double(py.numpy.array(freq.real)),double(py.numpy.array(freq.imag)));
            Vm = abs(Xm);
            F = size(Vm,1);
            T = size(Vm,2);
            
        case "STFT2"
            xm(40001:40192)=0;
            [Xm,f] = stft(xm,fs,'Window',hann(options_tf.WINDOWSIZE,'periodic'),...
                'OverlapLength',options_tf.WINDOWSIZE-options_tf.HOPSIZE,'FFTLength',...
                options_tf.FFTSIZE,'FrequencyRange','onesided');
            Vm = abs(Xm);
            F = size(Vm,1);
            T = size(Vm,2);
    end
    % storing information
    V{count+subcount}=Vm;
    V_pha{count+subcount}=angle(Xm);
    subcount=subcount+1;
    if subcount>max_examples
        return
    end
end
end