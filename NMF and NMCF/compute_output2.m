function xhat=compute_output2(reconstruction,TF,options_tf, V_m, V_m_pha, W_1, W_2, W_3, W_4, W_5, W_6, H_1_m, H_2_m, H_3_m, H_4_m, H_5_m, H_6_m,xm,g,fshifts,fs)
%% Purpose
% Compute separated heart, lung, noise supervised and noise unsupervised
% sounds
%% Inputs
% reconstruction= method of reconstructing separated signals in time domain
% TF= time-frequency representation
% options_tf= parameter settings related to time-frequency representation
% V_m= desired positive matrix with all recording examples
% V_m_pha= phase of all recording examples in V_mr t
% W_1, W_2, W_3, W_4, W_5, W_6= base matrices for heart, noise unsupervised, cry noise, 
% lung, stethoscope movement noise and respiratory support noise respectively
% H_1_m, H_2_m, H_3_m, H_4_m, H_5_m, H_6_m= activation matrices for the recording
% examples V_m for the related base matrices. 
% xm,g,fshifts= parameters for cochleagram time-frequency representation
% fs= sampling frequency
%% Outputs
% xhat= column matrix of separated sounds which contains heart, noise unsupervised, noise
% and lung respectively
for idx=1:length(V_m)
    %% ISTFT / RECONSTRUCTION
    switch reconstruction
        case "Synthesis"
            % get the mixture phase
            phi = V_m_pha{idx};
            switch TF
                case "STFT"
                    XmagHat = W_1*H_1_m{idx};
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,1,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                    
                    XmagHat = W_2*H_2_m{idx};
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,2,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                    
                    XmagHat = W_3*H_3_m{idx};
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,3,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                    
                    XmagHat = W_4*H_4_m{idx};
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,4,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                    
                    XmagHat = W_5*H_5_m{idx};
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,5,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                    
                    XmagHat = W_6*H_6_m{idx};
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,6,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                case "Q-transform"
                    XmagHat_full = W_1*H_1_m{idx};
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,1,idx) = real(icqt(XHat,g,fshifts))';
                    
                    XmagHat_full = W_2*H_2_m{idx};
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,2,idx) = real(icqt(XHat,g,fshifts))';
                    
                    XmagHat_full = W_3*H_3_m{idx};
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,3,idx) = real(icqt(XHat,g,fshifts))';
                    
                    XmagHat_full = W_4*H_4_m{idx};
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,4,idx) = real(icqt(XHat,g,fshifts))';
                    
                    XmagHat_full = W_5*H_5_m{idx};
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,5,idx) = real(icqt(XHat,g,fshifts))';
                    
                    XmagHat_full = W_6*H_6_m{idx};
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,6,idx) = real(icqt(XHat,g,fshifts))';
            end
        case "Filtering"
            % get the mixture phase
            phi = V_m_pha{idx};
            
            W=[W_1,W_2,W_3,W_4,W_5,W_6];
            H=[H_1_m{idx};H_2_m{idx};H_3_m{idx};H_4_m{idx};H_5_m{idx};H_6_m{idx}];
            WH = W * H;
            
            Mask_1= (W_1 * H_1_m{idx}) ./ WH;
            Mask_2= (W_2 * H_2_m{idx}) ./ WH;
            Mask_3= (W_3 * H_3_m{idx}) ./ WH;
            Mask_4= (W_4 * H_4_m{idx}) ./ WH;
            Mask_5= (W_5 * H_5_m{idx}) ./ WH;
            Mask_6= (W_6 * H_6_m{idx}) ./ WH;
            
            switch TF
                case "STFT"
                    % filter
                    XmagHat = V_m{idx}.*Mask_1;
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,1,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_2;
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,2,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_3;
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,3,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_4;
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,4,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_5;
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,5,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_6;
                    % create upper half of frequency before istft
                    XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,6,idx) = real(invmyspectrogram_modified(XHat,options_tf.HOPSIZE))';
                case "Q-transform"
                    % filter
                    XmagHat_full = V_m{idx}.*Mask_1;
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,1,idx) = real(icqt(XHat,g,fshifts))';
                    
                    % filter
                    XmagHat_full = V_m{idx}.*Mask_2;
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,2,idx) = real(icqt(XHat,g,fshifts))';
                    
                    % filter
                    XmagHat_full = V_m{idx}.*Mask_3;
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,3,idx) = real(icqt(XHat,g,fshifts))';
                    
                    % filter
                    XmagHat_full = V_m{idx}.*Mask_4;
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,4,idx) = real(icqt(XHat,g,fshifts))';
                    
                    % filter
                    XmagHat_full = V_m{idx}.*Mask_5;
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,5,idx) = real(icqt(XHat,g,fshifts))';
                    
                    % filter
                    XmagHat_full = V_m{idx}.*Mask_6;
                    % Multiply with phase
                    XHat = XmagHat_full.*exp(1i*phi);
                    xhat(:,6,idx) = real(icqt(XHat,g,fshifts))';
                case "Cochleagram"
                    xhat(:,1,idx) = synthesisFast(xm,Mask_1,options_tf.fRange,options_tf.cgL,fs);
                    xhat(:,2,idx) = synthesisFast(xm,Mask_2,options_tf.fRange,options_tf.cgL,fs);
                    xhat(:,3,idx) = synthesisFast(xm,Mask_3,options_tf.fRange,options_tf.cgL,fs);
                    xhat(:,4,idx) = synthesisFast(xm,Mask_4,options_tf.fRange,options_tf.cgL,fs);
                    xhat(:,5,idx) = synthesisFast(xm,Mask_5,options_tf.fRange,options_tf.cgL,fs);
                    xhat(:,6,idx) = synthesisFast(xm,Mask_6,options_tf.fRange,options_tf.cgL,fs);
                case "Q-transform2"
                    % filter
                    XmagHat = V_m{idx}.*Mask_1;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    XHat_py=py.numpy.array(real(XHat),pyargs('dtype','complex'));
                    XHat_py.imag=imag(XHat);
                    wave=py.librosa.icqt(XHat_py,pyargs(...
                        'sr',int32(fs),...
                        'bins_per_octave',int32(options_tf.bins),...
                        'sparsity',options_tf.sparsity,...
                        'hop_length',int32(options_tf.HOPSIZE)));
                    xhat(:,1,idx)=double(wave);
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_2;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    XHat_py=py.numpy.array(real(XHat),pyargs('dtype','complex'));
                    XHat_py.imag=imag(XHat);
                    wave=py.librosa.icqt(XHat_py,pyargs(...
                        'sr',int32(fs),...
                        'bins_per_octave',int32(options_tf.bins),...
                        'sparsity',options_tf.sparsity,...
                        'hop_length',int32(options_tf.HOPSIZE)));
                    xhat(:,2,idx)=double(wave);
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_3;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    XHat_py=py.numpy.array(real(XHat),pyargs('dtype','complex'));
                    XHat_py.imag=imag(XHat);
                    wave=py.librosa.icqt(XHat_py,pyargs(...
                        'sr',int32(fs),...
                        'bins_per_octave',int32(options_tf.bins),...
                        'sparsity',options_tf.sparsity,...
                        'hop_length',int32(options_tf.HOPSIZE)));
                    xhat(:,3,idx)=double(wave);
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_4;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    XHat_py=py.numpy.array(real(XHat),pyargs('dtype','complex'));
                    XHat_py.imag=imag(XHat);
                    wave=py.librosa.icqt(XHat_py,pyargs(...
                        'sr',int32(fs),...
                        'bins_per_octave',int32(options_tf.bins),...
                        'sparsity',options_tf.sparsity,...
                        'hop_length',int32(options_tf.HOPSIZE)));
                    xhat(:,4,idx)=double(wave);
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_5;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    XHat_py=py.numpy.array(real(XHat),pyargs('dtype','complex'));
                    XHat_py.imag=imag(XHat);
                    wave=py.librosa.icqt(XHat_py,pyargs(...
                        'sr',int32(fs),...
                        'bins_per_octave',int32(options_tf.bins),...
                        'sparsity',options_tf.sparsity,...
                        'hop_length',int32(options_tf.HOPSIZE)));
                    xhat(:,5,idx)=double(wave);
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_6;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    XHat_py=py.numpy.array(real(XHat),pyargs('dtype','complex'));
                    XHat_py.imag=imag(XHat);
                    wave=py.librosa.icqt(XHat_py,pyargs(...
                        'sr',int32(fs),...
                        'bins_per_octave',int32(options_tf.bins),...
                        'sparsity',options_tf.sparsity,...
                        'hop_length',int32(options_tf.HOPSIZE)));
                    xhat(:,6,idx)=double(wave);
                case "STFT2"
                    % filter
                    XmagHat = V_m{idx}.*Mask_1;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    xhat(:,1,idx) = real(istft(XHat,fs,'Window',hann(options_tf.WINDOWSIZE,'periodic'),...
                        'OverlapLength',options_tf.WINDOWSIZE-options_tf.HOPSIZE,'FFTLength',...
                        options_tf.FFTSIZE,'FrequencyRange','onesided'));
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_2;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    xhat(:,2,idx) = real(istft(XHat,fs,'Window',hann(options_tf.WINDOWSIZE,'periodic'),...
                        'OverlapLength',options_tf.WINDOWSIZE-options_tf.HOPSIZE,'FFTLength',...
                        options_tf.FFTSIZE,'FrequencyRange','onesided'));
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_3;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    xhat(:,3,idx) = real(istft(XHat,fs,'Window',hann(options_tf.WINDOWSIZE,'periodic'),...
                        'OverlapLength',options_tf.WINDOWSIZE-options_tf.HOPSIZE,'FFTLength',...
                        options_tf.FFTSIZE,'FrequencyRange','onesided'));
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_4;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    xhat(:,4,idx) = real(istft(XHat,fs,'Window',hann(options_tf.WINDOWSIZE,'periodic'),...
                        'OverlapLength',options_tf.WINDOWSIZE-options_tf.HOPSIZE,'FFTLength',...
                        options_tf.FFTSIZE,'FrequencyRange','onesided'));
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_5;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    xhat(:,5,idx) = real(istft(XHat,fs,'Window',hann(options_tf.WINDOWSIZE,'periodic'),...
                        'OverlapLength',options_tf.WINDOWSIZE-options_tf.HOPSIZE,'FFTLength',...
                        options_tf.FFTSIZE,'FrequencyRange','onesided'));
                    
                    % filter
                    XmagHat = V_m{idx}.*Mask_6;
                    % Multiply with phase
                    XHat = XmagHat.*exp(1i*phi);
                    xhat(:,6,idx) = real(istft(XHat,fs,'Window',hann(options_tf.WINDOWSIZE,'periodic'),...
                        'OverlapLength',options_tf.WINDOWSIZE-options_tf.HOPSIZE,'FFTLength',...
                        options_tf.FFTSIZE,'FrequencyRange','onesided'));
            end
        case "Best Mask"
            
            % get the mixture phase
            phi = V_m_pha{idx};
            
            Rec(:,:,1)=W_1 * H_1_m{idx};
            Rec(:,:,2)=W_2 * H_2_m{idx};
            Rec(:,:,3)=W_3 * H_3_m{idx};
            Rec(:,:,4)=W_4 * H_4_m{idx};
            Rec(:,:,5)=W_5 * H_5_m{idx};
            Rec(:,:,6)=W_6 * H_6_m{idx};
            
            for i=1:6
                % create masking filter
                Mask = logical(Rec(:,:,i)==max(Rec,[],3));
                
                switch TF
                    case "STFT"
                        % filter
                        XmagHat = V_m{idx}.*Mask;
                        
                        % create upper half of frequency before istft
                        XmagHat_full = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                        
                        % Multiply with phase
                        XHat = XmagHat_full.*exp(1i*phi);
                        
                        xhat(:,i,idx) = real(invmyspectrogram(XHat,options_tf.HOPSIZE))';
                    case "Q-transform"
                        % filter
                        XmagHat_full = V_m{idx}.*Mask;
                        
                        % Multiply with phase
                        XHat = XmagHat_full.*exp(1i*phi);
                        
                        xhat(:,i,idx) = real(icqt(XHat,g,fshifts))';
                    case "Cochleagram"
                        xhat(:,i,idx) = synthesisFast(xm{idx},Mask,options_tf.fRange,options_tf.cgL,fs)';
                end
            end
    end
end
end