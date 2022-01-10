function xhat=compute_output(reconstruction,TF,options_tf, V_m, V_m_pha, W_1, W_2, W_3, W_4, H_1_m, H_2_m, H_3_m, H_4_m)
%% Purpose
% Compute separated heart, lung, noise supervised and noise unsupervised
% sounds
%% Inputs
% reconstruction= method of reconstructing separated signals in time domain
% TF= time-frequency representation
% options_tf= parameter settings related to time-frequency representation
% V_m= desired positive matrix with all recording examples
% V_m_pha= phase of all recording examples in V_mr t
% W_1, W_2, W_3, W_4= base matrices for heart, noise unsupervised, noise
% and lung respectively
% H_1_m, H_2_m, H_3_m, H_4_m= activation matrices for the recording
% examples V_m for the related base matrices. 
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
            end
        case "Filtering"
            % get the mixture phase
            phi = V_m_pha{idx};
            
            W=[W_1,W_2,W_3,W_4];
            H=[H_1_m{idx};H_2_m{idx};H_3_m{idx};H_4_m{idx}];
            WH = W * H;
            
            Mask_1= (W_1 * H_1_m{idx}) ./ WH;
            Mask_2= (W_2 * H_2_m{idx}) ./ WH;
            Mask_3= (W_3 * H_3_m{idx}) ./ WH;
            Mask_4= (W_4 * H_4_m{idx}) ./ WH;
            
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
                case "Cochleagram"
                    xhat(:,1,idx) = synthesisFast(xm{idx},Mask_1,options_tf.fRange,options_tf.cgL,fs);
                    xhat(:,2,idx) = synthesisFast(xm{idx},Mask_2,options_tf.fRange,options_tf.cgL,fs);
                    xhat(:,3,idx) = synthesisFast(xm{idx},Mask_3,options_tf.fRange,options_tf.cgL,fs);
                    xhat(:,4,idx) = synthesisFast(xm{idx},Mask_4,options_tf.fRange,options_tf.cgL,fs);
            end
        case "Best Mask"
            % get the mixture phase
            phi = V_m_pha{idx};
            
            Rec(:,:,1)=W_1 * H_1_m{idx};
            Rec(:,:,2)=W_2 * H_2_m{idx};
            Rec(:,:,3)=W_3 * H_3_m{idx};
            Rec(:,:,4)=W_4 * H_4_m{idx};
            
            for i=1:4
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