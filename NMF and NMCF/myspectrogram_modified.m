function X = myspectrogram_modified(x,nfft,fs,window,noverlap,doplot,dbdown)
%% Paper Information
% Single-Channel Source Separation Tutorial Mini-Series
% https://ccrma.stanford.edu/~njb/teaching/sstutorial/ 
%% Purpose
%MYSPECTROGRAM Calculate spectrogram from signal.
%% Inputs
% NFFT is the FFT size used for each frame of A.  It should be a
% power of 2 for fastest computation of the spectrogram.
%
% Fs is the sampling frequency. Since all processing parameters are
% in units of samples, Fs does not effect the spectrogram itself,
% but it is used for axis scaling in the plot produced when
% MYSPECTROGRAM is called with no output argument (see below).
%
% WINDOW is the length M window function applied, IN ZERO-PHASE
% FORM, to each frame of A.  M cannot exceed NFFT.  For M<NFFT,
% NFFT-M zeros are inserted in the FFT buffer (for interpolated
% zero-phase processing).  The window should be supplied in CAUSAL
% FORM.
%
% NOVERLAP is the number of samples the sections of A overlap, if
% nonnegative.  If negative, -NOVERLAP is the "hop size", i.e., the
% number of samples to advance successive windows.  (The overlap is
% the window length minus the hop size.)  The hop size is called
% NHOP below.  NOVERLAP must be less than M.
%
% If doplot is nonzero, or if there is no output argument, the 
% spectrogram is displayed.
%
% When the spectrogram is displayed, it is ``clipped'' dbdown dB
% below its maximum magnitude.  The default clipping level is 100 
% dB down.
%
% Thus, MYSPECTROGRAM splits the signal into overlapping segments of
% length M, windows each segment with the length M WINDOW vector, in
% zero-phase form, and forms the columns of B with their zero-padded,
% length NFFT discrete Fourier transforms.
%
% With no output argument B, MYSPECTROGRAM plots the dB magnitude of
% the spectrogram in the current figure, using
% IMAGESC(T,F,20*log10(ABS(B))), AXIS XY, COLORMAP(JET) so the low
% frequency content of the first portion of the signal is displayed
% in the lower left corner of the axes.
%
% Each column of B contains an estimate of the short-term,
% time-localized frequency content of the signal A.  Time increases
% linearly across the columns of B, from left to right.  Frequency
% increases linearly down the rows, starting at 0.
%
% If A is a length NX complex signal, B is returned as a complex
% matrix with NFFT rows and
%      k = floor((NX-NOVERLAP)/(length(WINDOW)-NOVERLAP)) 
%        = floor((NX-NOVERLAP)/NHOP)
% columns.  When A is real, only the NFFT/2+1 rows are needed when
% NFFT even, and the first (NFFT+1)/2 rows are sufficient for
% inversion when NFFT is odd.
%
% See also: Matlab's SPECTROGRAM and Octave's STFT function.
% 02/04/02/jos: Created
% 02/12/04/jos: Added dbdown
% 07/23/08/jos: Changed name from SPECTROGRAM to MYSPECTROGRAM
%% Output
% X= spectrogram for the signal in vector x
% B = MYSPECTROGRAM(A,NFFT,Fs,WINDOW,NOVERLAP) calculates the 
%     spectrogram for the signal in vector A.  


%% Default Values
if nargin<7
    dbdown=100; 
end
if nargin<6
    doplot=0; 
end
if nargin<5
    noverlap=256; 
end
if nargin<4
    window=hamming(512);
end
if nargin<3
    fs=1; 
end
if nargin<2
    nfft=2048; 
end

% make sure it's a column
x = x(:);

M = length(window);
if (M<2) 
    error('myspectrogram: Expect complete window, not just its length'); 
end
if (M<2) 
    error('myspectrogram: Expect complete window, not just its length'); 
end
% zero-pad to fill a window:
if length(x)<M 
  x = [x;zeros(M-length(x),1)]; 
end
 % 0 if M even, 1 if odd
Modd = mod(M,2);
Mo2 = (M-Modd)/2;
% Make sure it's a column
w = window(:); 

if noverlap<0
  nhop = - noverlap;
  noverlap = M-nhop;
else
  nhop = M-noverlap;
end

nx = length(x);
%nframes = 1+floor((nx-noverlap)/nhop);
nframes = 1+ceil(nx/nhop);

% allocate output spectrogram
X = zeros(nfft,nframes); 

% zero-padding for each FFT
zp = zeros(nfft-M,1); 
xframe = zeros(M,1);
% input time offset = half a frame
xoff = 0 - Mo2; 
for m=1:nframes
%  M,Mo2,xoff,nhop
  if xoff<0
      % partial input data frame
    xframe(1:xoff+M) = x(1:xoff+M); 
  else
    if xoff+M > nx
      xframe = [x(xoff+1:nx);zeros(xoff+M-nx,1)];
    else
        % input data frame
      xframe = x(xoff+1:xoff+M); 
    end
  end
  % Apply window
  xw = w .* xframe; 
  xwzp = [xw(Mo2+1:M);zp;xw(1:Mo2)];
  X(:,m) = fft(xwzp);
  % advance input offset by hop size
  xoff = xoff + nhop; 
end

if (nargout==0) || doplot
  t = (0:nframes-1)*nhop/fs;
  f = 0.001*(0:nfft-1)*fs/nfft;
  Xdb = 20*log10(abs(X));
  Xmax = max(max(Xdb));
  % Clip lower limit to -dbdown dB so nulls don't dominate:
  clipvals = [Xmax-dbdown,Xmax];
  imagesc(t,f,Xdb,clipvals);
  % grid;
  axis('xy');
  colormap(jet);
  xlabel('Time (sec)');
  ylabel('Freq (kHz)');
end