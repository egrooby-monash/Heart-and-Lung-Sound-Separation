function [heartSound,lungSound]=singular_spectrum_analysis(mixture,fs1,l)
%% Paper Information
% Localizing Heart Sounds in Respiratory Signals Using Singular Spectrum Analysis
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5959957 
% Code based on online example https://www.mathworks.com/matlabcentral/fileexchange/58967-singular-spectrum-analysis-beginners-guide
%% Purpose
% localise heart sounds
%% Inputs
% X= mixutre
% fs1= sampling frequency
% l= window length
%% Outputs
% Separated heart sound (HS) and lung sound (LS) 
if nargin<3
    %window length
    l=200;
end

%% 1 Decomposition
%% 1a Embedding
% signal
s=mixture; 
% Resample to 4000Hz
s = resample(s,4000,fs1); 
fs=4000; 
% high pass filter 50Hz
[b,a] = butter(2,50/(fs/2),'high');
s = filtfilt(b,a,s);

% signal length
r=length(s); 
k=r-l+1;
% X=zeros(l,k);
% for m=1:l
%   X(m,:) = s((1:k)+m-1);
% end

X=zeros(k,l);
for m=1:l
  X(:,m) = s((1:k)+m-1);
end

%% 1b Singular value decomposition
% S=X*X';
S=X'*X; 
% S=X'*X/k; 

% get eigenvalues (lambda) and eigvectors (lambda)
[RHO,LAMBDA] = eig(S);
% extract the diagonal elements
LAMBDA = diag(LAMBDA);   
% sort eigenvalues
[LAMBDA,ind]=sort(LAMBDA,'descend'); 
% sort eigenvectors
RHO = RHO(:,ind);                    

% temporal principal components
%PC = X'*RHO;
PC=X*RHO; 

%% 2 Reconstruction
%% 2a Grouping
% removing noise and remaining is heart and lung

%FFT by column
Y=fft(RHO,[],1);
P2 = abs(Y/l);
P_RHO = P2(1:l/2+1,:);
f_RHO = (fs/l)*(0:(l/2))/l;

Y=fft(PC,[],1);
FFT_SIZE=size(Y,1); 
P2 = abs(Y/FFT_SIZE);
P_PC = P2(1:round(FFT_SIZE/2)+1,:);
f_PC = fs*(0:round(FFT_SIZE/2))/FFT_SIZE;

stop=find(cumsum(LAMBDA)/sum(LAMBDA)>0.95,1);
pairs=[];
freq_pairs_RHO=[]; 
freq_pairs_PC=[]; 
for i=1:stop-1
    for j=i+1:stop-1
        if abs(1-LAMBDA(i)/LAMBDA(j))<0.06 || abs(1-LAMBDA(j)/LAMBDA(i))<0.05
            [~, f1]=max(P_RHO(:,i));
            [~, f2]=max(P_RHO(:,j));
            freq1_RHO=f_RHO(f1);
            freq2_RHO=f_RHO(f2);
            
            [~, f1]=max(P_PC(:,i));
            [~, f2]=max(P_PC(:,j));
            freq1_PC=f_PC(f1);
            freq2_PC=f_PC(f2);
            if abs(1-freq1_RHO/freq2_RHO)<0.03 || abs(1-freq2_RHO/freq1_RHO)<0.03
                pairs=[pairs;i j];
                freq_pairs_RHO=[freq_pairs_RHO;freq1_RHO freq2_RHO];
                freq_pairs_PC=[freq_pairs_PC;freq1_PC freq2_PC];
            end
        else
            break
        end
    end
end

% % Grouping heart
% %need some rule on eigensum to select optimal pairs to differentiate
% %between heart and lung as in the paper is based on viewing a graph which
% %is not appropriate. But they concluded first 3 pairs is heart.
% if size(pairs,1)>=3
%     HS_loc=unique(pairs(1:3,:)); 
% else
%     HS_loc=unique(pairs); 
% end 
% 

% top 3 pairs and frequency less than 250Hz peak
acceptable_pairs=mean(freq_pairs_PC(1:min(3,size(pairs,1)),:),2)<250; 
HS_loc=unique(pairs(acceptable_pairs,:));


% Grouping lung
% lung locations are all the ones not heart
LS_loc=1:stop-1; 
for i=1:length(HS_loc)
    LS_loc(LS_loc==HS_loc(i))=[];
end
%% 2b Diagonal averging 
% invert projection
buf=PC(:,HS_loc)*RHO(:,HS_loc)'; 
buf=buf(end:-1:1,:);
HS=zeros(1,r);
% anti-diagonal averaging
for n=1:r 
    HS(n)=mean( diag(buf,-k+n) );
end

% invert projection
buf=PC(:,LS_loc)*RHO(:,LS_loc)'; 
buf=buf(end:-1:1,:);
LS=zeros(1,r);
% anti-diagonal averaging
for n=1:r 
    LS(n)=mean( diag(buf,-k+n) );
end

heartSound=zeros(r,1); 
lungSound=zeros(r,1); 
heartSound(end-length(HS)+1:end)=HS; 
lungSound(end-length(LS)+1:end)=LS; 

end