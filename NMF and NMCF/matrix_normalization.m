function [W]=matrix_normalization(W)
%% Paper Information
% Sparse NMF--half-baked or well done?
% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.698.2772&rep=rep1&type=pdf 
%% Purpose
% Normalise W 
%% Inputs
% Dictionary matrix W
%% Outputs
% Normalised matrix W

F=size(W,1);
Wn=ones(F)*W.^2;

W=W./Wn.^0.5;
end
