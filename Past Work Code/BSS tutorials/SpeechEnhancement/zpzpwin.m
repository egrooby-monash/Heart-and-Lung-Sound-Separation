function [output_sig] = zpzpwin(input_sig, window, padded_length)
% 420 hw2q1b
% function [output_sig] = zpzpwin(input_sig, window, padded_length)
% zero phase windows a signal and adds zeros to the end
% check that input_sig is the same size vector as the window:
if sum(size(input_sig) ~= size(window)) ...
| (length(input_sig)+1)/2 ~= round((length(input_sig)+1)/2) ... check odd
| min(size(window))~= 1 ...
| length(size(window)) ~= 2,
error('input_sig and window must be odd length matching-dimension vectors');
end % NOTE: code checking length (instead of size)
% also OK, but then must ALSO
% transpose in step below if mismatch
% OR: must refuse to allow either column or row vectors
% if not checking
% multiply window by input_sig:
temp = input_sig.*window;
% add zeros to end so that length(output_sig) is padded_length:
% add zeros in column or row as appropriate:
if size(input_sig,1) > size(input_sig,2) % column
output_sig = [temp((length(input_sig)+1)/2:end); ...
zeros(padded_length-length(input_sig),1); ...
temp(1:(length(input_sig)-1)/2)];
else
output_sig = [temp((length(input_sig)+1)/2:end), ...
zeros(1,padded_length-length(input_sig)), ...
temp(1:(length(input_sig)-1)/2)];
end
% now important last step: shift zeros to center