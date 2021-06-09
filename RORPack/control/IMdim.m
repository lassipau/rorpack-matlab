function dimZ0 = IMdim(freqs,p)
% function dimZ0 = IMdim(freqs,p)
%
% Compute the dimension of the internal model. 
% Parameters: 
% 'freqs' =  the (real) frequencies of the internal model. If 0 is a 
% frequency, it is required to be the first element of the vector.
% 'p' = number of outputs of the system
%
% 'dimZ0' = dimension of the internal model.

if freqs(1) == 0
    dimZ0 = (2*length(freqs)-1)*p;
else
    dimZ0 = 2*length(freqs)*p;
end