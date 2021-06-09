function [G1,G2] = ConstrIM(freqsReal,dimY)
% G1 = ConstrIMReal(freqs,dimY)
%
% Construct a real block-diagonal internal model
% freqs = Frequencies to be included in the controller, only real nonnegative
% frequencies, if zero frequency is included, it's the first element in the
% vector
% dimY = Number of copies of each frequency to be included in G1
% G1 = System matrix of the internal model


if freqsReal(1) == 0
    G1 = blkdiag(zeros(dimY),kron(diag(freqsReal(2:end)),[zeros(dimY),eye(dimY);-eye(dimY),zeros(dimY)]));
    G2 = [eye(dimY);kron(ones(length(freqsReal)-1,1),[eye(dimY);zeros(dimY)])];
else
    G1 = blkdiag(kron(diag(freqsReal),[zeros(dimY),eye(dimY);-eye(dimY),zeros(dimY)]));
    G2 = kron(ones(length(freqsReal),1),[eye(dimY);zeros(dimY)]);
end


% OLD IMPLEMENTATION:
%
% q = length(freqsReal);
% if freqsReal(1)==0, dimZ = dimY*(2*q-1); else, dimZ = dimY*2*q; end  
%   
% G1 = zeros(dimZ);
% 
% if freqsReal(1)==0
%   zoffset = dimY; 
%   nzfreqs = freqsReal(2:end);
%   
%   G2 = [eye(dimY);repmat([eye(dimY);zeros(dimY)],q-1,1)];
% 
% else
%   zoffset = 0;
%   nzfreqs = freqsReal;
%   
%   G2 = repmat([eye(dimY);zeros(dimY)],q,1);
% 
% end
% 
% for ind = 1:length(nzfreqs)
%   indran = zoffset+(ind-1)*2*dimY+(1:(2*dimY));
% 
%   G1(indran,indran) = nzfreqs(ind)*[zeros(dimY) eye(dimY);-eye(dimY) zeros(dimY)];
%   
% end
% 
