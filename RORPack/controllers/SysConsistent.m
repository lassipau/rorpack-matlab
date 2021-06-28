function Sys = SysConsistent(Sys,yref,wdist,freqsReal)
% function Sys = SysConsistent(Sys,yref,wdist,freqsReal)
%
% Checks the consistency of the RORPack system 'Sys' and the signals
% 'yref', 'wdist', and the list of (real) frequencies 'freqsReal'. If the
% disturbance input matrix Bd is missing, this is replaced with a zero
% matrix.

dimX = size(Sys.A,1);
dimU = size(Sys.B,2);
dimY = size(Sys.C,1);

if ~isfield(Sys,'D')
    Sys.D = zeros(dimY,dimU);
    fprintf('Feedthrough matrix D not given, defining as D=0.\n')
end

% Check that the dimensions of A, B, C, and D are consistent
checkA = isequal(size(Sys.A,2),dimX);
checkBC = isequal(size(Sys.B,1),dimX) && isequal(size(Sys.C,2),dimX);
checkD = isequal(size(Sys.D),[dimY,dimU]);

if ~checkA || ~checkBC || ~checkD
    error('Dimensions of A, B, C, and D are not consistent!')
end

% If the disturbance input matrix is not defined, define it as zero
if ~isfield(Sys,'Bd')
    if nargin >= 3 
        Sys.Bd = zeros(dimX,length(wdist(0)));
        fprintf('Disturbance input matrix Bd not given, defining as Bd=0.\n')
    else
        Sys.Bd = zeros(dimX,1);
        fprintf('Disturbance input matrix Bd not given, defining as Bd=0 for SCALAR disturbance.\n')
    end
end

% Check the possible matrices Cm and Dm (separate measured output of the
% system)
if isfield(Sys,'Cm')
    
    if ~isequal(size(Sys.Cm,2),dimX)
        error('Dimensions of Cm and A are not consistent!')
    end
    if isfield(Sys,'Dm')
        if ~isequal(size(Sys.Dm),[size(Sys.Cm,1),dimU])
            error('Dimensions of Cm and Dm are not consistent!')
        end
    else
        Sys.Dm = zeros(size(Sys.Cm,1),dimU);
        fprintf('Measurement feedthrough matrix Dm not given, defining as Dm=0.\n')
    end
    
end
        

% Check that the dimension of 'yref' is consistent with C.
if nargin >= 2 && ~isequal(size(yref(0)),[dimY,1])
    error('The dimensions of yref and C are not consistent!')
end

% Check that the dimension of 'wdist' is consistent with Bd.
if nargin >= 3 && ~isequal(size(wdist(0)),[size(Sys.Bd,2),1])
    error('The dimensions of wdist and Bd are not consistent!')
end

% Check that the frequencies are real and nonegative, do not contain
% repeated elements, and that zero (if present) is the first element of 
% the vector 'freqsReal'.
if nargin >= 4 
    
    if ~isreal(freqsReal)
        error('The frequencies in freqsReal are not real!')
    end
    if ~isequal(length(unique(freqsReal)),length(freqsReal))
        error('The vector freqsReal contains repeated frequencies!')
    end
    
    zero_ind = find(freqsReal==0);
    if ~isempty(zero_ind) && ~isequal(zero_ind,1)
        error('The zero frequency is not the first element of freqsReal!')
    end
    
    if any(freqsReal<0)
        error('The frequencies in freqsReal are not nonnegative!')
    end
end

fprintf('Consistency check passed.\n')