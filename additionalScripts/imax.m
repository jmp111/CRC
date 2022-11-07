%%  imax
% ------------------------------------------------------------------------
%   Input:
%   x: vector or matrix
%   
%   Output:
%   im: the index (vector) or indices (matrix) of the maxima
%   
% ------------------------------------------------------------------------
% J.M. Posma
% 25 February 2010
%%
function im=imax(x)
    [d,im]=max(x);
end