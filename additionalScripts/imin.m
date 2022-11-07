%%  imin
% ------------------------------------------------------------------------
%   Input:
%   x: vector or matrix
%   
%   Output:
%   im: the index (vector) or indices (matrix) of the minima
%   
% ------------------------------------------------------------------------
% J.M. Posma
% 25 February 2010
%%
function im=imin(x)
    [d,im]=min(x);
end