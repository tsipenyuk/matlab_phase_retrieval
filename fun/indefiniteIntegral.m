function fout = indefiniteIntegral(xin, fin)
% INDEFINITEINTEGRAL one-dimensional trapezoidal integral
%
%   fout = indefiniteIntegral(xin, fin) returns the indefinite
%   integral corresponding to the discretized function at xin with
%   values fin. The first entry of the integral is set to 0.
%
%   xin = vector containing points at which the function is
%   discretized
%
%   fin = vector containing function values at xin
%
%   fout = indefinite integral of fin discretized at xin
%
%   Arseniy Tsipenyuk, TUM M7
%   May 09, 2016
    fout = zeros(size(fin));
    for i = 1 : 1 : length(fin) - 1
        fout(i+1) = fout(i) +  (xin(i+1) - xin(i)) * (fin(i+1) + fin(i)) / 2;
    end
end

