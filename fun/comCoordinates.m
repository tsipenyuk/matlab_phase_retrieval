function x_com = comCoordinates(x_pcs, f_pcs)
% COMCOORDINATES Center-of-mass coordinates 
%
%   x_com = COMCOORDINATES(x_pcs, f_pcs) Given a piecewise linear
%   function, calculate the center-of-mass coordinates of the
%   pieces.
%
%   x_pcs = one-dimensional arrray of length n, n being any
%   integer, coordinates, at which the function is evaluated
%
%   f_pcs = one-dimensional array of length n, approximative values
%   of the function
%
%   x_com = one-dimensional array of length n-1, center-of-mass
%   coordinates of the mass pieces.
%
%   The calculation is performed using the following formula: for a
%   trapezoid defined by values f1 and f2 at the coordinates x1 and
%   x2,
%   $$
%   x_c = (x2^2 f2 - x1^2 f1) / 3 + (x1 + x2)(f1 x2 - f2 x1)/6.
%   $$
    
    x1 = x_pcs(1:end-1);
    x2 = x_pcs(2:end);
    f1 = f_pcs(1:end-1);
    f2 = f_pcs(2:end);
    
    x_com = (x2.^2 .* f2 - x1.^2 .* f1) / 3 + ...
            (x1 + x2) .* (f1 .* x2 - f2 .* x1) / 6;
end
