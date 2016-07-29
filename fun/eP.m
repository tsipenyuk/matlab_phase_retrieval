function res = eP(g)
% EP Non-negativity energy (error) functional 
%   Estimates how large is the negative part of the solution
    res =  g;
    res(g>=0) = 0;
    res = res.^2;
    res = sum(res(:));
end