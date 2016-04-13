function res = eM(g, A)
% EM Fourier modulus error functional 
%   Estimates how far is the Fourier modulus $|\hat g|$
%   from the desired Fourier modulus $A$.
    res = sum((abs(fftn(g)) - A).^2)
end