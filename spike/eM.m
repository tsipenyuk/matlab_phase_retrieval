function res = eM(g, A)
% EM Fourier modulus energy (error) functional 
%   Estimates how far is the Fourier modulus $|\hat g|$
%   from the desired Fourier modulus $A$.
    B = (abs(fftn(g)) - A).^2; 
    res = sum(B(:));
end