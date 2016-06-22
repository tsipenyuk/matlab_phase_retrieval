function res = eM(g, A)
% EM Fourier modulus energy (error) functional 
%   Estimates how far is the Fourier modulus $|\hat g|$
%   from the desired Fourier modulus $A$.

    G_modulus = abs(fftn(g));

    % Correct for the NaN-values in the data
    A_corr = A;
    A_corr(isnan(A)) = G_modulus(isnan(A));
    
    B = (G_modulus - A_corr).^2; 
    res = sum(B(:));
end