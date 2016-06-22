function res = pM(g, A)
% pM Replace Fourier modulus of g by A
%   Calculates the Fourier transform of g, replaces its
%   modulus by A, performs the backward Fourier transform.
    G = fftn(g);
    % Replace NaN-values in the measurement by current G-values
    A_mod = A;
    A_mod(isnan(A_mod)) = abs(G(isnan(A_mod))); 
    res = ifftn(A_mod .* exp(1i * angle(G)));
end