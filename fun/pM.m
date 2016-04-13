function res = pM(g, A)
% PM Replace Fourier modulus of g by A
%   Calculates the Fourier transform of g, replaces its
%   modulus by A, performs the backward Fourier transform.
    G = fftn(g);
    res = ifftn(A .* exp(1i * angle(G)));
end