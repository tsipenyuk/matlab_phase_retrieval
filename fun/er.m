function [g_new, error] = er(g, A)
% ER Error reduction algorithm
% 
% Perform an error-reduction algorithm step, using an array g as
% the input density and an array A (of the same size as g) as the
% square root of the measured intensity, so that 
% $$
%    |\hat g_{sol}| = A.
% $$
% for a phase problem solution $g_{sol}$. Calculate the modulus
% energy 
% $$
%          error = \int (|\hat g(k)| - A(k))^2 dk.
% $$
    g_new = pP(pM(g, A));
    error = eM(g_new, A);
end