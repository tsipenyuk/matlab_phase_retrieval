function [g_new, error] = wasserstein_1d(g, A, update_params)
% ER Wasserstein gradient flow algorithm (WGF-alg)
% 
% Perform one step of the WGF-alg, using an array g as
% the input density and an array A (of the same size as g) as the
% square root of the measured intensity.
% Calculate the modulus energy 
% $$
%          error = \int (|\hat g(k)| - A(k))^2 dk.
% $$
%
% INPUT
%
% g      one-dimensional array describing the density
% A      one-dimensional array describing the desired Fourier
%        modulus
%
% update_params    structure array containing following tuning
%                  details of the algorithm
%
% update_params.num_mass_pcs - number of mass pieces into which the
%                              density is split. These pieces are
%                              then moved according to the
%                              Wasserstein gradient. 
% update_params.step_size    - step size of the gradient descent.
%
% ALGORITHM
% The algorithm performs the following steps.
% 1. Split the mass of the density g into
% update_params.num_mass_pcs points. Store the coordinates of these
% points in the vector x.
% 2. Calculate the Wasserstein gradient at x.
% 3. Move points x by  (values stored in the gradient * step_size).
% 4. Use linear interpolation to calculate the new density g_new.


    g_new = pP(pM(g, A));
    error = eM(g_new, A);
end