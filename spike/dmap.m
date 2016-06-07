function [g_new, error] = dmap(g, A, beta)
% dmap - Difference map algorithm (optimized variant)
%
% Synopsis ([]s are optional)
%   [g_new, error] = hio(g, A, [beta])
%
% Description
%   Performs a phase retrieval update step and calculates the error
%   (energy) corresponding to optimized difference map algorithm 
%   (optimized version) as described in [Elser], p. 45.
%
% Inputs ([]s are optional)
%   (ndarray) g      Current  approximation to the solution of the
%                    phase retrieval problem
%   (ndarray) A      Phase retrieval data (square root of the
%                    measured intensity)
%   (scalar)  [beta = -1]
%                    Update parameter, one of the values used in
%                    in [Elser]. Setting beta=1 yields bauschke_hio.
%
% Outputs ([]s are optional)
%   (ndarray) g_new  Updated approximation to the solution of the 
%                    phase retrieval problem
%   (scalar)  error  Error (energy) corresponding to the updated 
%                    approximation
%
% Examples
%   %% 1D, two gaussians
%   x = [-20:0.2:20];
%   g_sol = exp(-x.^2);
%   shift = fix(length(g)/4);
%   g_sol = circshift(g, [0, shift]) + circshift(g, [0, -shift]);
%   A = abs(fftn(g_sol));
%   g_new = A;
%   E = [];
%   for i=1:1:200
%       [g_new, error] = dmap(g_new, A);
%       E = [E error];
%   end
%   plot(E);
%   plot(g_new);
%
%   %% 2D, two gaussians
%   [x1, x2] = ndgrid([-20:0.2:20], [-20:0.2:20]);
%   g = exp(-x1.^2 - x2.^2);
%   shift = fix(length(g)/4);
%   g_sol = circshift(g, [0, shift]) + circshift(g, [0, -shift]);
%   A = abs(fftn(g_sol));
%   g_new = A;
%   E = [];
%   for i=1:1:200
%       [g_new, error] = dmap(g_new, A);
%       E = [E error];
%   end
%
% See also
%   er
%   bio
%   hio_bauschke
%   hio_fienup
%   
% Requirements
%   pM (modulus projection)
%   pP (non-negative projection)
%
% References
%   V. Elser, “Phase retrieval by iterated projections,”
%       J. Opt. Soc. Am. A., vol. 20, pp. 40–55, 2003.
%   doc/phase_retrieval_algorithms.pdf
%
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase Retrieval Sandbox root folder.
%
% Changes
%   2016-06-01  First Edition
    if nargin == 2
        beta = -1; % you may try other values, [Elser] 
                   % uses 1 (HIO with hio.beta=1), 0.5 
                   % and 0.3 among others.
    end

    pM_g = pM(g, A);
    pP_g = pP(g);
    % Double check! Suspicious! Verify beta=1 yielding HIO!
    g_new = g - pP_g - pM_g - (1 + beta) * pP(pM_g) + ...
                              (1 - beta) * pM(pP_g, A);
    error = eM(g_new, A);
end