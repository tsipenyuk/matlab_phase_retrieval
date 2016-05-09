function [K, RHS, uh_k, uh_kmq, uh_kpq, quot_k, phase_k] = wfpr_1d(k, q, k_grid, u_hat, sqrtI, h)
% WFK_1D Wasserstein Fredholm Phase Retrieval --- kernel and RHS
%
%   [K, RHS, ~, ~, ~, ~, ~] = WFK_1D(k, q, k_grid, u_hat, sqrtI, h).
%   Returns the kernel and the right-hand side of the Fredholm
%   equation corresponding to implicit Euler applied to the phase
%   retrieval modulus energy in the quadratic Wasserstein space. This
%   function is constructed to be used with the Fie library, see also FIE.
%
%   k = an array of the type that is returned by NDGRID
%   function. Describes a coordinate in an N-dimensional space.
%    
%   q = an array of the type that is returned by NDGRID
%   function. Describes a coordinate over which the Fredholm integral
%   is taken.
%
%   k_grid = a one-dimensional vector. Describes the
%   coordinates at which functions u_hat and sqrtI are
%   evaluated. (Mathematically obsolete, since u_hat and sqrtI are
%   fftshifted on the Fourier side anyway. Kept explicit to allow for
%   numerical tweaking.)  
%
%   u_hat = Fourier transform of the current approximation. MUST
%   vanish at the boundaries to ensure correct interpolation
%   behaviour.  
%
%   sqrtI = Fourier transform modulus of the solution (square root of
%   the measured intensity). MUST varish at the boundaries to ensure
%   correct interpolation behaviour.  
%
%   h = scalar number, step length of the gradient step.  
%
%   K = an array of the same size as k and q, returns the values of
%   the Fredholm kernel.  
%
%   RHS = an array of the same size as k, returns the values of
%   the Fredholm equation right-hand side.  
%
%   All other output variables serve the debugging purposes and shall
%   be removed on release.
%
%   Arseniy Tsipenyuk, TUM M7
%   May 09, 2016

    % The term $\sqrt{I} / |\hat u|$ is divergent for small values of
    % $|\hat u|$. To ensure numerical stability, perform a cut-off for
    % large values of the quotient
    TOL = 1e6;
    one_div_u = ones(size(u_hat)) ./ abs(u_hat);
    one_div_u(one_div_u > TOL) = TOL;
    quot = sqrtI .* one_div_u;
    
    % Linear interpolation of the density FT: The interpolation must be
    % two-dimensional to ensure the compability with the Fie library,
    % so we blow up the results in a dummy dimension
    dummy = 1:1:2;
    [KG, dG] = ndgrid(k_grid, dummy);
    [U_int, ~] = ndgrid(u_hat, dummy);
    [I_int, ~] = ndgrid(sqrtI, dummy);
    [Q_int, ~] = ndgrid(quot, dummy);

    % Interpolation
    u_hat_func = griddedInterpolant(KG, dG, U_int, 'linear', 'none');
    sqrtI_func = griddedInterpolant(KG, dG, I_int, 'linear', 'none'); 
    quot_func  = griddedInterpolant(KG, dG, Q_int, 'linear', 'none');

    % Transpose the variables to ensure optimal interpolant evaluation
    % efficiency in Fie (which uses meshgrids, not ndgrids)
    uh_k = u_hat_func(k', q');
    % It may be bad to calculate these quantities after interpolation.
    uh_kmq = u_hat_func(k' - q', q');
    uh_kpq = u_hat_func(k' + q', q');
    sqrtI_k = sqrtI_func(k'); % Unnecessary evaluation may be optimized
    quot_k = quot_func(k', q');
    
    % Set values beyond the interpolation grid to 0
    uh_k(isnan(uh_k)) = 0;
    uh_kmq(isnan(uh_kmq)) = 0;
    uh_kpq(isnan(uh_kpq)) = 0;
    sqrtI_k(isnan(sqrtI_k)) = 0;
    quot_k(isnan(quot_k)) = 0;
    
    phase_k = exp(1j * angle(uh_k));

    K_part = 1j .* q' .* (uh_kmq .* (1 - quot_k) + ...
             phase_k .* quot_k .* ...
               (conj(phase_k) .* uh_kmq + phase_k .* conj(uh_kpq)));
    K = -1j .* h .* k' .* K_part;
    RHS = -1j .* k .* (uh_k - sqrtI_k .* exp(1j .* angle(uh_k)));
    
    % Transpose the result back to ensure that K has the same size as k
    % and q.
    K = K'; 
    RHS = RHS';
end