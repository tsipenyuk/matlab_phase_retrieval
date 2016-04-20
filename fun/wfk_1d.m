function [K, uh_k, uh_kmq, uh_kpq, quot_k, phase_k] = wfk_1d(k, q, k_grid, u_hat, sqrtI, h)
% WFK_1D Wasserstein Fredholm Kernel
%
% Fredholm Kernel of the  Wasserstein gradient corresponding to the
% phase retrieval modulus energy.
 
    % For the interpolation of u_hat and sqrtI to be correct, they
    % must vanish at the beginning and the end:
    u_hat_0 = fftshift(u_hat);
    sqrtI_0 = fftshift(sqrtI);

    % The term $\sqrt{I} / |\hat u|$ is divergent for small values of
    % $|\hat u|$. To ensure numerical stability, perform a cut-off
    % for large values of the quotient
    TOL = 1000;
    quot = sqrtI_0 ./ abs(u_hat_0);
    quot(quot > TOL) = TOL;
    quot(quot < -TOL) = -TOL;
    
    % Linear interpolation of the density FT:
    % The interpolation must be two-dimensional to ensure the
    % compability with Fie library, so we blow up the results in a
    % dummy dimension
    dummy = 1:1:2;
    [KG, dG] = ndgrid(k_grid, dummy);
    [U_int, ~] = ndgrid(u_hat_0, dummy);
    [I_int, ~] = ndgrid(sqrtI_0, dummy);
    [Q_int, ~] = ndgrid(quot, dummy);

    % Interpolation
    u_hat_func = griddedInterpolant(KG, dG, U_int, 'linear', 'linear');
    sqrtI_func = griddedInterpolant(KG, dG, I_int, 'linear', 'linear'); 
    quot_func  = griddedInterpolant(KG, dG, Q_int, 'linear', 'linear');

    % Dependence on the 2nd argument is the dummy dependence here
    % Transpose to ensure optimal interpolant evaluation efficiency 
    % in Fie (which uses meshgrids, not ndgrids)
    uh_k = u_hat_func(k', q');
    % It may be bad to calculate these quantities after
    % interpolation.
    uh_kmq = u_hat_func(k' - q', q');
    uh_kpq = u_hat_func(k' + q', q');

    
    quot_k = quot_func(k', q');
    phase_k = exp(1j * angle(uh_k));

    K_part = 1j .* q' .* (uh_kmq .* (1 - quot_k) + ...
         phase_k .* quot_k .* (conj(phase_k) .* uh_kmq + ...
                               phase_k .* conj(uh_kpq)));
    K = - 1j .* h .* k' .* K_part;
    K = K';
end