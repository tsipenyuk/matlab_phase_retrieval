function RHS = wfrhs_1d(k, k_grid, u_hat, sqrtI)
% WFK_1D Wasserstein Fredholm equation right-hand side
%
%  Fredholm right-hand side of the  Wasserstein gradient
%  corresponding to the phase retrieval modulus energy.
%  
%  Warning. Bad coding style: partly repeats the code from wfpr_1d.
%  Think about possible redesign for future releases.
%
%  Arseniy Tsipenyuk, TUM M7
%  May 10th, 2016

    % Linear interpolation of the density FT:
    % External interpolation 'none' defaults to NaN values. We want
    % to use external interpolation to 0 instead.
    u_hat_func = griddedInterpolant(k_grid, u_hat, 'linear', 'none');
    sqrtI_func = griddedInterpolant(k_grid, sqrtI, 'linear', 'none');

    uh_k = u_hat_func(k);
    sqrtI_k = sqrtI_func(k);
    
    uh_k(isnan(uh_k)) = 0;
    sqrtI_k(isnan(sqrtI_k)) = 0;

    RHS = -1j .* k .* (uh_k - sqrtI_k .* exp(1j .* angle(uh_k)));
end
