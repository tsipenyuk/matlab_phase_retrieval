function RHS = wfrhs_1d(k, k_grid, u_hat, sqrtI)
% WFK_1D Wasserstein Fredholm equation right-hand side
%
% Fredholm right-hand side of the  Wasserstein gradient
% corresponding to the phase retrieval modulus energy.
%
% Warning. Bad coding style: partly repeats the code from wfk_1d.
% Think about possible redesign for future releases.

    % For the interpolation of u_hat and sqrtI to be correct, they
    % must vanish at the beginning and the end:
    u_hat_0 = fftshift(u_hat);
    sqrtI_0 = fftshift(sqrtI);

    % Linear interpolation of the density FT:
    % External interpolation 'none' defaults to NaN values. We want
    % to use external interpolation to 0 instead.
    u_hat_func = griddedInterpolant(k_grid, u_hat_0, 'linear', 'linear');
    sqrtI_func = griddedInterpolant(k_grid, sqrtI_0, 'linear', 'linear'); 

    uh_k = u_hat_func(k);
    sqrtI_k = sqrtI_func(k);
    RHS = -1j .* k .* (uh_k - sqrtI_k .* exp(1j .* angle(uh_k)));
end
