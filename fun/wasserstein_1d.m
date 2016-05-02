function [g_new, error, w_hat_s, w_com] = wasserstein_1d(g, sqrtI, update_params)
% ER Wasserstein gradient flow algorithm (WGF-alg)
% 
% Perform one step of the WGF-alg, using an array g as
% the input density and an array A (of the same size as g) as the
% square root of the measured intensity.
% Calculate the modulus energy 
% $$
%          error = \int (|\hat g(k)| - sqrtI(k))^2 dk.
% $$
%
% INPUT
%
% g      one-dimensional array describing the density
% sqrtI  one-dimensional array describing the desired Fourier
%        modulus
%
% update_params    structure array containing following tuning
%                  details of the algorithm
%
% update_params.x_list       - spacial coordinates of the density g
% update_params.k_list       - Fourier coordinates of \hat g
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
    
    % Right ends of mass pieces
    [x_pcs, F_pcs] = splitMass(update_params.num_mass_pcs, ...
                               update_params.x_list, g);
    
    M = sum(g);
    disp(['Total mass: ' num2str(M)]);

    % Center-of-mass coordinates of mass pieces
    x_com = comCoordinates(update_params.x_list(1), x_pcs);
    
    % Calculate the Wasserstein gradient at x_com
    % Use Fie library to solve the corresponding Fredholm equation
    g_hat = fft(g);
    function res = K(k_, q_)
        res = wfk_1d(k_, q_, update_params.k_list, ...
                     g_hat, sqrtI, update_params.h);
    end
    function res = RHS(k_)
        res = wfrhs_1d(k_, update_params.k_list, g_hat, sqrtI);
    end

    % Call the Fie library
    [sol, errest, cond] = Fie(1, update_params.k_list(1), ...
                              update_params.k_list(end), ...
                              1, @K, @RHS, 1e-6, 1e-3); 

    % Interpolate the solution at the Fourier space points
    w_hat_s = ntrpFie(sol, linspace(1, length(g_hat), ...
                                    length(g_hat)));
    % Due to internal interpolation properties, see wfk_1d and
    % wfrhs_1d, we need to back-shift the result
    w_hat = ifftshift(w_hat_s);
    w_list = real(ifft(w_hat));
    
    w_func = griddedInterpolant(update_params.x_list, w_list, 'linear', ...
                                'none');
    % Calculate the velocity at the center-of-mass points
    w_com = w_func(x_com);
    x_com_new = x_com + update_params.h * w_com;
    x_rpt_new = rptCoordinates(update_params.x_list(1), x_com_new);
    
    % Calculate the new density
    x_tmp = [update_params.x_list(1) x_rpt_new];
    m = F_pcs(1);
    g_rpt_new = zeros(size(x_rpt_new));
    for i = 1:1:length(g_rpt_new)
        g_rpt_new(i) = m / (x_tmp(i + 1) - x_tmp(i));
    end
    
    % Leftmostpoint is a dummy to ensure non-NaN values for the
    % leftmost part of the interpolation
    [sort_x_rpt_new, sort_index] = sort(x_rpt_new);
    sort_g_rpt_new = g_rpt_new(sort_index);
    g_func_new = griddedInterpolant(sort_x_rpt_new, ...
                                    sort_g_rpt_new, ...
                                    'linear', 'none');
    g_new = g_func_new(update_params.x_list);

    % If the values of g_new lie beyond the boundaries of
    % interpolation, set them in a way that ensures the 
    % overall mass conservation
    M_new = sum(g_new(~isnan(g_new)));
    if M >= M_new
        len = sum(isnan(g_new));
        small_mass_fill = (M - M_new) / len;
        g_new(isnan(g_new)) = small_mass_fill; 
    else
        disp(['Something''s wrong - new mass is too big...' num2str(M) ...
              ' vs. new mass ' num2str(M_new)]);
    end

    error = eM(g_new, sqrtI);
    
    % Debugging info 
        figure;
        hold on;
        subplot(2,1,1)
        plot(x_com, F_pcs, 'xb'); 
        plot(x_com_new, F_pcs, 'og');
        plot(x_rpt_new, F_pcs, '-r');
        yyaxis right
        plot(x_rpt_new, g_rpt_new, '-b');
        plot(update_params.x_list, g_new, '.g');
        subplot(2,1,2);
        plot(real(w_hat));
        plot(imag(w_hat));
        hold off;
        % Debugging info -- end
end


%==== Function used to calculate the cumulative distribution function ====
function fout = myCdf(xin, fin)
    fout = fin;
    for i = 1 : 1 : length(fin) - 1
        fout(i+1) = fout(i) + fout(i+1) * (xin(i+1) - xin(i));
    end
    fout = fout - fout(1);
end % myCdf
%==========================================================================
