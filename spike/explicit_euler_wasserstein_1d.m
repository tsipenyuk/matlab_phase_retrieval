function [g_new, error, w_hat_s, w_com] = explicit_euler_wasserstein_1d(g, sqrtI, update_params)
% ER Wasserstein gradient flow algorithm (WGF-alg)
% 
% Perform one step of the WGF-alg, using an array g as
% the input density and an array A (of the same size as g) as the
% square root of the measured intensity.
% Calculate the modulus energy 
% $$
%          error = \int (|\hat g(k)| - sqrtI(k))^2 dk.
% $$
% Use the explicit Euler update.
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
    
    g_func = griddedInterpolant(update_params.x_list, g, ...
                                'linear', 'none');
    pM_g_func = griddedInterpolant(update_params.x_list, pM(g, sqrtI), ...
                                   'linear', 'none');
    
    g_com = g_func(x_com);
    pM_g_com = real(pM_g_func(x_com));
    
    w_com = - gradient(g_com - pM_g_com);
    
    x_com_new = x_com + update_params.h * w_com;
    x_rpt_new = rptCoordinates(update_params.x_list(1), x_com_new);
    
    % Calculate the new density
    x_tmp = sort([update_params.x_list(1) x_rpt_new]);
    m = F_pcs(1);
    g_rpt_new = zeros(size(x_rpt_new));
    for i = 1:1:length(g_rpt_new)
        g_rpt_new(i) = m / (x_tmp(i + 1) - x_tmp(i));
    end
    
    % Leftmostpoint is a dummy to ensure non-NaN values for the
    % leftmost part of the interpolation
    g_func_new = griddedInterpolant(x_rpt_new, ...
                                    g_rpt_new, ...
                                    'linear', 'none');
    g_new = g_func_new(update_params.x_list);
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
%        figure;
%        hold on;
%        plot(x_com, F_pcs, 'xb'); 
%        plot(x_com_new, F_pcs, 'og');
%        plot(x_rpt_new, F_pcs, '-r');
%        yyaxis right
%        plot(x_rpt_new, g_rpt_new, '.b');
%        plot(update_params.x_list, g_new, '-g');
%        hold off;
%%    % Debugging info -- end
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
