function [g_new, error, w_hat_s, w_com] = wasserstein_1d(g, sqrtI, update_params)
% ER Wasserstein gradient flow algorithm (WGF-alg)
% 
%   Perform one step of the WGF-alg, using an array g as
%   the input density and an array A (of the same size as g) as the
%   square root of the measured intensity.
%   Calculate the modulus energy 
%   $$
%            E[g] = \int (|\hat g(k)| - sqrtI(k))^2 dk.
%   $$
%   
%   INPUT
%   
%   g     = one-dimensional array describing the density
%
%   sqrtI = one-dimensional array describing the desired Fourier
%   modulus
%   
%   update_params =  structure array containing following tuning
%                    details of the algorithm
%   
%   update_params.x_list       = spacial coordinates of the density g
%   update_params.k_list       = Fourier coordinates of \hat g
%   update_params.num_mass_pcs = number of mass pieces into which the
%                                density is split. These pieces are
%                                then moved according to the
%                                Wasserstein gradient. 
%   update_params.step_size    = step size of the gradient descent.
%   
%   ALGORITHM
%   The algorithm performs the following steps.
%   1. Split the mass of the density g into
%   update_params.num_mass_pcs points. Store the coordinates of these
%   points in the vector x.
%   2. Calculate the Wasserstein gradient at x.
%   3. Move points x by  (values stored in the gradient * step_size).
%   4. Use linear interpolation to calculate the new density g_new.
    
    % Right ends of mass pieces
    [x_pcs, F_pcs] = splitMass(update_params.num_mass_pcs, ...
                               update_params.x_list, g);
    M = F_pcs(end);

    % Center-of-mass coordinates of mass pieces
    %x_com = comCoordinates(update_params.x_list(1), x_pcs);
    
    % Calculate the Wasserstein gradient at x_com
    % Use Fie library to solve the corresponding Fredholm equation

    % Define the FIE-compatible kernel and RHS functions
    g_hat = fft(g);
    function res = K(k_, q_)
        [res,~,~,~,~,~,~] = wfpr_1d(k_, q_, update_params.k_list, ...
                                    fftshift(g_hat), fftshift(sqrtI), ...
                                    update_params.h);
    end
    function res = RHS(k_)
        res = wfrhs_1d(k_, update_params.k_list, ...
                        fftshift(g_hat), fftshift(sqrtI));
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
    w_com = w_func(x_rpt);
    x_com_new = x_com + update_params.h * w_com;
    % Wrap the values around (everything that is too far on the
    % right, comes back on the left side)
    x_com_new_bd = makeModular(update_params.x_list(1), ...
                               update_params.x_list(end), ...
                               x_com_new);
    % Sort new coordinates in ascending order
    x_com_new_s = sort(x_com_new_bd);
    x_rpt_new = rptCoordinates(update_params.x_list(1), ...
                               x_com_new_s);
    if any(x_rpt_new ~= sort(x_rpt_new))
        disp('Awww...')
        
        figure
        hold on
        plot(x_rpt_new, 'x')
        plot(sort(x_rpt_new), 'o')
        hold off
    end

    % Our function will be interpolated between x_list(1) and
    % x_list(end) --- pad current coordinates, if needed
    if x_rpt_new(1) > update_params.x_list(1)
        x_rpt_new = [update_params.x_list(1) x_rpt_new];
        F_rpt = [0 F_pcs];
    end
    if x_rpt_new(end) < update_params.x_list(end)
        x_rpt_new = [x_rpt_new update_params.x_list(end)];
        F_rpt = [F_pcs M];
    end
    
    F_func_new = griddedInterpolant(x_rpt_new, ...
                                    F_rpt, ...
                                    'linear', 'none');
    figure
    plot(x_com_new)
    plot(x_com_new_s)
    plot(x_rpt_new)
    
    % Calculate the new density
    m = F_pcs(1);
    g_rpt_new = zeros(size(x_rpt_new));
    g_rpt_new(1) = 2 * m / (x_rpt_new(2) - update_params.x_list(1));
    g_rpt_new(end) = m / (update_params.x_list(end) - x_rpt_new(end-1));
    for i = 2:1:(length(g_rpt_new) - 1)
        g_rpt_new(i) = 2 * m / (x_rpt_new(i + 1) - ...
                                x_rpt_new(i - 1));
    end
    if any(g_rpt_new < 0)
        disp('Density less than 0!');
    end

    
    % Leftmostpoint is a dummy to ensure non-NaN values for the
    % leftmost part of the interpolation
    g_func_new = griddedInterpolant(x_rpt_new, ...
                                    g_rpt_new, ...
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
        hold on;
        subplot(2,1,1)
        plot(x_com, F_pcs, 'xb'); 
        plot(x_com_new, F_pcs, 'og');
        plot(x_rpt_new, F_pcs, '-r');
        yyaxis right
        plot(x_rpt_new, g_rpt_new, '-b');
        plot(update_params.x_list, g_new, '.g');
        subplot(2,1,2);
        plot(real(w_hat_s));
        plot(imag(w_hat_s));
        hold off;
        % Debugging info -- end
end
