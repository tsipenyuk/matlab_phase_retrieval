function [g_new, error] = phaseUpdate(g, A, h)
% phaseUpdate - Gradient flow in the phase space
%
% Synopsis ([]s are optional)
%   [g_new, error] = phaseUpdate(g, A, h)
%
% Description
%  Calculate the Fourier phase $\phi$ of the current approximation
%  $g$. Then, use gradient descent method to change the Fourier
%  phase by minimizing the non-negativity functional
%  $$
%  E[g[\phi]] = \int_{\{g[\phi]<0\}} g[\phi](x)^2 \, dx, 
%  $$
%  where
%  $$
%  g[\phi] = \frac 1 (2\pi)^d \int e^{ikx} A(k) e^{i\phi(k)} dk.
%  $$
%  
%  In other words, this algorithm performs the following updates:
%  1) calculate phase $\phi$;
%  2) update phase as follows:
%     $\phi_{new} = \phi - h * \frac{delta}{delta\phi}E[phi]$;
%  3) symmetrize $\phi$ to reduce numerical errors;*)
%  4) calculate new density $g_{new}$.
%
%  *) The phase must be symmetrical, since $g$ and $g_new$ must be
%  real-valued. 
%
% Inputs
%   (ndarray) g  Current approximation
%   (ndarray) A  Desired Fourier modulus (square root of the
%                measured intensity)
%   (double)  h  Gradient descent step length
%
% Outputs
%   (ndarray) g_new  New approximation
%   (double)  error  Error functional of the new approximation
%                    $E[g_{new}]$, cf. Description
%
% Examples
%   %Example description
%   x = 1;
%
% See also
%   er
%   
% References
%   J. R. Fienup, “Phase retrieval algorithms: a comparison,” 
%       Applied Optics, vol. 21, p. 2762-2763, end of Section IV,
%       1982.
%
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase Retrieval Sandbox root folder.
%
% Changes
%   2016-07-29 First version
    if nargin == 2
        h = 0.05;
    end
    
    phi = fftshift(angle(fftn(g)));

    % Calculate the gradient
    g_id = 1.0 * (g < 0);
    g_mod = fftshift(ifftn(g_id .* (abs(g) - g)));
    ksi = cos(phi) .* imag(g_mod) + sin(phi) .* real(g_mod);
    phi_grad = ksi .* fftshift(A);

    % Update the phase
    phi_new = phi - h * phi_grad;

    % Symmetrize the new phase
    sz = size(phi_new);
    phi_flip = -phi_new(end:-1:1);
    phi_flip = reshape(phi_flip, sz);
    phi_new_is = ifftshift(0.5 * (phi_new + phi_flip));
    
    g_new = real(ifftn(A .* exp(1j * phi_new_is)));
    % Calculate the error
    % To be comparable with modulus errors, one must account for
    % the scaling factor in the discrete Plancherel thm. Cf., for
    % example,
    % https://en.wikipedia.org/wiki/Discrete_Fourier_transform,
    % section 'The Planchele thm. and Parseval thm.' (accessed 29-07-2016).
    error = eP(g_new)*numel(g_new);
end