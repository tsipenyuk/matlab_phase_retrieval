function g = initialGuess(sqrtI)
% initialGuess - yields some non-negative real-valued density
%
% Synopsis
%   g = initialGuess(sqrtI)
%
% Description
%   Generate a non-negative real-valued density using the square
%   root of the intensity measurement. Specifically, set NaN-values
%   to 0 and calculate the inverse FT of the autocorrelation.
%
% Inputs
%   (array) sqrtI ndgrid or meshgrid array, square root of the
%                 measured intensity
%
% Outputs
%   (array) g     density guess
%
% Examples
%   %% 3D, two gaussians
%   [H, K, L] = ndgrid([-20:0.2:20], [-20:0.2:20],  [-20:0.2:20]);
%   g = exp(-x1.^2 - x2.^2 - x3.^2);
%   shift = fix(length(g)/4);
%   F = circshift(g, [0, shift, shift]) + ...
%       circshift(g, [0, -shift, 0]);
%   plot_slice(H,K,L,F,[-1.2 .8 2],2,[-2 -.2])
%
% Authors
%   Arseniy Tsipenyuk <tsipenyu(at)ma.tum.de>
%
% License
%   See Phase nRetrieval Sandbox root folder.
%
% Changes
%   2016-06-15  First Edition

    sqrtI_patch = sqrtI;
    sqrtI_patch(isnan(sqrtI)) = 0;
    g = sqrt(pP(ifftn(sqrtI_patch.^2)));
end