% Two gaussian peaks density

% Coordinate grid
g = mGrid('nPoints', [200 200], 'xMin', [0 0], 'xMax', [10 10], ...
          'kMin', [0 0], 'kMax', [200 200]);
x = g.xMeshgrid;

% means and widths of gaussians
m1 = {2 2};
m2 = {6 4};
w1 = {0.2 0.2};
w2 = {0.1 0.4};


density = 0.5 * exp(-(x{1} - m1{1}).^2 ./ w1{1} -(x{2} - m1{2}).^2 ./ w1{2} )...
          + exp(-(x{1} - m2{1}).^2 ./ w2{1} -(x{2} - m2{2}).^2 ./ w2{2} );
m = molecule('density', density, 'grid', g);
m.contour();

