% Two gaussian peaks density

% Coordinate grid
g = mGrid('nPoints', 200, 'xMin', 0, 'xMax', 10);
x = g.xMeshgrid;

% means and widths of gaussians
m1 = 2;
m2 = 6;
w1 = 0.2;
w2 = 0.1;


density = 0.5 * exp(-(x - m1).^2 / w1) + exp(-(x - m2).^2 / w2);
m = molecule('density', density, 'grid', g);
disp('Created an instance "m" of class molecule.')
m.plot();

