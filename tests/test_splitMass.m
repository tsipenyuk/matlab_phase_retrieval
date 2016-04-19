% test molecule

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

x_list = m.grid.xAxes{1};
f_list = m.density;

% Calculate the CDF for plotting purposes
F_list = f_list;
for i = 1 : 1 : length(f_list) - 1
    F_list(i+1) = F_list(i) + F_list(i+1) * (x_list(i+1) - x_list(i));
end

% Total mass of the density
M = F_list(end);

num_mass_pieces = 10;
x_pcs = splitMass(num_mass_pieces, x_list, f_list);
F_pcs = linspace(M / num_mass_pieces, M, num_mass_pieces);

x_com = comCoordinates(0, x_pcs);
F_com = comCoordinates(0, F_pcs);

% Plot results
figure
hold on
plot(x_list, f_list);
plot(x_list, F_list);
plot(x_pcs, F_pcs, 'cx');
plot(horzcat(0,x_pcs), horzcat(0,F_pcs), 'c--');
plot(x_com, F_com, 'bo');


legend('f', '\int f', 'Mass pcs: right end coords', '--//-- (lin. interpolation)', ['Mass ' ...
                    'pcs: center-of-mass coords.'], 'location', 'northwest')
hold off

fig = gcf;
saveas(fig, 'tests/test_splitMass.png');
disp(['Script test_splitMass finished successfully.'])
disp(['Plotted results were ' ...
      'saved in the file test_splitmass.png.']);
