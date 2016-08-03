data = cif2mat('data/1hgv/1hgv-sf.cif')
r = data.refln.fiber_coordinate;
z = data.refln.fiber_layer;
F = data.refln.fiber_F_meas_au;

% Plot F where F~=0
c = (F~=0); 
figure; 
scatter(r(c), ... % Radial coordinate in Fourier space
        z(c), ... % 
        50, ...   % Plot -- scatter circle radius
        F(c), ... % Plot -- scatter circle color
        'filled');