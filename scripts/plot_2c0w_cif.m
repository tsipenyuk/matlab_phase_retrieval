data = cif2mat('data/2c0w/2c0w-sf.cif')
r = data.refln.pdbx_fiber_coordinate;
z = data.refln.pdbx_fiber_layer;
F = data.refln.pdbx_fiber_F_meas_au;

% Plot F where F~=0
c = (F~=0); 
figure; 
scatter(r(c), ... % Radial coordinate in Fourier space
        z(c), ... % 
        10, ...   % Plot -- scatter circle radius
        F(c), ... % Plot -- scatter circle color
        'filled');