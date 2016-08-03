% The following script plots fiber diffraction data for 2c0w
% filamentous bacteriophage, see 
% http://www.rcsb.org/pdb/explore/explore.do?structureId=2c0w
%
% Although the data file was slightly modified to ensure correct
% parsing using cif2mat, no data was altered.
data = cif2mat('2c0w-sf.cif')
r = data.refln.pdbx_fiber_coordinate;
z = data.refln.pdbx_fiber_layer;
F = data.refln.pdbx_fiber_F_meas_au; % Square root of measured intensity

% Plot F where F~=0
c = (F~=0); 
figure; 
scatter(r(c), z(c), ...  
        10, ...   % scatter circle radius
        F(c), ... % scatter circle color
        'filled');

% Uncomment the following lines to plot the molecular structure
%pdbData = pdbread('2c0w.pdb');
%molviewer(pdbData);