% The following script plots fiber diffraction data for 1hgv
% filamentous bacteriophage, see 
% http://www.rcsb.org/pdb/explore/explore.do?structureId=1hgv
%
% Although the data file was slightly modified to ensure correct
% parsing using cif2mat, no data was altered.
data = cif2mat('1hgv-sf.cif')
r = data.refln.fiber_coordinate;
z = data.refln.fiber_layer;
F = data.refln.fiber_F_meas_au;

% Plot F where F~=0
c = (F~=0); 
figure; 
scatter(r(c), z(c), ...  
        50, ...   % scatter circle radius
        F(c), ... % scatter circle color
        'filled');

% Uncomment the following lines to plot the molecular structure
%pdbData = pdbread('1hgv.pdb');
%molviewer(pdbData);
