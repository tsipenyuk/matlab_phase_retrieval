% Instantiate error-reduction algorithm assuming that an instance
% of molecule called 'm' is already specified.

update_params.x_list = m.grid.xAxes{1};
update_params.k_list = m.grid.kAxes{1};
update_params.h = 1e-1;
update_params.num_mass_pcs = 200;

% Set solution
w1 = phaseRetrievalAlgorithm(@wasserstein_1d, m, update_params);
disp('Created an instance "w1" of class phaseRetrievalAlgorith.m')

% Set initial approximation to be flat, equal 1
w1 = w1.set_density(ones(size(m.density)) * sum(m.density) / length(m.density))