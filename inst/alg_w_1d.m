% Instantiate error-reduction algorithm assuming that an instance
% of molecule called 'm' is already specified.
update_params.bla = 1
w1 = phaseRetrievalAlgorithm(@wasserstein_1d, m, update_params);
disp('Created an instance "w1" of class phaseRetrievalAlgorith.m')