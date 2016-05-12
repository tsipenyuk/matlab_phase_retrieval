% Test the integration and derivative routines
f = randomPositivePeaks_1d(5); % Random function
N_pts = 30;
x_list = linspace(0, 1, N_pts);

f_list = circshift(f(x_list), [N_pts/2, 0]); % Remove discontinuities on
                                      % the boundaries

F_list = indefiniteIntegral(x_list, f_list);
g_list = derivative_v1(x_list, F_list);
H_list = indefiniteIntegral(x_list, g_list);
disp(['Old mass: ' num2str(F_list(end)) '; new mass: ' ...
      num2str(H_list(end)) '; difference: ' ...
      num2str(F_list(end) - H_list(end))])

figure
hold on
plot(x_list, f_list, 'o')
plot(x_list, F_list)
plot(x_list, g_list, 'x')
plot(x_list, H_list, 'o')
hold off