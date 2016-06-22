function [x, y, z, f] = pA(g, N)
% pA Singleton atomicity projection

    [X, Y, Z] = ndgrid(1:size(g,1),1:size(g,2),1:size(g,3));
    F = g;
    
    x_full = reshape(X, [], 1);
    y_full = reshape(Y, [], 1); 
    z_full = reshape(Z, [], 1);
    f_full = reshape(F, [], 1);

    x = zeros([1, N]);
    y = zeros([1, N]);
    z = zeros([1, N]);
    f = zeros([1, N]);
    
    for iVal = 1:N
        % Find the voxel with the maximal value; record it, remove
        % it from future iterations
        [~, argmax] = max(f_full);
        x(iVal) = x_full(argmax);
        y(iVal) = y_full(argmax);
        z(iVal) = z_full(argmax);
        f(iVal) = f_full(argmax);
        f_full(argmax) = NaN;
    end
end