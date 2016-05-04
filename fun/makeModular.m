function x_new = makeModular(x_min, x_max, x_list)
% MAKEMODULAR 'Wrap around' entries of a list between two values.
%
% Caveat. Is VERY inefficient if maximal values of x_list are much
% greater than x_max, or if minimal values of x_list are much less
% than x_min.
    if x_min >= x_max
        error('Input error: x_min must be smaller than x_max.')
    end
    
    while any(x_list > x_max)
        x_list(x_list > x_max) = x_list(x_list > x_max) - x_max + x_min;
    end
    
    while any(x_list < x_min)
        x_list(x_list < x_min) = x_list(x_list < x_min) ...
            + x_max - x_min;
    end
    
    x_new = x_list;
end
