function res = positivityProjection(g)
% POSITIVITY_PROJECTION Replace negative entries of g by 0
    res = abs(g .* (g >= 0));
end