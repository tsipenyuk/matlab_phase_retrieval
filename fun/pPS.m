function res = pPS(g, S)
% pS Positivity and support projection
%    Replace all data of g by zero at S
    res = real(abs(g .* (g >= 0)));
    res(~S) = 0;
end