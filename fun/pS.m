function res = pS(g, S)
% pS Support projection
%    Replace all data of g by zero at S
    res = zeros(size(g));
    res(S) = g(S);
end