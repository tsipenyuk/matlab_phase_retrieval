function f = randomPositivePeaks_1d(n);
% RANDOMPOSITIVEPEAKS_1d Generate a function from n Legendre polynomials
%    
% Auxilary routine to generate functions with some peaks between 0
% and 1 and decaying towards zero else.
    coeff = randi([-50,50],1,n);
    coeff = coeff / 10;
    function y = func(x)
        y = zeros(size(x));
        for i = 1:1:n
            y = y + coeff(i) * legendreP(i, x);
        end
        y = exp(-abs(y));
    end
    f = @func;
end