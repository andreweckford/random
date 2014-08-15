function y = randig(lambda,mu,rows,cols,tolerance)

% RANDIG  Returns a scalar, vector, or matrix of inverse Gaussian
% distributed random variables.
%
%   c = RANDIG(lambda,mu) returns a scalar random variable.
%   c = RANDIG(lambda,mu,n) returns an n-by-n matrix.
%   c = RANDIG(lambda,mu,m,n) returns an m-by-n matrix.
%   c = RANDIG(lambda,mu,m,n,tol) returns an m-by-n matrix, and specifies
%   the tolerance in inverting the IG CDF (default: 1e-4)
%
% In each case the inverse Gaussian parameters are given by lambda and mu.
% For a description see:
% http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution

    function y = cdf_ig(x,lambda,mu)
        
        y = (1 + erf(sqrt(lambda./(2*x)).*(x/mu - 1)))./2;
        y = y + exp(2*lambda/mu)*(1 + erf(-1*sqrt(lambda./(2*x)).*(x/mu + 1)))./2;
        
    end

default_tol = 1e-4;

switch nargin
    case 4
        tolerance = default_tol;
    case 3
        cols = rows;
        tolerance = default_tol;
    case 2
        cols = 1;
        rows = 1;
        tolerance = default_tol;
end

% uniformly distributed random variable
z = rand(rows,cols);

y = zeros(rows,cols);

for r = 1:rows
    
    for c = 1:cols
        
        lohi = [0 1];
        
        % invert the IG CDF
        % exponential search until something greater than z is found
        
        while (cdf_ig(lohi(2),lambda,mu) < z(r,c))
            lohi(2) = lohi(2) * 2;
        end
        
        % search the range until found
        while (abs(lohi(2) - lohi(1)) > tolerance)
            mid = (lohi(1) + lohi(2))/2;
            if (cdf_ig(mid,lambda,mu) < z(r,c))
                lohi(1) = mid;
            else
                lohi(2) = mid;
            end
        end
        
        y(r,c) = (lohi(1) + lohi(2))/2;
        
    end
    
end

end