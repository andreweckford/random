function y = randlevy(s,rows,cols)

% RANDLEVY  Returns a scalar, vector, or matrix of Lévy distributed random variables.
%
%   c = RANDLEVY(s) returns a scalar random variable.
%   c = RANDLEVY(s,n) returns an n-by-n matrix
%   c = RANDLEVY(s,m,n) returns an m-by-n matrix.
%
% In each case the scale parameter is given by s.
% For a description see:
% http://en.wikipedia.org/wiki/Levy_distribution

switch nargin
    case 2
        cols = rows;
    case 1
        cols = 1;
        rows = 1;
end

y = s./ (2*(erfcinv(rand(rows,cols)).^2));

