function [val, suffix] = num2eng(num, sigdig)
%NUM2ENG Convert numbers to a engineering output using SI prefixes.
%   [V, P] = NUM2ENG(X) converts a number X into string representation in
%   which string V carries the numerical value and string P the SI prefix
%   to form a numerical value corresponding to the input number X.
%
%   [V, P] = NUM2ENG(X, S) converts a number X into string representation
%   in which string V carries the numerical value and string P the SI
%   prefix to form a numerical value corresponding to the input number X.
%   The output number V is rounded to the number of significant digit
%   specified by S.
%
%   Example:
%       [V, E] = num2eng(10 ^ pi, 4) produces outputs
%
%       V =
%
%       1.385
%
%
%       E =
%
%       k
%
% if input does not exist or is empty, throw an exception

suffix={'y', 'z', 'a', 'f', 'p', 'n', '\mu', 'm', '', 'k', 'M', 'G', ...
    'T', 'P', 'E', 'Z', 'Y'};

coef = - 24 : 3 : 24;
expon = 10 .^ coef;
grth = sign(abs(num) - expon);
in = find(grth == 1, 1, 'last');
suffix = suffix{in};

mcoef = expon(in);

val = num / mcoef;

if nargin == 1
    val = num2str(val);
    return
end

if sigdig < 1
    sigdig = 1;
end
expon = [1, 10, 100];
grth = sign(abs(val) - expon);
in = find(grth == 1, 1, 'last');
val = round(val * 10 ^ (sigdig - in)) / 10 ^ (sigdig - in);
val = num2str(val);
if strfind(val, '.')
    val = horzcat(val, repmat('0', 1, sigdig - length(val) + 1));
else
    if sigdig > length(val)
        val = horzcat(val, '.', repmat('0', 1, sigdig - length(val)));
    end
end