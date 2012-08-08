function variable = check_value(process, variable, varname, def_val, lims)

% Function that checks if imported numeric value makes sense
% 'variable' is the variable to check
% 'def_val' is a three value vector telling which is the default value if
% the 'variable' is completely wrong, if the 'variable is less then the
% lower limit, is the 'variable' is greater than the upper limit
% 'lims' is a two value vector telling the lower and upper limits
if isnan(variable) || length(variable) > 1
    error_msg(process, 1, {['Using ' varname ' of ' ...
        num2str(def_val(1))], ['Could not acquire ' varname]});
    variable = def_val(1);
elseif variable < lims(1)
    error_msg(process, 1, {['Using ' varname ' of ' ...
        num2str(def_val(2))], ...
        [varname ' is too short (' num2str(variable) ')']});
    variable = def_val(2);
elseif variable > lims(2)
    error_msg(process, 1, ...
        {['Using ' varname ' of ' num2str(def_val(3))], ...
        [varname ' is too long (' num2str(min_uniq_length) ')']});
    variable = def_val(3);
end