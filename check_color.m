function col = check_color(process, col, col_name, def_col)
if length(col) ~= 3
    col = def_col';
	error_msg(process, 1, {['Could not acquire ' col_name ], ...
        ['Using ' col_name ' of [' num2str(def_col(1)) ', ' ...
        num2str(def_col(2)) ', ' num2str(def_col(3)) ']']});
else
    def_val = repmat(def_col, 1, 3);
    lim_val = [0, 1; 0, 1; 0, 1];
    for i = 1 : 3
        col(i) = check_value(process, col(i), ...
            ['col(' num2str(i) ')'], def_val(i, :), lim_val(i, :));
    end
end