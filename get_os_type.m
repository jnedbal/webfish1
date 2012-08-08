function OS = get_os_type(process)

[st, w] = unix('uname');
error_msg(process, st, {'Failed getting the system type', w});
if st ~= 0 || length(w) < 1
    OS = 0;
else
    os_type = w(1 : end - 1);
    if strcmpi(os_type, 'darwin')
        OS = 1;     % MAC OS X computer
	elseif strcmpi(os_type, 'linux')
        OS = 2;     % Linux computer
	else
        OS = 0;     % Unknown operating system
    end
end
