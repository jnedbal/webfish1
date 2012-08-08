function fatal_msg(process, text)
global log_file
global tasks
global webserver
global fatal_error

if fatal_error == 2
    return
end

fatal_error = 2;
if process > 0
    [fid, w] = fopen([log_file(1 : end - 3) 'fatal'], 'w');
    if fid == -1
        st = 1;
        error_msg(process, st, {'Failed writing fatal error file', w});
    else
        fprintf(fid, 'Fatal Error!\n');
        fprintf('Fatal Error!\n');
        for i = 1 : length(text)
            fprintf(fid, [text{i} '\n']);
            fprintf([text{i} '\n']);
        end
        fclose(fid);
    end
end

write_log(process, '<span style="color:#FF0000">Fatal Error!</span>');
for i = 1 : length(text)
    write_log(process, ['<span style="color:#FF0000">' text{i} '</span>']);
end
producehtml(process, '</div><h1>Fatal Error!</h1>', 'Processing Your Query');
producehtml(process, '<h1>Your project has been removed</h1>', ...
    'Processing Your Query');
producehtml(process, ['<h1>Check the log above for cause of the error ' ...
    'and try resubmitting</h1>'], ...
    'Processing Your Query');
if ~isequal(webserver.html, '.')
    mailout(process);
end

if isempty(tasks{1})
    return
end
file = dir(log_file);
if ~isempty(file)
    unix_cmd = ['mv ' log_file ' done/'];
    [st, w] = unix(unix_cmd);
    error_msg(0, st, {'Failed storing log file', w});
end
fname = ['tasks/' tasks{1} '.*'];
file = dir(fname);
if ~isempty(file)
    unix_cmd = ['rm ' fname];
    [st, w] = unix(unix_cmd);
    error_msg(0, st, {'Failed deleting project input files', w});
end
%fname = ['tasks/' outtasks{1} '/*'];
fname = ['tasks/' tasks{1} '/*'];
file = dir(fname);
if ~isempty(file)
    unix_cmd = ['rm -r ' fname];
    [st, w] = unix(unix_cmd);
    error_msg(0, st, {'Failed deleting project output files', w});
end
%fname = ['tasks/' outtasks{1}];
fname = ['tasks/' tasks{1}];
file = dir(fname);
if ~isempty(file)
    unix_cmd = ['rmdir ' fname];
    [st, w] = unix(unix_cmd);
    error_msg(process, st, {'Failed deleting project directory', w});
end