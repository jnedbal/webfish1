function config_out = read_config(hit_string, process)
global config_file
global folder

%% Read value from the config file
unix_cmd = ...
    ['grep ''^<' hit_string '>'' -i ' config_file ' > ' folder 'tmp'];
[st, w] = unix(unix_cmd);
error_msg(process, st, {['Failed reading file ' config_file], w});
[fid, w] = fopen([folder 'tmp']);
    if fid ~= -1
        st = 0;
    else
        st = 1;
    end
    error_msg(process, st, {['Failed loading config file line: "' ...
        hit_string '"'], w});
    z = 0;
    if feof(fid) == 1
        error_msg(process, st, {['Failed loading config file line' ...
            hit_string], w});
    end
    config_out = {''};
    while feof(fid) == 0
        z = z + 1;
        tline = fgetl(fid);
        tp1 = strfind(tline, ['<' hit_string '>']);
        tp2 = length(hit_string) + 2;
        tp3 = strfind(tline, ['</' hit_string '>']);
        config_out{z} = tline(tp1 + tp2 : tp3 - 1);
        write_log(process, [hit_string ': ' config_out{z}]);
    end
    if z == 1
        config_out = cell2mat(config_out);
    end
fclose(fid);
delete([folder 'tmp']);