function primer_par(process)
global last_process
global folder
global primer3_dir
global primer_opt_size
global primer_min_size
global primer_max_size
global primer_opt_temp
global primer_min_temp
global primer_max_temp
global primer_max_tm_diff
global primer_salt
global primer_gc_clamp
global primer_min_gc
global primer_max_gc

%% Load saved target sequences
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'Analysis parameters loaded');
    write_log(process);
    return
end

%% Import target sequences from fasta file
primer3_dir = read_config('primer3_dir', process);
if primer3_dir(end) ~= '/'
    primer3_dir = [primer3_dir '/'];
end
primer_opt_size = ...
    round(str2double(read_config('primer_opt_size', process)));
primer_min_size = ...
    round(str2double(read_config('primer_min_size', process)));
primer_max_size = ...
    round(str2double(read_config('primer_max_size', process)));
primer_opt_temp = str2double(read_config('primer_opt_temp', process));
primer_min_temp = str2double(read_config('primer_min_temp', process));
primer_max_tm_diff = ...
    str2double(read_config('primer_max_tm_diff', process));
primer_max_temp = str2double(read_config('primer_max_temp', process));
primer_salt = ...
    round(str2double(read_config('primer_salt_corrections', process)));
primer_gc_clamp = ...
    round(str2double(read_config('primer_gc_clamp', process)));
primer_min_gc = ...
    round(str2double(read_config('primer_min_gc', process)));
primer_max_gc = ...
    round(str2double(read_config('primer_max_gc', process)));

%% Check values
primer_opt_size = check_value(process, primer_opt_size, ...
    'primer_opt_size', [18, 13, 32], [12, 35]);
primer_min_size = check_value(process, primer_min_size, ...
    'primer_min_size', [15, 12, 33], [12, 35]);
primer_max_size = check_value(process, primer_max_size, ...
    'primer_max_size', [21, 14, 35], [12, 35]);
if primer_min_size >= primer_opt_size
    primer_min_size = primer_opt_size - 1;
end
if primer_max_size <= primer_opt_size
    primer_max_size = primer_opt_size + 1;
end

primer_opt_temp = check_value(process, primer_opt_temp, ...
    'primer_opt_temp', [60, 51, 71], [30, 80]);
primer_min_temp = check_value(process, primer_min_temp, ...
    'primer_min_temp', [57, 50, 70], [30, 80]);
primer_max_temp = check_value(process, primer_max_temp, ...
    'primer_max_temp', [63, 52, 72], [30, 80]);
primer_max_tm_diff = check_value(process, primer_max_tm_diff, ...
    'primer_max_tm_diff', [1, 4, 100], [1, 100]);
if primer_min_temp >= primer_opt_temp
    primer_min_temp = primer_opt_temp - 0.1;
end
if primer_max_temp <= primer_opt_temp
    primer_max_temp = primer_opt_temp + 0.1;
end

primer_salt = check_value(process, primer_salt, ...
    'primer_salt', [1, 1, 1], [0, 2]);

primer_gc_clamp = check_value(process, primer_gc_clamp, ...
    'primer_gc_clamp', [1, 0, 3], [0, primer_min_size]);

primer_min_gc = check_value(process, primer_min_gc, ...
    'primer_min_gc', [30, 10, 70], [0, 99]);
primer_max_gc = check_value(process, primer_max_gc, ...
    'primer_max_gc', [70, 30, 90], [2, 100]);
if primer_min_gc >= primer_max_gc
    primer_min_gc = primer_max_gc - 1;
end


%% Save output
save([folder 'mats/fp' num2str(process) '.mat'], 'primer3_dir', ...
    'primer_opt_size', 'primer_min_size', 'primer_max_size', ...
    'primer_opt_temp', 'primer_min_temp', 'primer_max_temp', ...
    'primer_salt', 'primer_gc_clamp', 'primer_min_gc', 'primer_max_gc', ...
    'primer_max_tm_diff');
write_log(process, 'Imported analysis parameters');

%% Process finished
write_log(process);