function import_target(process)
global last_process
global folder
global target_seq
global target_name
global target_file
global blast_length
global fatal_error

%% Load saved target sequences
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'Loaded target sequence');
    write_log(process);
    return
end

%% Import target sequences from fasta file
target_file = read_config('tar_seq', process);
file = dir(target_file);
if isempty(file)
    fatal_error = 1;
    fatal_msg(process, {['Target file "' target_file '" does not exist']});
    return
end
[target_name target_seq] = import_fasta_sequence(target_file, process);
if fatal_error > 0; return; end

tsl = sum(cellfun(@length, target_seq));
if tsl > blast_length
	fatal_error = 1;
    fatal_msg(process, {['Target sequence is too long (' num2str(tsl) ...
        ' bp)'], ['Maximum target sequence length is ' ...
        num2str(blast_length) ' bp']});
    return
end

save([folder 'mats/fp' num2str(process) '.mat'], 'target_name', ...
    'target_seq', 'target_file');
write_log(process, 'Imported target sequence');

%% Process finished
write_log(process);