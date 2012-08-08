function genes_load(process)
global last_process
global folder
global genes_seq
global genes_name
global genes_file
global fatal_error

%% Load saved target sequences
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'Loaded gene sequences');
    write_log(process);
    return
end

%% Import target sequences from fasta file
genes_file = read_config('genes', process);
files = dir(genes_file);
if isempty(files)
    write_log(process, 'No genes file provided');
    save([folder 'mats/fp' num2str(process) '.mat'], 'genes_file');
    return
end

[genes_name genes_seq] = import_fasta_sequence(genes_file, process);
if fatal_error > 0; return; end

%% Save output
save([folder 'mats/fp' num2str(process) '.mat'], 'genes_file', ...
    'genes_name', 'genes_seq');
write_log(process, 'Imported gene sequences');

%% Process finished
write_log(process);