function anal_par(process)
global last_process
global folder
global min_uniq_length
global min_nonuniq_length
global max_nonuniq_gap
global min_pcr_size
global max_pcr_size
global min_rep_pcr_size
global max_rep_pcr_size
global pcr_overlap
global blast_length
global enzymes

%% Load saved target sequences
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'Analysis parameters loaded');
    write_log(process);
    return
end

%% Import target sequences from fasta file
min_uniq_length = ...
    uint32(str2double(read_config('min_uniq_length', process)));
min_nonuniq_length = ...
    uint32(str2double(read_config('min_nonuniq_length', process)));
max_nonuniq_gap = ...
    uint32(str2double(read_config('max_nonuniq_gap', process)));
min_pcr_size = round(str2double(read_config('min_pcr_size', process)));
max_pcr_size = round(str2double(read_config('max_pcr_size', process)));
min_rep_pcr_size = ...
    round(str2double(read_config('min_rep_pcr_size', process)));
max_rep_pcr_size = ...
    round(str2double(read_config('max_rep_pcr_size', process)));
pcr_overlap = round(str2double(read_config('pcr_overlap', process)));
enzymes = read_config('enzymes', process);
if ~isempty(enzymes)
    enzymes = textscan(enzymes, '%s', 'delimiter', ',');
    enzymes = enzymes{1};
end


%% Check values
min_pcr_size = check_value(process, min_pcr_size, 'min_pcr_size', ...
    [4000, 500, 4000], [500, blast_length]);
max_pcr_size = check_value(process, max_pcr_size, 'max_pcr_size', ...
    [8000, min_pcr_size + 100, 8000], [min_pcr_size, blast_length]);
min_uniq_length = check_value(process, min_uniq_length, ...
    'min_uniq_length', [4000, max([1000, min_pcr_size]), 4000], ...
    [max([1000, min_pcr_size]), blast_length]);
min_nonuniq_length = check_value(process, min_nonuniq_length, ...
    'min_nonuniq_length', [4000, max([1000, min_pcr_size]), 4000], ...
    [max([1000, min_pcr_size]), blast_length]);
max_nonuniq_gap = check_value(process, max_nonuniq_gap, ...
    'max_nonuniq_gap', [1000, 100, 4000], [100, blast_length]);
pcr_overlap = check_value(process, pcr_overlap, ...
    'pcr_overlap', [300, 0, round(min_pcr_size / 2)], ...
    [0, min_pcr_size]);

%% Save output
save([folder 'mats/fp' num2str(process) '.mat'], 'min_uniq_length', ...
    'min_nonuniq_length', 'max_nonuniq_gap', 'min_pcr_size', ...
    'max_pcr_size', 'min_rep_pcr_size', 'max_rep_pcr_size', ...
    'enzymes', 'pcr_overlap');
write_log(process, 'Imported analysis parameters');

%% Process finished
write_log(process);