function genes_aln(process)
global last_process
global folder
global genes_seq
global genes_name
global chromosome
global chr
global genes_pos
global fatal_error

%% Check if there is any genes file
if isempty(genes_name)
    return
end

%% Load saved target sequences
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'Loaded gene sequences');
    write_log(process);
    return
end

%% Get chromosome name of the target
c = zeros(length(chromosome), 1);
for i = 1 : length(chromosome)
    c(i) = chromosome{i}(1);
end
c = median(c);

%% Align genes to the target chromosome
write_log(process, 'Align genes to the target chromosome');
exec_blast(process, genes_seq, genes_name, chr(c));
delete([folder 'target.fa']);
if fatal_error > 0; return; end

%% Load the alignment data
write_log(process, 'Read the alignment output');
[gbit_score, gquery_from, gquery_to, ghit_from, ghit_to, gchromosome] = ...
    load_blast(process, genes_name, chr(c), 1);
if fatal_error > 0; return; end

%% Generate the positions of the genes on the chromosome
write_log(process, 'Acquire the gene positions on the chromosome');
genes_pos = get_1st_alns(process, genes_name, genes_seq, gquery_from, gquery_to, ghit_from, ghit_to, ...
    gchromosome, chr(c));
if fatal_error > 0; return; end

%% Save output
write_log(process, 'Genes positions acquired');
save([folder 'mats/fp' num2str(process) '.mat'], 'genes_pos');

%% Process finished
write_log(process);