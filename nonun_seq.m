function nonun_seq(process)
global last_process
global folder
global aln_pos
global chromosome
global query_from
global query_to
global hit_from
global hit_to
global target_name
global nonuniq_bin

%% Load saved target sequences
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'Loaded vector of non-unique sequences');
    write_log(process);
    return
end

%% Searching for non-unique sequences within the target region
nonuniq_bin = cell(1, length(target_name));
for i = 1 : length(target_name)
    write_log(process, ...
        ['Generating map of non-unique sequences in "' ...
        target_name{i} '"']);
    nonuniq_bin{i} = false(1, aln_pos(i, 2) - aln_pos(i, 1) + 1);
    for j = 2 : length(query_from{i})
        h = sort([query_from{i}(j), query_to{i}(j)]);
        nonuniq_bin{i}(h(1) : h(2)) = true;
    end
    in = find(((hit_from{i} >= uint32(aln_pos(i, 1)) & ...
        hit_from{i} <= uint32(aln_pos(i, 2))) | ...
        (hit_to{i} >= uint32(aln_pos(i, 1)) & ...
        hit_to{i} <= uint32(aln_pos(i, 2)))) & ...
        chromosome{i} == chromosome{i}(1));
    for j = 2 : length(in)
        h = sort([hit_from{i}(in(j)), hit_to{i}(in(j))]) - aln_pos(i, 1) + 1;
        h(1) = max([1, h(1)]);
        h(2) = min([length(nonuniq_bin{i}), h(2)]);
        nonuniq_bin{i}(h(1) : h(2)) = 1;
    end
end

%% Save output
write_log(process, 'Identified non-unique sequences');
save([folder 'mats/fp' num2str(process) '.mat'], 'nonuniq_bin');

%% Process finished
write_log(process);