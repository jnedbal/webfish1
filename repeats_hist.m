function repeats_hist(process)
global last_process
global folder
global aln_pos
global query_from
global query_to
global hit_from
global hit_to
global target_name
global nonuniq_vec
global nonuniq_reg
global chromosome

%% Load saved target sequences
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'Histograms of repeats count generated');
    write_log(process);
    return
end

%% Produce histograms of number of repeats confined to non-unique regions
nonuniq_vec = cell(1, length(target_name));
for i = 1 : length(target_name)
    write_log(process, ['Generating histograms of repeats within '...
        'non-unique regions of target region ''' target_name{i} '''']);
    nonuniq_vec{i} = uint32(zeros(1, aln_pos(i, 2) - aln_pos(i, 1) + 1));
    in = find(chromosome{i} == chromosome{i}(1));
	qf = int32(query_from{i}(in));
    qt = int32(query_to{i}(in));
	hf = int32(hit_from{i}(in)) - int32(aln_pos(i, 1) - 1);
    ht = int32(hit_to{i}(in)) - int32(aln_pos(i, 1) - 1);
    for j = 1 : size(nonuniq_reg{i}, 2)
        in = find(qf >= nonuniq_reg{i}(1, j) & ...
            qf <= nonuniq_reg{i}(2, j) & ...
            qt >= nonuniq_reg{i}(1, j) & ...
            qt <= nonuniq_reg{i}(2, j) & ...
            ((hf >= nonuniq_reg{i}(1, j) & ...
            hf <= nonuniq_reg{i}(2, j)) | ...
            (ht >= nonuniq_reg{i}(1, j) & ...
            ht <= nonuniq_reg{i}(2, j))));
        for k = 1 : length(in)
            h = sort([qf(in(k)), qt(in(k))]);
            nonuniq_vec{i}(h(1) : h(2)) = nonuniq_vec{i}(h(1) : h(2)) + 1;
            h = sort([hf(in(k)), ht(in(k))]);
            h(1) = max([1, h(1)]);
            h(2) = min([length(nonuniq_vec{i}), h(2)]);
            nonuniq_vec{i}(h(1) : h(2)) = nonuniq_vec{i}(h(1) : h(2)) + 1;
        end
    end
end

%% Make regions of histograms which are similar to distant regions zero
for i = 1 : length(target_name)
    write_log(process, ['Delete histogram reqions that are not unique ' ...
        'in target sequence ''' target_name{i} '''']);
    %in = find(chromosome{i} == chromosome{i}(1));
	qf = int32(query_from{i});
    qt = int32(query_to{i});
	hf = int32(hit_from{i}) - int32(aln_pos(i, 1) - 1);
    ht = int32(hit_to{i}) - int32(aln_pos(i, 1) - 1);
    for j = 1 : size(nonuniq_reg{i}, 2)
        % Find indices of queries that overlap with the non-unique region
        inq = (qf >= nonuniq_reg{i}(1, j) & ...
            qf <= nonuniq_reg{i}(2, j)) | ...
            (qt >= nonuniq_reg{i}(1, j) & ...
            qt <= nonuniq_reg{i}(2, j));
        % Find indices of queries that lie outside of the non-unique region
        ink = (qf < nonuniq_reg{i}(1, j) | ...
            qf > nonuniq_reg{i}(2, j)) & ...
            (qt < nonuniq_reg{i}(1, j) | ...
            qt > nonuniq_reg{i}(2, j));
        % Find indices of hits that overlap with the non-unique region
        inh = (((hf >= nonuniq_reg{i}(1, j) & ...
            hf <= nonuniq_reg{i}(2, j)) | ...
            (ht >= nonuniq_reg{i}(1, j) & ...
            ht <= nonuniq_reg{i}(2, j)))) & ...
            chromosome{i} == chromosome{i}(1);
        % Find indices of hits on different chromosome
        inc = chromosome{i} ~= chromosome{i}(1);
        % Find indices of hits on the same chromosome but outside of the
        % non-unique region
        ins = ((hf < nonuniq_reg{i}(1, j) | ...
            hf > nonuniq_reg{i}(2, j)) & ...
            (ht < nonuniq_reg{i}(1, j) | ...
            ht > nonuniq_reg{i}(2, j))) & ...
            chromosome{i} == chromosome{i}(1);
        in = find((inq & inc) | (inq & ins) | (inh & ink));
        for k = 1 : length(in)
            h = sort([qf(in(k)), qt(in(k))]);
            nonuniq_vec{i}(h(1) : h(2)) = 0;
        end
    end
end

%% Save output
write_log(process, 'Histogram of numbers of similar regions done');
save([folder 'mats/fp' num2str(process) '.mat'], 'nonuniq_vec');

%% Process finished
write_log(process);