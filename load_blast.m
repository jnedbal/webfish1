function [bit_score, query_from, query_to, hit_from, hit_to, ...
    chromosome] = load_blast(process, target_name, chr, genes)
global folder
global fatal_error
global blast_cutoff

bit_score = cell(length(target_name), 1);
query_from = bit_score;
query_to = bit_score;
hit_from = bit_score;
hit_to = bit_score;
chromosome = bit_score;
for i = 1 : length(chr)
    write_log(process, ['Import BLAST alignment for "' chr{i} '"']);
    for j = 1 : length(target_name)
        % Select alignments with bit score of more than 200
        check_file_presence(process, ...
            [folder chr{i} '.' num2str(j) '.out']);
        if fatal_error > 0; return; end
        unix_cmd = ['awk -F"\t" '' $NF > ' num2str(blast_cutoff) ' '' ' ...
            folder chr{i} '.' num2str(j) '.out > ' folder chr{i} '.bss'];
        [st, w] =  unix(unix_cmd);
        if st ~= 0
            fatal_error = 1;
            fatal_msg(process, ...
                {['Failed reading bit scores of alignments from "' ...
                chr{i} '"'], w});
            return
        end
        % Import bit score, query from, query to, hit from, hit to
        bit_score{j} = vertcat(bit_score{j}, import_blast(process, 12, ...
            [folder chr{i} '.bss'], [folder 'tmp'], '%f32'));
        query_from{j} = vertcat(query_from{j}, import_blast(process, 7, ...
            [folder chr{i} '.bss'], [folder 'tmp'], '%u32'));
        query_to{j} = vertcat(query_to{j}, import_blast(process, 8, ...
            [folder chr{i} '.bss'], [folder 'tmp'], '%u32'));
        hit_from{j} = vertcat(hit_from{j}, import_blast(process, 9, ...
            [folder chr{i} '.bss'], [folder 'tmp'], '%u32'));
        hit_to{j} = vertcat(hit_to{j}, import_blast(process, 10, ...
            [folder chr{i} '.bss'], [folder 'tmp'], '%u32'));
        delete([folder chr{i} '.bss'])
        % Get number of loaded hits
        nr_sel_lines = length(hit_to{j}) - length(chromosome{j});
        chromosome{j} = vertcat(chromosome{j}, ...
            uint8(ones(nr_sel_lines, 1) * i));
        % Get number of hits
        unix_cmd = ['wc -l ' folder chr{i} '.' num2str(j) ...
            '.out | awk ''{print $1}'' > ' folder 'tmp'];
        [st, w] =  unix(unix_cmd);
        if st ~= 0
            fatal_error = 1;
            fatal_msg(process, {['Failed acquiring number of bit ' ...
                'scores of alignments from "' chr{i} '"'], w});
            return
        end
        [fid, w] = fopen([folder 'tmp'], 'r');
        if fid == -1
            fatal_error = 1;
            fatal_msg(process, {'Failed loading tmp file', w});
            return
        end
        nr_lines = str2double(fgetl(fid)) - 4;
        fclose(fid);
        if nr_lines <= nr_sel_lines && nr_lines > 1 && nargin < 4
            fatal_error = 1;
            fatal_msg(process, {['Not enough alignments for "' ...
                regexprep(target_name{j}, '\\', '\\\') '" on "' chr{i} ...
                '"']});
            return
        end
        delete([folder 'tmp']);
        delete([folder chr{i} '.' num2str(j) '.out']);
    end
end

%% Sort alignments
for j = 1 : length(target_name)
    [bit_score{j}, in] = sort(bit_score{j}, 'descend');
    query_from{j} = query_from{j}(in);
    query_to{j} = query_to{j}(in);
    hit_from{j} = hit_from{j}(in);
    hit_to{j} = hit_to{j}(in);
    chromosome{j} = chromosome{j}(in);
end