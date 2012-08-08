function gen_un_frags(process)
global last_process
global target_name
global uniq_frag_pos
global uniq_frag_pos_over
global uniq_frag_seq
global uniq_frag_seq_over
global uniq_frag_name
global uniq_reg
global folder
global target_seq
global min_pcr_size
global max_pcr_size
global log_file
global pcr_overlap
global enzymes
global fatal_error

%% Check if fragments are produced
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'PCR fragments for unique sequences are ready');
    write_log(process);
    return
end

alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';

%% Define unique fragments positions
uniq_frag_pos = cell(length(target_name), 1);
uniq_frag_pos_over = cell(length(target_name), 1);
uniq_frag_seq = cell(length(target_name), 1);
uniq_frag_seq_over = cell(length(target_name), 1);

%% Open a file to write the fragment sequences
[fid, w] = fopen([log_file(1 : end - 3) 'unique_segs.txt'], 'w');
if fid == -1
    fatal_error = 1;
    fatal_msg(process, {'Failed writing file with unique fragments', w});
    return
end

%% Split contiguous regions into fragments
for i = 1 : length(target_name)
    fprintf(fid, [repmat('*', 1, length(target_name{i}) + 8) '\n']);
    fprintf(fid, ['*   ' target_name{i} '   *\n']);
    fprintf(fid, [repmat('*', 1, length(target_name{i}) + 8) '\n\n']);
    m = 0;  % Variable numbering unique regions in given target
    for j = 1 : size(uniq_reg{i}, 2)
        % Length of given unique region
        uniq_reg_len = uniq_reg{i}(2, j) - uniq_reg{i}(1, j) + 1;
        % Process long enough unique regions
        if uniq_reg_len >= min_pcr_size
            % Number of fragments to cover unique region
            nr_pcrs = ceil(uniq_reg_len / max_pcr_size);
            for k = 1 : nr_pcrs
                m = m + 1;
                % Starting position of PCR fragment
                uniq_frag_pos{i}(1, m) = ...
                    round(uniq_reg_len * (k - 1) / nr_pcrs) + ...
                    uniq_reg{i}(1, j);
                % Ending position of PCR fragment
                uniq_frag_pos{i}(2, m) = ...
                    round(uniq_reg_len * k / nr_pcrs) - 1 + ...
                    uniq_reg{i}(1, j);
                % Starting position of PCR fragment with overlap; it is
                % trimmed so that it does not span further than the unique
                % region
                uniq_frag_pos_over{i}(1, m) = max([uniq_reg{i}(1, j), ...
                    uniq_frag_pos{i}(1, m) - pcr_overlap]);
                % Ending position of PCR fragment with overlap; it is
                % trimmed so that it does not span further than the unique
                % region
                uniq_frag_pos_over{i}(2, m) = min([uniq_reg{i}(2, j), ...
                    uniq_frag_pos{i}(2, m) + pcr_overlap]);
                % Sequence of the PCR fragment
                uniq_frag_seq{i}{m} = ...
                    target_seq{i}((uniq_frag_pos{i}(1, m) : ...
                    uniq_frag_pos{i}(2, m)));
                % Sequence of the PCR fragment with overlap in lower case
                % letters
                uniq_frag_seq_over{i}{m} = ...
                    lower(target_seq{i}((uniq_frag_pos_over{i}(1, m) : ...
                    uniq_frag_pos_over{i}(2, m))));
                % Find coordinates of non-overlapping PCR fragment
                in = (uniq_frag_pos{i}(1, m) - ...
                    uniq_frag_pos_over{i}(1, m)) + ...
                    (1 : (uniq_frag_pos{i}(2, m) - ...
                    uniq_frag_pos{i}(1, m) + 1));
                % Capitalize letters of the non-overlapping part of the PCR
                % fragment
                uniq_frag_seq_over{i}{m}(in) = ...
                    upper(uniq_frag_seq_over{i}{m}(in));
                % Give fragment a name
                uniq_frag_name{i}{m} = ['U' num2str(j) alph(k)];
                fprintf(fid, ['Target sequence: ' target_name{i} '\n']);
                fprintf(fid, ['Segment name: ' uniq_frag_name{i}{m} '\n']);
                fprintf(fid, ['Segment from: ' ...
                    num2str(uniq_frag_pos{i}(1, m)) ' bp\n']);
                fprintf(fid, ['Segment to: ' ...
                    num2str(uniq_frag_pos{i}(2, m)) ' bp\n']);
                fprintf(fid, ['Segment length: ' ...
                    num2str(uniq_frag_pos{i}(2, m) - ...
                    uniq_frag_pos{i}(1, m) + 1) ' bp\n']);
                for n = 1 : length(enzymes)
                    [e sites] = ...
                        rebasecuts(uniq_frag_seq_over{i}{m}, enzymes{n});
                    fprintf(fid, ['Number of ' enzymes{n} ' sites: ' ...
                        num2str(length(sites)) '\n']);
                end
                fprintf(fid, ['Segment sequence: ' uniq_frag_seq{i}{m} ...
                    '\n']);
                fprintf(fid, ['Segment sequence with overhang: ' ...
                    uniq_frag_seq_over{i}{m} '\n\n']);
            end
        end
    end
end
fclose(fid);

%% Save output
save([folder 'mats/fp' num2str(process) '.mat'], 'uniq_frag_pos', ...
    'uniq_frag_pos_over', 'uniq_frag_seq', 'uniq_frag_seq_over', ...
    'uniq_frag_name');
write_log(process, 'Generated PCR fragments for unique sequences');

%% Process finished
write_log(process);