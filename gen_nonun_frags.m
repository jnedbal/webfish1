function gen_nonun_frags(process)
global last_process
global target_name
global nonuniq_frag_pos
global nonuniq_frag_pos_over
global nonuniq_frag_seq
global nonuniq_frag_seq_over
global nonuniq_frag_name
global nonuniq_reg
global folder
global target_seq
global min_rep_pcr_size
global max_rep_pcr_size
global hist_pos
global pcr_overlap

%% Check if fragments are produced
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'PCR fragments for non-unique sequences are ready');
    write_log(process);
    return
end

alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';

%% Define unique fragments positions
nonuniq_frag_pos = cell(length(target_name), 1);
nonuniq_frag_pos_over = cell(length(target_name), 1);
nonuniq_frag_seq = cell(length(target_name), 1);
nonuniq_frag_seq_over = cell(length(target_name), 1);

%% Split histogram contiguous regions into fragments
for i = 1 : length(target_name)
    m = 0;  % Variable numbering unique regions in given target
    % Index of each segment within the non-unique regions
    n = zeros(1, size(nonuniq_reg{i}, 2));
    %rep_hits{i} = cell(length(hist_map{i}), 1);
    for j = 1 : size(hist_pos{i}, 2)
        % Length of given histogram stretch
        hist_map_len = hist_pos{i}(2, j) - hist_pos{i}(1, j) + 1;
        % Histogram position in non-unique region
        ind = find(sum(sign(mean(hist_pos{i}(:, j)) - ...
            nonuniq_reg{i})) == 0);
        % Process long enough histogram stretches only
        if hist_map_len >= min_rep_pcr_size
            % Number of fragments to cover given histogram stretch
            nr_pcrs = ceil(hist_map_len / max_rep_pcr_size);
            for k = 1 : nr_pcrs
                m = m + 1;
                % Starting position of PCR fragment
                nonuniq_frag_pos{i}(1, m) = ...
                    round(hist_map_len * (k - 1) / nr_pcrs) + ...
                    hist_pos{i}(1, j);
                % Ending position of PCR fragment
                nonuniq_frag_pos{i}(2, m) = ...
                    round(hist_map_len * k / nr_pcrs) - ...
                    1 + hist_pos{i}(1, j);
                % Starting position of PCR fragment with overlap; it is
                % trimmed so that it does not span further than the unique
                % region
                nonuniq_frag_pos_over{i}(1, m) = ...
                    max([hist_pos{i}(1, j), nonuniq_frag_pos{i}(1, m) - ...
                    pcr_overlap]);
                % Ending position of PCR fragment with overlap; it is
                % trimmed so that it does not span further than the unique
                % region
                nonuniq_frag_pos_over{i}(2, m) = ...
                    min([hist_pos{i}(2, j), nonuniq_frag_pos{i}(2, m) + ...
                    pcr_overlap]);
                % Sequence of the PCR fragment
                nonuniq_frag_seq{i}{m} = ...
                    target_seq{i}((nonuniq_frag_pos{i}(1, m) : ...
                    nonuniq_frag_pos{i}(2, m)));
                % Sequence of the PCR fragment with overlap in lower case
                % letters
                nonuniq_frag_seq_over{i}{m} = ...
                    lower(target_seq{i}((nonuniq_frag_pos_over{i}(1, m) ...
                    : nonuniq_frag_pos_over{i}(2, m))));
                % Capitalize letters of the non-overlapping part of the PCR
                % fragment
                in = (nonuniq_frag_pos{i}(1, m) - ...
                    nonuniq_frag_pos_over{i}(1, m)) + ...
                    (1 : (nonuniq_frag_pos{i}(2, m) - ...
                    nonuniq_frag_pos{i}(1, m) + 1));
                % Capitalize letters of the non-overlapping part of the PCR
                % fragment
                nonuniq_frag_seq_over{i}{m}(in) = ...
                    upper(nonuniq_frag_seq_over{i}{m}(in));
                % Give fragment a name
                %nonuniq_frag_name{i}{m} = ['N' num2str(k) alph(j)];
                n(ind) = n(ind) + 1;
                nonuniq_frag_name{i}{m} = ['N' num2str(ind) alph(n(ind))];
            end
        end
    end
end

%% Delete fragments outside of BACs
nonuniq_frag_pos{1}(:, 8 : 10) = [];
nonuniq_frag_pos_over{1}(:, 8 : 10) = [];
nonuniq_frag_seq{1}(:, 8 : 10) = [];
nonuniq_frag_seq_over{1}(:, 8 : 10) = [];
nonuniq_frag_name{1}(:, 8 : 10) = [];


%% Save output
save([folder 'mats/fp' num2str(process) '.mat'], 'nonuniq_frag_pos', ...
    'nonuniq_frag_pos_over', 'nonuniq_frag_seq', ...
    'nonuniq_frag_seq_over', 'nonuniq_frag_name');
write_log(process, 'Generated PCR fragments for non-unique sequences');

%% Process finished
write_log(process);