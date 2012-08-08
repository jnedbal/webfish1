function report_nonun_frags(process)
global last_process
global log_file
global target_name
global nonuniq_frag_name
global nonuniq_frag_pos
global rep_hits
global nonuniq_frag_seq
global nonuniq_frag_seq_over
global aln_pos
global sel_nonun_frags
global fatal_error
global enzymes
global folder

%% Check if fragments are produced
if last_process >= process
    write_log(process, 'Report about non-unique segments finished');
    write_log(process);
    return
end

%% Open a file to write the fragment sequences
[fid, w] = fopen([log_file(1 : end - 3) 'nonunique_segs.txt'], 'w');
if fid == -1
    fatal_error = 1;
    fatal_msg(process, ...
        {'Failed writing file with non-unique fragments', w});
    return
end

%% Print the report
for i = 1 : length(target_name)
    fprintf(fid, [repmat('*', 1, length(target_name{i}) + 8) '\n']);
    fprintf(fid, ['*   ' target_name{i} '   *\n']);
    fprintf(fid, [repmat('*', 1, length(target_name{i}) + 8) '\n\n']);
    if isempty(nonuniq_frag_name)
        continue
    end
    AB = false(size(nonuniq_frag_seq{i}, 2), aln_pos(i, 2) - aln_pos(i, 1));
    for m = 1 : length(nonuniq_frag_seq{i})
        fprintf(fid, ['Target sequence: ' target_name{i} '\n']);
        fprintf(fid, ['Segment name: ' nonuniq_frag_name{i}{m} '\n']);
        fprintf(fid, ['Segment from: ' ...
            num2str(nonuniq_frag_pos{i}(1, m)) ' bp\n']);
        fprintf(fid, ['Segment to: ' ...
            num2str(nonuniq_frag_pos{i}(2, m)) ' bp\n']);
        fprintf(fid, ['Segment length: ' ...
            num2str(nonuniq_frag_pos{i}(2, m) - ...
            nonuniq_frag_pos{i}(1, m) + 1) ' bp\n']);
        aln_bin = false(1, aln_pos(i, 2) - aln_pos(i, 1));
        aln_bin(nonuniq_frag_pos{i}(1, m) : nonuniq_frag_pos{i}(2, m)) ...
            = true;
        if ~isempty(rep_hits{i}{m})
            for k = 1 : size(rep_hits{i}{m}, 2)
                s = sort(rep_hits{i}{m}(1 : 2, k));
                aln_bin(s(1) : s(2)) = 1;
            end
        end
        aln_len = length(find(aln_bin));
        AB(m, :) = aln_bin;
        fprintf(fid, ['Total length covered by the segment: ' ...
            num2str(aln_len) ' bp\n']);
        for n = 1 : length(enzymes)
            [e sites] = ...
                rebasecuts(nonuniq_frag_seq_over{i}{m}, enzymes{n});
            fprintf(fid, ['Number of ' enzymes{n} ' sites: ' ...
                num2str(length(sites)) '\n']);
        end
        fprintf(fid, ['Segment sequence: ' nonuniq_frag_seq{i}{m} '\n']);
        fprintf(fid, ['Segment sequence with overhang: ' ...
            nonuniq_frag_seq_over{i}{m} '\n\n']);
    end
    for n = 1 : min(5, length(nonuniq_frag_seq{i}))
        write_log(process, ['Checking combinations of ' num2str(n) ...
            ' non-unique segments']);
        combs = nchoosek(1 : length(nonuniq_frag_seq{i}), n);
        aln_len = zeros(size(AB, 1), 1);
        for p = 1 : size(combs, 1)
            aln_bin = false(1, size(AB, 2));
            for q = 1 : size(combs, 2)
                aln_bin = aln_bin | AB(combs(p, q), :);
            end
            aln_len(p) = length(find(aln_bin));
        end
        [v, max_len] = max(aln_len);
        fprintf(fid, ['Total length covered by ' num2str(n) ...
            ' segments: ' num2str(aln_len(max_len)) '\n']);
        for q = 1 : size(combs, 2)
            fprintf(fid, ['Segment name: ' ...
                nonuniq_frag_name{i}{combs(max_len, q)} '\n']);
            sel_nonun_frags{i}{n}(q) = combs(max_len, q);
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);

%% Save output
save([folder 'mats/fp' num2str(process) '.mat'], 'sel_nonun_frags');
write_log(process, 'Wrote report about non-unique segments');

%% Process finished
write_log(process);