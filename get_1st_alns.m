function [aln_pos, aln_chr, aln_len, tar_pos] = ...
    get_1st_alns(process, target_name, target_seq, query_from, ...
    query_to, hit_from, hit_to, chromosome, chr)

aln_chr = cell(length(target_name), 1);
tar_pos = zeros(length(target_name), 2);
aln_pos = zeros(length(target_name), 2);
aln_len = zeros(length(target_name), 1);
for j = 1 : length(target_name)
    tar_pos(j, 1 : 2) = sort([query_from{j}(1), query_to{j}(1)]);
    aln_pos(j, 1 : 2) = sort([hit_from{j}(1), hit_to{j}(1)]);
    aln_len(j) = aln_pos(j, 2) - aln_pos(j, 1) + 1;
    aln_chr{j} = horzcat(upper(chr{chromosome{j}(1)}(1)), ...
        lower(chr{chromosome{j}(1)}(2:end)));
    aln_chr{j} = regexprep(aln_chr{j}, '_', ' ');
    if length(target_seq{j}) ~= tar_pos(j, 2) - tar_pos(j, 1) + 1 || ...
            length(target_seq{j}) ~= aln_len(j)
        error_msg(process, 1, {['Length of target sequence "' ...
            regexprep(target_name{j}, '\\', '\\\') ...
            '" and first alignment do not match']});
    end
    write_log(process, ['Localized target sequence "' ...
        regexprep(target_name{j}, '\\', '\\\') '" on "' aln_chr{j} '"']);
end