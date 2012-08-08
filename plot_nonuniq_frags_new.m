function plot_nonuniq_frags_new(ylims, i, k)
global nonuniq_frag_pos
global nonuniq_frag_name
global aln_pos
global rep_hits
global en_nonunseg
global en_mistarget
global nonunseg_ypos
global nonunseq_color
global mistarget_color
global fontname
global segment_fsize
global sel_nonun_frags

if en_nonunseg == false
    return
end

% Plot non-unique fragments
ycoef = nonunseg_ypos(1);
ycoeft = nonunseg_ypos(1);
y_coef_step = ...
    min([abs(diff(nonunseg_ypos)) / (size(nonuniq_frag_pos{i}, 2) - 1), ...
    0.02]);

if k <= length(sel_nonun_frags{i})
    used_frag_pos = nonuniq_frag_pos{i}(:, sel_nonun_frags{i}{k});
    used_frag_name = nonuniq_frag_name{i}(sel_nonun_frags{i}{k});
    used_rep_hits = rep_hits{i}(sel_nonun_frags{i}{k});
else
    used_frag_pos = nonuniq_frag_pos{i};
    used_frag_name = nonuniq_frag_name{i};
    used_rep_hits = rep_hits{i};
end
for j = 1 : size(used_frag_pos, 2)
    if en_mistarget == true
        aln_bin = false(1, aln_pos(i, 2) - aln_pos(i, 1));
        aln_bin(used_frag_pos(1, j) : ...
            used_frag_pos(2, j)) = true;
        if ~isempty(used_rep_hits{j})
            for k = 1 : size(used_rep_hits{j}, 2)
                s = sort(used_rep_hits{j}(1 : 2, k));
                aln_bin(s(1) : s(2)) = 1;
            end
        end
        aln_len = length(find(aln_bin));
    end
    c = used_frag_pos(1 : 2, j) + aln_pos(i, 1);
    frag_mid = mean(c);
	plot(c, ylims(2) * [ycoef, ycoef], '-', 'LineWidth', 2, ...
        'Color', nonunseq_color)
    yc = ylims(2) * (ycoef + [- 0.005, 0.005]);
  	plot([c(1), c(1)], yc, 'b-', 'LineWidth', 0.1, ...
        'Color', nonunseq_color);
	plot([c(2), c(2)], yc, 'b-', 'LineWidth', 0.1, ...
        'Color', nonunseq_color);
    [v, e] = num2eng(c(2) - c(1) + 1, 2);
    [v1, e1] = num2eng(aln_len, 2);
    if en_mistarget == true
        text_label = ['\color{' nonunseq_color '}' ...
            used_frag_name{j} ' (' v ' ' e 'bp/\color{' ...
            mistarget_color '}' v1 ' ' e1 'bp\color{' nonunseq_color ...
            '}) \rightarrow'];
    else
        text_label = ['\color{' nonunseq_color '}' ...
            used_frag_name{j} ' (' v ' ' e 'bp) \rightarrow'];
    end
	text(frag_mid, ylims(2) * ycoeft, text_label, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'Rotation', 90, 'FontName', fontname, 'FontSize', segment_fsize);
    if en_mistarget == true
        if ~isempty(used_rep_hits{j})
            for k = 1 : size(used_rep_hits{j}, 2)
                plot(used_rep_hits{j}(1 : 2, k) + aln_pos(i, 1), ...
                    ylims(2) * [ycoef, ycoef], '-', 'LineWidth', 1, ...
                    'Color', mistarget_color)
            end
        end
    end
    ycoef = ycoef + y_coef_step;
end