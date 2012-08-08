function plot_uniq_frags(ylims, i)
global uniq_frag_pos
global uniq_frag_name
global aln_pos
global en_unseg
global unseg_ypos
global unseq_color
global fontname
global segment_fsize

if en_unseg == false
    return
end
% Plot unique fragments
for j = 1 : size(uniq_frag_pos{i}, 2)
    c = [uniq_frag_pos{i}(1, j), uniq_frag_pos{i}(2, j)] + aln_pos(i, 1);
    frag_mid = mean(c);
	plot(c, ylims(2) * [unseg_ypos, unseg_ypos], '-', 'LineWidth', 2, ...
        'Color', unseq_color)
    yc = ylims(2) * (unseg_ypos + [- 0.005, 0.005]);
  	plot([c(1), c(1)], yc, '-', 'LineWidth', 0.1, 'Color', unseq_color);
	plot([c(2), c(2)], yc, '-', 'LineWidth', 0.1, 'Color', unseq_color);
    [v, e] = num2eng(c(2) - c(1) + 1, 2);
    text_label = [uniq_frag_name{i}{j} ' (' v ' ' e 'bp) \rightarrow'];
	text(frag_mid, ylims(2) * unseg_ypos, text_label, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'Color', unseq_color, 'Rotation', 90, 'FontName', fontname, ...
        'FontSize', segment_fsize);
end