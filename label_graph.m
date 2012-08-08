function label_graph(ylims, i, ax)
global aln_chr
global aln_pos
global organism
global target_name
global xaxis_sigdig
global en_hist
global hist_color
global axis_num_fsize
global axis_fsize
global title_fsize
global xdir_xaxis
global fontname

xlims = aln_pos(i, 1 : 2);
xlim_lab = cell(1, 2);
for j = 1 : 2
	[v, c] = num2eng(aln_pos(i, j), xaxis_sigdig);
    xlim_lab{j} = horzcat(v, ' ', c, 'bp');
end
org = horzcat(upper(organism(1)), lower(organism(2 : end)));
in = strfind(org, '_');
org(in) = ' ';
org(in + 1) = upper(org(in + 1));
xaxis = ['Position on ' aln_chr{i} ' of ' org];
if en_hist == true
	yaxis = 'Number of Similar Hits';
    yticks = ylims;
else
	yaxis = '';
    yticks = [];
end
tit = [' Similar and Unique Regions in Target Region ''' target_name{i} ...
    ''' '];
if xdir_xaxis == false
    xdir = 'normal';
else
    xdir = 'reverse';
end
set(ax(i), 'XLim', xlims, 'YLim', ylims, 'XTick', xlims, 'YTick', ...
    yticks, 'XTickLabel', xlim_lab, 'Box', 'on', 'XDir', xdir, ...
    'FontSize', axis_num_fsize, 'FontName', fontname);
xlabel(xaxis, 'FontSize', axis_fsize, 'FontName', fontname);
ylabel(yaxis, 'FontSize', axis_fsize, 'FontName', fontname, ...
    'Color', hist_color);
title(tit, 'FontSize', title_fsize, 'FontName', fontname);