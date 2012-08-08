function graph_par(process)
global last_process
global folder
global xaxis_sigdig
global yaxis_sigdig
global im_pos
global dpi_res
global gformat
global hist_color
global yel_color
global shade_yel_color
global white_color
global shade_white_color
global genes_color
global unseq_color
global nonunseq_color
global mistarget_color
global title_fsize
global axis_fsize
global axis_num_fsize
global genes_fsize
global segment_fsize
global region_fsize
global fontname
global genes_ypos
global unseg_ypos
global nonunseg_ypos
global unreg_ypos
global nonunreg_ypos
global rev_xaxis
global en_hist
global en_area
global en_genes
global en_unseg
global en_nonunseg
global en_unreg
global en_nonunreg
global en_mistarget

%% Load saved target sequences
if last_process >= process
    load([folder 'mats/fp' num2str(process) '.mat']);
    write_log(process, 'Analysis parameters loaded');
    write_log(process);
    return
end

%% Import target sequences from fasta file
xaxis_sigdig = str2double(read_config('xaxis_sig_digits', process));
yaxis_sigdig = str2double(read_config('yaxis_sig_digits', process));
im_pos = ...
    str2num(read_config('image_position', process)); %#ok<ST2NM>
dpi_res = str2double(read_config('dpi_res', process));
gformat = lower(read_config('graph_format', process));
hist_color = str2num(read_config('hist_color', process)); %#ok<ST2NM>
yel_color = str2num(read_config('nonunique_color', process)); %#ok<ST2NM>
shade_yel_color = ...
    str2num(read_config('shade_nonunique_color', process)); %#ok<ST2NM>
white_color = str2num(read_config('unique_color', process)); %#ok<ST2NM>
shade_white_color = ...
    str2num(read_config('shade_unique_color', process)); %#ok<ST2NM>
genes_color = lower(read_config('genes_color', process));
unseq_color = ...
    str2num(read_config('unique_segment_color', process)); %#ok<ST2NM>
nonunseq_color = read_config('non-unique_segment_color', process);
mistarget_color = read_config('mistarget_color', process);
title_fsize = str2double(read_config('title_size', process));
axis_fsize = str2double(read_config('axis_size', process));
axis_num_fsize = str2double(read_config('axis_num_size', process));
genes_fsize = str2double(read_config('genes_size', process));
segment_fsize = str2double(read_config('segment_size', process));
region_fsize = str2double(read_config('region_size', process));
fontname = lower(read_config('font_name', process));

genes_ypos = str2double(read_config('genes_ypos', process));
unseg_ypos = str2double(read_config('unseg_ypos', process));
nonunseg_ypos = str2num(read_config('nonunseg_ypos', process)); %#ok<ST2NM>
unreg_ypos = str2double(read_config('unreg_ypos', process));
nonunreg_ypos = str2double(read_config('nonunreg_ypos', process));
rev_xaxis = str2double(read_config('rev_xaxis', process));
en_hist = str2double(read_config('plot_hist', process));
en_area = str2double(read_config('plot_area', process));
en_genes = str2double(read_config('plot_genes', process));
en_unseg = str2double(read_config('plot_unseg', process));
en_nonunseg = str2double(read_config('plot_nonunseg', process));
en_unreg = str2double(read_config('plot_unreg', process));
en_nonunreg = str2double(read_config('plot_nonunreg', process));
en_mistarget = str2double(read_config('plot_mistarget', process));


%% Check values
xaxis_sigdig = check_value(process, xaxis_sigdig, 'xaxis_sigdig', ...
    [5, 5, 5], [1, 10]);
yaxis_sigdig = check_value(process, yaxis_sigdig, 'yaxis_sigdig', ...
    [1, 1, 1], [1, 10]);
if length(im_pos) ~= 4
    im_pos = [0.05, 0.1, 0.4, 0.8];
	error_msg(process, 1, {'Using im_pos of [0.05, 0.1, 0.4, 0.8]', ...
        'Could not acquire im_pos'});
else
    def_val = [0, 0, 0; 0, 0, 0; 1, 1, 1; 1, 1, 1];
    for i = 1 : length(im_pos)
        im_pos(i) = check_value(process, im_pos(i), ...
            ['impos(' num2str(i) ')'], def_val(i, :), [0, 1]);
    end
end
dpi_res = check_value(process, dpi_res, 'dpi_res', ...
    [150, 50, 1200], [50, 1200]);
gformats = {'-dbmpmono', '-dbmp16m', '-dbmp256', '-dbmp', '-dmeta', ...
    '-deps', '-depsc', '-deps2', '-depsc2', '-dhdf', '-dill', '-djpeg', ...
    '-dpbm', '-dpbmraw', '-dpcxmono', '-dpcx24b', '-dpcx256', ...
    '-dpcx16', '-dpdf', '-dpgm', '-dpgmraw', '-dpng', '-dppm', ...
    '-dppmraw', '-dtiff', '-dtiffn', '-tiff'};
if ~cellfun(@(x) strcmp(x, gformat), gformats)
	error_msg(process, 1, {'Using gformat ''-dpng''', ...
        'Could not acquire gformat'});
    gformat = '-dpng';
end
hist_color = check_color(process, hist_color, 'hist_color', [0; 0.7; 0]);
yel_color = check_color(process, yel_color, 'nonunique_color', [1; 1; 0]);
shade_yel_color = check_color(process, shade_yel_color, ...
    'shade_nonunique_color', [0.75;0.75;0]);
white_color = check_color(process, white_color, 'unique_color', [1; 1; 1]);
shade_white_color = check_color(process, shade_white_color, ...
    'shade_unique_color', [0.75; 0.75; 0.75]);
colors = {'yellow', 'magenta', 'cyan', 'red', 'green', 'blue', 'white', ...
    'black'};
if ~cellfun(@(x) strcmp(x, genes_color), colors)
	error_msg(process, 1, {'Using genes_color ''red''', ...
        'Could not acquire genes_color'});
    genes_color = 'red';
end
unseq_color = check_color(process, unseq_color, 'unique_segment_color', ...
    [1; 0; 1]);
if ~cellfun(@(x) strcmp(x, nonunseq_color), colors)
	error_msg(process, 1, {'Using nonunique_segment_color ''blue''', ...
        'Could not acquire nonunique_segment_color'});
    nonunseq_color = 'blue';
end
if ~cellfun(@(x) strcmp(x, mistarget_color), colors)
	error_msg(process, 1, {'Using mistarget_color ''black''', ...
        'Could not acquire mistarget_color'});
    mistarget_color = 'black';
end

title_fsize = check_value(process, title_fsize, 'title_fsize', ...
    [14, 6, 20], [2, 30]);
axis_fsize = check_value(process, axis_fsize, 'axis_fsize', ...
    [12, 6, 20], [2, 30]);
axis_num_fsize = check_value(process, axis_num_fsize, 'axis_num_fsize', ...
    [10, 6, 20], [2, 30]);
genes_fsize = check_value(process, genes_fsize, 'genes_fsize', ...
    [8, 5, 20], [2, 30]);
segment_fsize = check_value(process, segment_fsize, 'segment_fsize', ...
    [6, 5, 20], [2, 30]);
region_fsize = check_value(process, region_fsize, 'region_fsize', ...
    [8, 5, 20], [2, 30]);
fontnames = {'arial'};
if ~cellfun(@(x) strcmp(x, fontname), fontnames)
	error_msg(process, 1, {'Using fontname ''arial''', ...
        'Could not acquire fontname'});
    fontname = 'arial';
end

genes_ypos = check_value(process, genes_ypos, 'genes_ypos', ...
    [0.1, 0, 1], [0, 1]);
unseg_ypos = check_value(process, unseg_ypos, 'unseg_ypos', ...
    [0.1, 0, 1], [0, 1]);
if length(nonunseg_ypos) ~= 2
    nonunseg_ypos = [0.4, 0.4];
	error_msg(process, 1, {'Using nonunseg_ypos of [0.4, 0.4]', ...
        'Could not acquire nonunseg_ypos'});
else
    def_val = [0.4, 0.4, 0.4; 0.4, 0.4, 0.4];
    for i = 1 : length(nonunseg_ypos)
        nonunseg_ypos(i) = check_value(process, nonunseg_ypos(i), ...
            ['impos(' num2str(i) ')'], def_val(i, :), [0, 1]);
    end
end
unreg_ypos = check_value(process, unreg_ypos, 'unreg_ypos', ...
    [0.1, 0, 1], [0, 1]);
nonunreg_ypos = check_value(process, nonunreg_ypos, 'nonunreg_ypos', ...
    [0.1, 0, 1], [0, 1]);

rev_xaxis = check_logical(process, rev_xaxis, 'rev_xaxis', false);
en_hist = check_logical(process, en_hist, 'plot_hist', true);
en_hist = check_logical(process, en_hist, 'en_hist', true);
en_area = check_logical(process, en_area, 'plot_area', true);
en_genes = check_logical(process, en_genes, 'plot_genes', true);
en_unseg = check_logical(process, en_unseg, 'plot_unseg', true);
en_nonunseg = check_logical(process, en_nonunseg, 'plot_nonunseg', true);
en_unreg = check_logical(process, en_unreg, 'plot_unreg', true);
en_nonunreg = check_logical(process, en_nonunreg, 'plot_nonunreg', true);
en_mistarget = check_logical(process, en_mistarget, 'plot_mistarget', true);

%% Save output
save([folder 'mats/fp' num2str(process) '.mat'], 'xaxis_sigdig', ...
    'yaxis_sigdig', 'im_pos', 'dpi_res', 'gformat', 'hist_color', ...
    'yel_color', 'shade_yel_color', 'white_color', 'shade_white_color', ...
    'genes_color', 'unseq_color', 'nonunseq_color', 'mistarget_color', ...
    'title_fsize', 'axis_fsize', 'axis_num_fsize', 'genes_fsize', ...
    'segment_fsize', 'region_fsize', 'fontname', 'genes_ypos', ...
    'unseg_ypos', 'nonunseg_ypos', 'unreg_ypos', 'nonunreg_ypos', ...
    'rev_xaxis', 'en_hist', 'en_area', 'en_genes', 'en_unseg', ...
    'en_nonunseg', 'en_unreg', 'en_nonunreg', 'en_mistarget');
write_log(process, 'Imported analysis parameters');

%% Process finished
write_log(process);