function gen_pix(process)
global last_process
global aln_pos
global target_name
global nonuniq_vec
global hist_color
global en_hist
global dpi_res
global gformat
global sel_nonun_frags
global log_file

%% Check if graphics is ready
if last_process >= process
    write_log(process, 'Graphics generated');
    write_log(process);
    return
end

%% Open a new figure window and set its size and position
close all;
fig = zeros(1, length(target_name));
ax = zeros(1, length(target_name));
for i = 1 : length(target_name)
    for k = 1 : length(sel_nonun_frags{i}) + 1
    	write_log(process, ['Generate Graphics for ''' target_name{i} '''']);
        % Generate the figure
        [fig(i), ax(i)] = new_figure;
        hold on
        % Generate ylims
        ylims = ylimgen(i);
        % Find coordinates where any transition between unique/non-unique
        % region or mapped area appears and plot them
        plot_areas(ylims, i);
        % Plot histogram of number of repeats confined to each non-unique
        % region
        if en_hist == true
            area((1 : length(nonuniq_vec{i})) + aln_pos(i, 1), ...
            nonuniq_vec{i}, 'FaceColor', hist_color, 'LineStyle', 'none');
        end
        % Plot Genes Names
        plot_genes(aln_pos(i, 1 : 2), ylims);
        % Plot Unique and Non-Unique Regions
        plot_regions(ylims, i);
        % Plot Fragments
        plot_uniq_frags(ylims, i)
        plot_nonuniq_frags_new(ylims, i, k)
        % Label Graph
        label_graph(ylims, i, ax);
        if k <= length(sel_nonun_frags{i})
            t = num2str(k);
        else
            t = 'all';
        end
        print(gcf, ['-r' num2str(dpi_res)], gformat, ...
            [log_file(1 : (end - 4)) '_' num2str(i) '_' t]);
        close all;
    end
end

%% Save output
write_log(process, 'Generated graphic output');

%% Process finished
write_log(process);