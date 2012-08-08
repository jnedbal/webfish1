function [fig, ax] = new_figure
global im_pos

%% Open a new figure window and set its size and position
fig = figure('Visible', 'off');
units = get(gcf, 'units');
set(gcf, 'units', 'normalized', 'outerposition', im_pos);
set(gcf, 'units', units);
ax = axes;