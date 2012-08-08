function plot_areas(ylims, i)
global aln_pos
global uniq_map
global nonuniq_map
global nonuniq_reg
global en_area
global shade_yel_color
global shade_white_color
global yel_color
global white_color

if en_area == false
    return
end

%% Find coordinates where any transition between unique/non-unique region
%% or mapped area appears
separ = sort(horzcat(uniq_map{i}(1, :), uniq_map{i}(2, :), ...
	nonuniq_map{i}(1, :), nonuniq_map{i}(2, :), ...
    nonuniq_reg{i}(1, :), nonuniq_reg{i}(2, :)));

%% Check each region separated by the coordinates and deide which color it
%% should be plotted with
for j = 1 : length(separ) - 1
    % Check if the separators are not next to each other
    if separ(j + 1) - separ(j) > 1
        % mid is a number in the middle of the area confined by the
        % separators
        mid = mean(separ(j : j + 1));
        % color index
        colin = 0;
        % if mid is located in a unique map select color index 1
        if find(sum(sign(uniq_map{i} - mid)) == 0)
            colin = 1;
        end
        % if mid is located in a non-unique map increase color index by 2
        if find(sum(sign(nonuniq_map{i} - mid)) == 0)
            colin = colin + 2;
        end
        % if mid is located in a non-unique region increase color index by
        % 4
        if find(sum(sign(nonuniq_reg{i} - mid)) == 0)
            colin = colin + 4;
        end
        % Generate coordinates for the fill area
        X = [separ(j), separ(j + 1), separ(j + 1), separ(j)] + aln_pos(i, 1);
        Y = [0, 0, ylims(2), ylims(2)];
        if colin == 0
            C = [0, 0, 0];
        elseif colin == 1
            C = white_color;            % white
        elseif colin == 2
            C = yel_color;              % yellow
        elseif colin == 3
            C = [0, 0, 0];
        elseif colin == 4
            C = [0, 0, 0];
        elseif colin == 5
            C = shade_white_color;      % grey
        elseif colin == 6
	        C = shade_yel_color;        % yellow/grey
        elseif colin == 7
	        C = [0, 0, 0];
        end
        % Generate colored area
        fill(X, Y, C, 'EdgeColor', 'none')
    end
end