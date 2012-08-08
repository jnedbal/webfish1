function plot_genes(aln, ylims)
global genes_name
global genes_pos
global genes_ypos
global genes_color
global fontname
global genes_fsize
global en_genes

if en_genes == false
    return
end

% Plot Genes Names
for j = 1 : length(genes_name)
    c = genes_pos(j, 1 : 2);
    gene_mid = mean(c);
    if gene_mid >= aln(1) && gene_mid <= aln(2)
        plot(c, ylims(2) * [genes_ypos, genes_ypos], '-', ...
            'LineWidth', 2, 'Color', genes_color);
        yc = ylims(2) * (genes_ypos + [- 0.005, 0.005]);
    	plot([c(1), c(1)], yc, '-', 'LineWidth', 0.1, ...
            'Color', genes_color);
    	plot([c(2), c(2)], yc, '-', 'LineWidth', 0.1, ...
            'Color', genes_color);
    	text_label = ' \leftarrow';
    	text(gene_mid, ylims(2) * genes_ypos, text_label, ...
            'HorizontalAlignment', 'left', 'Rotation', 90, ...
            'VerticalAlignment', 'middle', 'Color', genes_color, ...
            'FontName', fontname, 'FontSize', genes_fsize);
    	text_label = ['    ' genes_name{j}];
    	text(gene_mid, ylims(2) * genes_ypos, text_label, ...
            'HorizontalAlignment', 'left', 'Rotation', 90, ...
            'VerticalAlignment', 'middle', 'Color', genes_color, ...
            'FontName', fontname, 'FontSize', genes_fsize);
    end
end