# PieParty v1.3 - Visualizing cells from scRNA-seq data as pie charts
### PieParty is a visualization tool allowing to represent every cell in single-cell sequencing plots (UMAP, tSNE, ect.) as a pie chart. Each slice in a pie chart represents the expression of a single gene. Custom gene lists and coloring can be applied. 

<img src="https://github.com/harbourlab/PieParty/blob/master/testis.png" width="300">


_Requirements: Python 3. Please do not use Python 2.<br>
Python packages required: matplotlib

Usage: 
<pre>
python PieParty.py -g expression_file.csv -c cell_coordinates.csv -l genelist1.csv genelist2.csv
</pre>

### Options (required):

-g _gene expression file.csv_<br>
For the expression matrix, you will need a matrix of RNA counts with the cell IDs as columns and           gene names as rows.    The expression matrix column names should match the cell IDs in the UMAP/tSNE coordinate table. In Seurat, this count matrix can be generated from accessing the data under Seurat_Object@assays$RNA@data.

-c _cell coordinates.csv_ <br>
For the UMAP/tSNE coordinate table, you will need a column for each reduction (UMAP or tSNE) dimension (column 1 = UMAP_1, column 2 = UMAP_2). Row names should be the cell IDs, matching the expression matrix. In Seurat, this data can be accessed under Seurat_Object@reductions$umap@cell.embeddings.

-l _genelist.csv_ <br>
lists of genes that should be plotted (csv). One list is minimum. If more lists are provided, separate coloring can be applied. Gene lists need to be csv files with one column, and one gene per row.



### Options (not required):

-o _output_file_name.png_ <br> Default = "output.png"

-color _(hex colors or matplotlib colormaps)_ <br>
If multiple gene lists are used same amount of colors should be provided. e.g. "-color autumn" for one gene list, or "-color @FF00FF @FFFF00" for two gene lists. The latter example are hex codes. Colormap names can be found here https://matplotlib.org/tutorials/colors/colormaps.html. Default is "viridis". In case a colormap is chosen, PieParty will auto-assign colors in order according to the genes in the gene lists provided.

-p proportionalize (True or False) <br>
If two or more gene lists are provided, PieParty can normalize for difference in gene amounts in the lists. This is useful in many cases, as lists of e.g. 4 macrophage markers, and 102 keratinocyte markers will produce pie charts that are predomonatelly filled with slices from the keratinocyte markers, although overall each individual marker may be expressed less as the 4 macrophage markers. Default is True.

-ct percentage cutoff (float) <br>
Cutoff percentage of expression in a pie chart a gene has to meet to be included. Default is 1%, meaning if a gene is not expressed at least 1% of the summed expression of all other genes in the pie it is excluded.

-ce expression cutoff (float) <br>
If desired, a expression cutoff can be applied. Default is 0. 

-lc lighten colors (True or False) <br>
This setting will lighten the colors of the pie slices according to their expression value. Default is True.

-gc lighten colors based on global (True or False) <br>
If this is set to True, PieParty will use color intesities reflecting the expression compared to ALL other genes in the single cell dataset. This might lead to very light colored plots in many cases. Default is False.

-pr plot resolution (int) <br>
Default is 13000, which results in high-resolution plots and reccomended. 2400 is still good enough in most cases and decreases computation time if wished.

