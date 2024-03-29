# PieParty v1.7.3 - Visualizing cells from scRNA-seq data as pie charts
### PieParty is a visualization tool allowing to represent every cell in single-cell sequencing plots (UMAP, tSNE, ect.) as a pie chart. Each slice in a pie chart represents the expression of a single gene. Custom gene lists and coloring can be applied. 

<img src="https://github.com/harbourlab/PieParty/blob/master/testis.png" width="300">

<img src="https://github.com/harbourlab/PieParty/blob/master/Screen Shot 2020-09-18 at 8.03.00 PM.png" width="800">

<img src="https://github.com/harbourlab/PieParty/blob/master/Screen Shot 2020-09-18 at 8.03.17 PM.png" width="800">

<img src="https://github.com/harbourlab/PieParty/blob/master/2_t6_giant_pies.png" width="300">



Please cite:
Kurtenbach, S., Dollar, J. J., Cruz, A. M., Durante, M. A., Decatur, C. L., and Harbour, J. W.; PieParty: visualizing cells from scRNA-seq data as pie charts; Life Sci Alliance; 5 4 2021



_Requirements: Python 3. Please do not use Python 2.<br>
Python packages required: matplotlib_

Usage: 
<pre>
python PieParty.py -g expression_file.csv -c cell_coordinates.csv -l genelist1.csv genelist2.csv
</pre>

PieParty plots depict every cell in single cell plots such as UMAP and tSNE as pie charts, where every slice of the pie chart corresponds to the relative gene expression of one gene. A set of genes can be defined to be plotted, and custom coloring can be defined. Different sets of genes can be used, e.g. one list with macrophage markers, and one list with cancer cell markers. Every list can be assigned a unique color or color palette, e.g. macrophage markers in yellow, and cancer cells in red.


This will generate the following outputs:
1) The .png picture containing the PieParty figure
2) "labels.svg" which lists all genes and which colors were assigned

Note:
While running the PieParty script the command line will issue a warning, as the generated PNG file is very large. The warning is normal and can be ignored. "DecompressionBombWarning: Image size (169000000 pixels) exceeds limit of
89478485 pixels, could be decompression bomb DOS attack."

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
If two or more gene lists are provided, PieParty can normalize for difference in gene amounts in the lists. This is useful in many cases, as lists of e.g. 4 macrophage markers, and 102 keratinocyte markers will produce pie charts that are predomonatelly filled with slices from the keratinocyte markers, although overall each individual marker may be expressed less as the 4 macrophage markers. The normalized expression value is calculated by multiplication of the expression value with the number of total genes in all lists, divided by the number of genes in the list of the respective gene (N_(gene_list)). Default is True.

-ct percentage cutoff (float) <br>
Cutoff percentage of expression in a pie chart a gene has to meet to be included. Default is 1%, meaning if a gene is not expressed at least 1% of the summed expression of all other genes in the pie it is excluded.

-ce expression cutoff (float) <br>
If desired, a expression cutoff can be applied. Default is 0. 

-lc lighten colors (True or False) <br>
This setting will lighten the colors of the pie slices according to their expression value. Default is True.

-gc lighten colors based on global (True or False) <br>
If multiple gene lists are used as inputs, it might be desired that PieParty does not use the same scale to lighten the colors of pie slices. For instance, if gene list A has high expressing genes, and gene list B has only low expressing genes, the colors will be very faint for the second list in general. In this case -gc can be set to False, to increase color saturation. PieParty will then use the maximum gene expression value for each list to determinte color saturation. However, be aware that if this is set to False a direct comparison of gene expression is not possible. Default is True.

-pr plot resolution (int) <br>
Default is 13000, which results in high-resolution plots and reccomended. 2400 is still good enough in most cases and decreases computation time if wished.

-gp plot giant pies (True, False) <br>
Set this to "True" to plot one giant pie per cell cluster. The size of the pie chart will correlate with the number of cells in the cluster. 
Plotting giant pies requires a csv file defining which cell (1st column) is in which cluster (2nd column). The cluster file has to be defined with "-cf". I reccoment trying this without the lighten colors function "-lc False", as depending on your dataset the giant pies will be light in color.
The coordinates file can be exported as a csv file in R with “write.csv(SeuratObject@reductions$umap@cell.embeddings, file = "SeuratObject_UMAPcoordinates.csv”)”

-ps pie scaling (float) <br>
Scaling factor by which the pie sizes are multiplied with. 


Usage example:
<pre>
python PieParty.py -g expression_file.csv -c cell_coordinates.csv -l genelist1.csv genelist2.csv -gp True -cf clusters.csv -lc False
</pre>


### Example data files:
You can download a small example dataset from this github page (example dataset folder), which includes three files neccessary to run PieParty. Execute PieParty with:
<pre>
python PieParty.py -c coordinates_anon.csv -l example_gene_list.csv -g expression_anon.csv
</pre>
This should run in about 5 minutes and produce the following output png containing 183 cells:

<img src="https://github.com/harbourlab/PieParty/blob/master/example_output.png" width="300">


