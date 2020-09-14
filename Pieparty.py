

## (I)  the basics.. import, make matrices..

print("PieParty Version 1.5 starting")

import matplotlib.pyplot as plt
import matplotlib as pl
import os
import csv
from PIL import Image
import argparse
import sys
import copy
import colorsys
import numpy as np

def lighten_color(color, amount=0): # percentage lightening. Put number between 0 and 1. 1 is no change, 0 is white
    try:
        c = pl.colors.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*pl.colors.to_rgb(c))
    return pl.colors.to_hex(colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2]))

def write_to_file(row):
    with open(output_filename, "a") as f:
        f.write(row)
        f.write("\n")

def pick_color(color, percentile, lighten):
    if color[0] == "@":
        return_color = ("#" + color[1:])
    else:
        cmap = plt.cm.get_cmap(color)
        rgba = cmap(percentile)
        return_color = pl.colors.to_hex(rgba)
    if lighten < 1:  # 1 is max..
        return(lighten_color(return_color, lighten))
    else:
        return(return_color)

def csv_to_array(filename):
    with open(filename, 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        return list(csv_reader)

def plotpie(sizes, colors, coords):
    plt.pie(sizes, colors=colors)
    plot_output_filename = "temp_plot_" + output_filename
    plt.savefig(plot_output_filename, dpi=600, format="png", transparent=True)
    plt.clf()
    piechart = Image.open(plot_output_filename)
    piechart.thumbnail(size, Image.ANTIALIAS)
    background.paste(piechart, (coords[0], coords[1]), piechart)
    piechart.close()
    os.remove(plot_output_filename)

parser = argparse.ArgumentParser(description='PieParty_args')
parser.add_argument('-g', '--matrix_expression_file', help='gene expression file', required=True, type=str)
parser.add_argument('-c', '--matrix_coordinates', help='coordinate file for all cells', required=True, type=str)
parser.add_argument('-l', '--gene_lists', help='lists of genes to plot. Can be one or multiple', required=True, nargs='+', type=str)
parser.add_argument('-o', '--output_filename', help='output file name', required=False, type=str, default="output.png")
parser.add_argument('-color', '--colors', help='list of colors (one per gene lost) or matplotlib color spectrum. default is "viridis", alternatively hex colors can be provided e.g. "@FFFFFF"', required=False, type=str, nargs='+', default=["viridis"])
parser.add_argument('-p', '--proportionalize', help='True or False. normalize for number of genes in lists, default True', required=False, type=str, default='True')
parser.add_argument('-ct', '--cutoff_percentage', help='percentage cutoff to be included in pie charts. Default 1%', required=False, type=float, default=1.0)
parser.add_argument('-ce', '--cutoff_expression_value', help='Cutoff of expression value can be set. Default 0', required=False, type=float, default=0)
parser.add_argument('-lc', '--lighten_colors', help='True or False. Lighten colors in pie charts based on gene expression. Default True.', required=False, type=str, default='True')
parser.add_argument('-gc', '--lighten_colors_based_on_global', help='True or False. Lighten colors based on global gene expression, default False', required=False, type=str, default='False')
parser.add_argument('-pz', '--plot_resolution', help='Plot resolution in pixels. Default 13000', required=False, type=int, default=13000)
parser.add_argument('-gp', '--giant_pies', help='Plot one pie chart per cluster', required=False, type=str, default='False')
parser.add_argument('-cf', '--cluster_file', help='File with cluster information for every cell', required=False, type=str)
args = vars(parser.parse_args())

expression_file = args['matrix_expression_file'] #"20200429_scUM_Aggregate_expression.csv"
coordinates_file = args['matrix_coordinates'] #"20200429_scUM_Aggregate_TSNE_coordinates_Class2only.csv"
genes_files = args['gene_lists'] #["Class2_up_genes.csv", "Class1_up_genes.csv"]
cluster_information_file = args['cluster_file']
giant_pies = args['giant_pies']
colors = args['colors'] #["@FF0000", "@0000FF"]
output_filename = args['output_filename']#"lighten_binary_class2.png"
proportionalize = args['proportionalize'] #True  # in case unequal amounts of genes are present in lists, and should be weighted the same (will devide epresion by sample size in file)
expression_cutoff_proportional = args['cutoff_percentage'] # 1  # percentage. needs to be this percentage of total pie chart to be included
expression_cutoff_value = args['cutoff_expression_value'] #0
lighten_colors = args['lighten_colors'] # lightens the colors depending on their expression
global_expression_colors = args['lighten_colors_based_on_global'] # lightens the color based on all genes in the single cell gene expression dataset.
plot_size_in_pixels = args['plot_resolution']  # 2400 is good enough in most cases

if giant_pies == 'True':
    if cluster_information_file is not None:
        cluster_information = csv_to_array(cluster_information_file)
    else:
        sys.exit("Error: If giant pie plotting is chosen a cluster file has to be provided.")

if os.path.exists(output_filename):
    os.remove(output_filename)

###################################### this prepares all the matrices used #################################################
temp_gene_lists = [csv_to_array(i) for i in genes_files]
gene_lists = []
for x, j in enumerate(temp_gene_lists):
    gene_lists.append([i[0] for i in j])
if len(colors) != len(gene_lists):
    sys.exit("Error: Number of colors must be equal to number of gene lists")
matrix_expression = csv_to_array(expression_file)
coordinates = csv_to_array(coordinates_file)

print("... Gathered all ingredients, mixing the dough now")


# In[3]:


###### if big pies should be plotted, we are making new "fake" matrices here to feed into the pipeline #############

if giant_pies == 'True':
    # 1) make a list of all clusters
    unique_clusters = []
    for x, i in enumerate(cluster_information):
        if x > 0: # exclude header in first row
            if str(i[1]) not in unique_clusters:
                unique_clusters.append(str(i[1])) # str allows for various cluster names

    # 2) Get coordinates of all cells per cluster, to make an average, keeps sorting as in unique_clusters
    giant_pies_coordinates = [[[],[]] for i in range(len(unique_clusters))] # list of all x and y coordinates for each cluster [[[x, x, ..], [y, y, ..]], ...]
    cell_IDs = [i[0] for i in coordinates]
    all_cells_in_clusters = [[] for i in range(len(unique_clusters))]  # list of cell_IDs in clusters
    for x, i in enumerate(cluster_information): # check all cells in this file and assign spot in giant_pies_coordinates
        if x > 0:
            if i[0] in cell_IDs:
                if i[0] != "":
                    giant_pies_coordinates[unique_clusters.index(i[1])][0].append(float(coordinates[cell_IDs.index(i[0])][1]))    # x-coord, keeps original sorting..
                    giant_pies_coordinates[unique_clusters.index(i[1])][1].append(float(coordinates[cell_IDs.index(i[0])][2]))    # y-coord
                    all_cells_in_clusters[unique_clusters.index(i[1])].append(i[0])

    # 3) make average coordinates for each cluster, keeps sorting as i unique_clusters
    average_coordinates = [['', 'tSNE_1', 'tSNE_2']]
    for cluster, i in enumerate(giant_pies_coordinates):
        temp = []
        temp.append(str(cluster))
        for x, j in enumerate(i):
            temp.append(str(np.average(j)))
        average_coordinates.append(temp)

    # 4) Make new fake expression matrix
    new_fake_expression_matrix = []
    temp = copy.copy(unique_clusters)
    temp.insert(0,"")
    new_fake_expression_matrix.append(temp) # add all clusters (new fake cells), and a empty space on top

    for x, gene_column in enumerate(matrix_expression): # go through all genes in original matrix
        if x > 0: # first column is cell IDs.
            current_gene = gene_column[0]
            if any(current_gene in sublist for sublist in gene_lists):
                temp_expr_values_for_one_gene = []  # final to be added to new_fake_gene_expression_matrix
                temp_expr_values_for_one_gene.append(current_gene)
                for cells_list_in_cluster in all_cells_in_clusters:   # go through each list of cells (in each cluster)
                    temp_expr_values_for_current_cluster = []
                    for num, cell_ID in enumerate(cells_list_in_cluster):
                        cell_index_in_original_expr_file = matrix_expression[0].index(cell_ID)
                        temp_expr_values_for_current_cluster.append(float(gene_column[cell_index_in_original_expr_file]))
                    temp_expr_values_for_one_gene.append(np.average(temp_expr_values_for_current_cluster))
                new_fake_expression_matrix.append(temp_expr_values_for_one_gene)

# this is where the magic happens
    matrix_expression = new_fake_expression_matrix
    coordinates = average_coordinates

    print("... prepared the large bowls for the big pies")


# In[4]:


if lighten_colors == 'True':
    expression_values = [[] for i in range(len(genes_files))]
    for gene_row in matrix_expression:
        for list_nr, gene_list in enumerate(gene_lists):
            if gene_row[0] in gene_list:
                expression_values[list_nr].append(gene_row[1:])  # adds all expression values without gene name

    list_all_genes = []  # one list for each gene file containing all expression values
    for sublist in expression_values:
        temp_list = []
        for gene in sublist:
            temp_list += gene
        list_all_genes.append(temp_list)
    if global_expression_colors == 'True':
        temp_list = []
        for i in list_all_genes:
            temp_list += i
        list_all_genes = [copy.copy(temp_list)]

    cutoffs = []  # these are the cutoff values (top 10%). Everything above that value is max color for the respective list
    for gene_list in list_all_genes:
        gene_list_floats = [float(i) for i in gene_list]
        sorted_gene_list = sorted(gene_list_floats)
        while sorted_gene_list[0] == 0:
            sorted_gene_list.pop(0)
        cutoffs.append(sorted_gene_list[(int(len(sorted_gene_list)*0.9))])

grey_pies = []  # [[[piezises], [colors], [x, y], ...]
colored_pies = []  # [[[piezises], [colors], [x, y], ...]


### get_x_and_y_min_and_max
x_values = []
y_values = []
for i in range(len(coordinates)):
    if i > 0:
        x_values.append(float(coordinates[i][1]))
        y_values.append(float(coordinates[i][2]))
x_min = min(x_values)
x_max = max(x_values)
y_min = min(y_values)
y_max = max(y_values)
x_stretch = x_max - x_min
y_stretch = y_max - y_min
x_correction_factor = plot_size_in_pixels/x_stretch
y_correction_factor = plot_size_in_pixels/y_stretch

coordinates_sans_1 = coordinates[1:]  # exclude first row (header)
expression_sans_1 = matrix_expression[1:]  # exclude first row (header)
all_genes_list = [a[0] for a in expression_sans_1]  # make list of all genes


#make empty background
img = Image.new('RGB', (plot_size_in_pixels, plot_size_in_pixels), (255, 255, 255))
img.save(output_filename, "PNG")
size = 30*plot_size_in_pixels/2400, 30*plot_size_in_pixels/2400  # size of pies


# iterate through all cells
for x, i in enumerate(coordinates_sans_1):
    current_cell_ID = i[0]
    pie_sizes = []
    pie_colors = []
    pie_number_of_genes_in_list = []  # important if data needs to be proportionalized
    cell_ID_location = matrix_expression[0].index(current_cell_ID)

    for file_nr, matrix_order in enumerate(gene_lists):  # iterate through gene lists
        for CTA_nr, current_CTA in enumerate(matrix_order):
            try:
                expression_row = all_genes_list.index(current_CTA)
            except:
                sys.exit('''The gene "''' + current_CTA + '''" was not found in expression dataset''')
            expression_value = float(expression_sans_1[expression_row][cell_ID_location])
            pie_sizes.append(expression_value)
            percentile = (CTA_nr / (len(matrix_order))) * 0.999  # color range percentile has to be from 0 - < 1

            if lighten_colors == 'True':
                if global_expression_colors == 'True':
                    lighten = expression_value / float(cutoffs[0])
                else:
                    lighten = expression_value/float(cutoffs[file_nr])
            else:
                lighten = 1
            pie_colors.append(pick_color(colors[file_nr], percentile, lighten))
            pie_number_of_genes_in_list.append(len(matrix_order))

### make all the cutoffs
    original_pie_sizes = copy.deepcopy(pie_sizes)
    original_number_of_genes_in_list = copy.deepcopy(pie_number_of_genes_in_list)
    if expression_cutoff_proportional != 0:
        pie_sizes2 = []
        pie_colors2 = []
        pie_number_of_genes_in_list2 = []
        percentage_cutoff = sum(pie_sizes)*expression_cutoff_proportional/100
        for r, t in enumerate(pie_sizes):
            if t > percentage_cutoff:
                pie_sizes2.append(t)
                pie_colors2.append(pie_colors[r])
                pie_number_of_genes_in_list2.append(pie_number_of_genes_in_list[r])
        pie_sizes = copy.deepcopy(pie_sizes2)
        pie_colors = copy.deepcopy(pie_colors2)
        pie_number_of_genes_in_list = copy.deepcopy(pie_number_of_genes_in_list2)
    if expression_cutoff_value != 0:
        pie_sizes2 = []
        pie_colors2 = []
        pie_number_of_genes_in_list2 = []
        for r, t in enumerate(pie_sizes):
            if t > percentage_cutoff_value:
                pie_sizes2.append(t)
                pie_colors2.append(pie_colors[r])
                pie_number_of_genes_in_list2.append(pie_number_of_genes_in_list[r])
        pie_sizes = copy.deepcopy(pie_sizes2)
        pie_colors = copy.deepcopy(pie_colors2)
        pie_number_of_genes_in_list = copy.deepcopy(pie_number_of_genes_in_list2)
    if proportionalize == 'True':
        pie_sizes2 = []
        for r, t in enumerate(pie_sizes):
            pie_sizes2.append(t/pie_number_of_genes_in_list[r])
        pie_sizes = copy.deepcopy(pie_sizes2)
    if sum(pie_sizes) > 0:
        pie_sizes2 = []
        pie_colors2 = []
        pie_number_of_genes_in_list2 = []
        for r, t in enumerate(pie_sizes):
            if t > 0:  # exclude all without expression
                pie_sizes2.append(t)
                pie_colors2.append(pie_colors[r])
                pie_number_of_genes_in_list2.append(pie_number_of_genes_in_list[r])
        pie_sizes = copy.deepcopy(pie_sizes2)
        pie_colors = copy.deepcopy(pie_colors2)
        pie_number_of_genes_in_list = copy.deepcopy(pie_number_of_genes_in_list2)
        while sum(pie_sizes) < 1:  # to plot pie charts the sum has to be at least 1, else the pie is incomplete
            pie_sizes = [s*2 for s in pie_sizes]
        pie_coordinates = [int((float(i[1]) - x_min) * x_correction_factor), int(plot_size_in_pixels - (float(i[2]) - y_min) * y_correction_factor)]  # x , y
        colored_pies.append([pie_sizes, pie_colors, pie_coordinates])
    else:
        pie_coordinates = [int((float(i[1]) - x_min) * x_correction_factor), int(plot_size_in_pixels - (float(i[2]) - y_min) * y_correction_factor)]
        grey_pies.append([[1], ["#CCCCCC"], pie_coordinates])



print("... Done crunching the ingredients. Baking the pies now!")
print("... Note: Resolution of PieParty images is very high, you might get a system warning (bomb DOS attack warning). Trust me.. it's fine.")
background = Image.open(output_filename)
for i in grey_pies:
    plotpie(i[0], i[1], i[2])
for i in colored_pies:
    plotpie(i[0], i[1], i[2])
background.save(output_filename)
print("... All done! " + str(x+1) + " pies were baked. Enjoy!")


# In[10]:


# Make the labels file

if os.path.exists("labels.svg"):
    os.remove("labels.svg")


def write_to_file(row):
    with open("labels.svg", "a") as f:
        f.write(row)
        f.write("\n")

def draw_rect(x_coord, y_coord, color, width, hight1, opacity):
    return '''<rect x="''' + str(x_coord) + '''" opacity="''' + str(opacity) + '''" y="''' + str(y_coord - hight1) + '''" fill="''' + color + '''" width="''' + str(width) + '''" height="''' + str(hight1) + '''"/>'''


write_to_file('''<svg viewBox="0 0 320 1000" xmlns="http://www.w3.org/2000/svg">''')


y_coord = 20
for y, gene_list in enumerate(gene_lists):
    y_coord += 15
    for x, gene in enumerate(gene_list):
        percentile = (x / (len(gene_list)))
        color = pick_color(colors[y], percentile, 1)
        write_to_file(draw_rect(20, y_coord, color, 10, 10, "100"))
        if lighten_colors == 'True':
            color = pick_color(colors[y], percentile, 0.75)
            write_to_file(draw_rect(30, y_coord, color, 10, 10, "100"))
            color = pick_color(colors[y], percentile, 0.25)
            write_to_file(draw_rect(40, y_coord, color, 10, 10, "100"))
            write_to_file('''<text text-anchor="start" font-family="Arial" x="55" y="''' + str(y_coord-2) + '''" font-size="8" >''' + gene + '''</text>''')
        else:
            write_to_file('''<text text-anchor="start" font-family="Arial" x="35" y="''' + str(y_coord-2) + '''" font-size="8" >''' + gene + '''</text>''')

        y_coord += 15
write_to_file("</svg>")


# In[ ]:





# In[ ]:
