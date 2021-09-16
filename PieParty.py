#!/usr/bin/env python
# coding: utf-8

# In[1]:


color_greys_as_clusters = True  # this will be used to color each cluster uniquely down the road


# for testing purposes in e.g. Jupyter
#import ipywidgets as widgets   # this is for the progress bar
#expression_file = "20200429_scUM_Aggregate_expression.csv"
#coordinates_file = "20200506_scUM_Aggregate_TSNE_coordinates.csv"
#genes_files = ["Class1_up_genes.csv", "Class2_up_genes.csv"]
#cluster_information_file = "clusters.csv"
#giant_pies = 'True'
#output_filename = "output.png"
#colors = ["winter", "autumn"]
#lighten_colors = 'False'
#global_expression_colors = 'False'
#plot_size_in_pixels = 11800
#expression_cutoff_proportional = 1
#expression_cutoff_value = 0
#proportionalize = 'True'
#color_greys_as_clusters = 'True'


## (I)  the basics.. import, make matrices..
print("PieParty Version 1.7.2 starting")
print("... Note: You might get a system warning (bomb DOS attack warning), which is due to the high resolution of the output image. Can be ignored")

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

def plotpie(sizes, colors, coords, piechart_size):
    plt.pie(sizes, colors=colors)
    plot_output_filename = "temp_plot_" + output_filename
    plt.savefig(plot_output_filename, dpi=600, format="png", transparent=True)
    plt.clf()
    piechart = Image.open(plot_output_filename)
    piechart.thumbnail(piechart_size, Image.ANTIALIAS)  # size is the pie chart size
    #background.paste(piechart, (coords[0], coords[1]), piechart)
    background.paste(piechart, (int(float(coords[0]) - (piechart_size[0] * 0.510)), int(float(coords[1]) - (piechart_size[1] * 0.331))), piechart)
    piechart.close()
    os.remove(plot_output_filename)

def plotpie_grey(coords):
    piechart.thumbnail(size, Image.ANTIALIAS)  # size is the pie chart size
    background.paste(piechart, (int(float(coords[0])-(size[0]*0.510)), int(float(coords[1])-(size[1]*0.331))), piechart)         
    
def write_to_labels_file(row):
    with open(labels_file_name, "a") as f:
        f.write(row)
        f.write("\n")

def draw_rect(x_coord, y_coord, color, width, hight1, opacity):
    return '''<rect x="''' + str(x_coord) + '''" opacity="''' + str(opacity) + '''" y="''' + str(y_coord - hight1) + '''" fill="''' + color + '''" width="''' + str(width) + '''" height="''' + str(hight1) + '''"/>'''

def make_the_labels_file():
    if os.path.exists(labels_file_name):
        os.remove(labels_file_name)
    write_to_labels_file('''<svg viewBox="0 0 320 1000" xmlns="http://www.w3.org/2000/svg">''')
    y_coord = 20

    for y, gene_list in enumerate(gene_lists):
        y_coord += 15
        for x, gene in enumerate(gene_list):
            percentile = (x / (len(gene_list)))  # this picks the base color on a spectrum
            color = pick_color(colors[y], percentile, 1)
            write_to_labels_file(draw_rect(20, y_coord, color, 10, 10, "100"))
            if lighten_colors == 'True':
                color = pick_color(colors[y], percentile, 0.66)
                write_to_labels_file(draw_rect(30, y_coord, color, 10, 10, "100"))
                color = pick_color(colors[y], percentile, 0.25)
                write_to_labels_file(draw_rect(40, y_coord, color, 10, 10, "100"))
                write_to_labels_file('''<text text-anchor="start" font-family="Arial" x="55" y="''' + str(y_coord-2) + '''" font-size="8" >''' + gene + '''</text>''')
            else:
                write_to_labels_file('''<text text-anchor="start" font-family="Arial" x="35" y="''' + str(y_coord-2) + '''" font-size="8" >''' + gene + '''</text>''')
            y_coord += 15
    write_to_labels_file("</svg>")

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
parser.add_argument('-gc', '--lighten_colors_based_on_global', help='True or False. Lighten colors based on global gene expression, default: True', required=False, type=str, default='True')
parser.add_argument('-pz', '--plot_resolution', help='Plot resolution in pixels. Default 11800', required=False, type=int, default=11800)
parser.add_argument('-gp', '--giant_pies', help='Plot one pie chart per cluster', required=False, type=str, default='False')
parser.add_argument('-cf', '--cluster_file', help='File with cluster information for every cell', required=False, type=str)
parser.add_argument('-ps', '--pie_sizes_scaling_factor', help='multiplies the standard pie size by this factor', required=False, type=float, default=1)
args = vars(parser.parse_args())

expression_file = args['matrix_expression_file'] #"20200429_scUM_Aggregate_expression.csv"
pie_size_scaling_factor = args['pie_sizes_scaling_factor']
coordinates_file = args['matrix_coordinates'] #"20200429_scUM_Aggregate_TSNE_coordinates_Class2only.csv"
genes_files = args['gene_lists'] #["Class2_up_genes.csv", "Class1_up_genes.csv"]
cluster_information_file = args['cluster_file']
giant_pies = args['giant_pies']
colors = args['colors']  # ["@FF0000", "@0000FF"]
output_filename = args['output_filename']#"lighten_binary_class2.png"
if output_filename[-4:] != ".png":
    output_filename += ".png"
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
        os.remove(plot_output_filename)
        sys.exit("Error: If giant pie plotting is chosen a cluster file has to be provided.")

if os.path.exists(output_filename):
    os.remove(output_filename)

labels_file_name = "labels_" + output_filename[:-4] +".svg"
###################################### this prepares all the matrices used #################################################
temp_gene_lists = [csv_to_array(i) for i in genes_files]
gene_lists = []
for x, j in enumerate(temp_gene_lists):
    gene_lists.append([i[0] for i in j])
if len(colors) != len(gene_lists):
    os.remove(plot_output_filename)
    sys.exit("Error: Number of colors must be equal to number of gene lists")

matrix_expression = csv_to_array(expression_file)
coordinates = csv_to_array(coordinates_file)
coordinates_sans_1 = coordinates[1:]  # exclude first row (header)

### get_x_and_y_min_and_max
x_values = []
y_values = []
for i in range(len(coordinates_sans_1)):
    x_values.append(float(coordinates_sans_1[i][1]))
    y_values.append(float(coordinates_sans_1[i][2]))
x_min = min(x_values)
x_max = max(x_values)
y_min = min(y_values)
y_max = max(y_values)
x_stretch = x_max - x_min
y_stretch = y_max - y_min
x_correction_factor = plot_size_in_pixels/x_stretch
y_correction_factor = plot_size_in_pixels/y_stretch

#make empty background
img = Image.new('RGB', (int(plot_size_in_pixels*1.1), int(plot_size_in_pixels*1.1)), (255, 255, 255))  # *1.1 to add a white border
img.save(output_filename, "PNG")

small_pie_size = pie_size_scaling_factor*30*plot_size_in_pixels/2400, 30*plot_size_in_pixels/2400  # size of pies

print("... Gathered all ingredients, mixing the dough now")


# In[2]:


###### if big pies should be plotted, we are making new "fake" matrices here to feed into the pipeline #############

if giant_pies == 'True':
    # 1) make a list of all clusters
    unique_clusters = []
    for x, i in enumerate(cluster_information):
        if x > 0: # exclude header in first row
            if str(i[1]) not in unique_clusters:
                unique_clusters.append(str(i[1])) # str allows for various cluster names           
#check

    # 2) Get coordinates of all cells per cluster, to make an average, keeps sorting as in unique_clusters
    all_cells_coorinates_by_cluster = [[[],[]] for i in range(len(unique_clusters))]  # list of all x and y coordinates for each cluster [[[x, x, ..], [y, y, ..]], ...]
    cell_IDs = [i[0] for i in coordinates]
    all_cells_in_clusters = [[] for i in range(len(unique_clusters))]  # list of cell_IDs in clusters
    for x, i in enumerate(cluster_information): # [cellid, cluster] check all cells in this file and assign spot in all_cells_coorinates_by_cluster
        if x > 0:
            if i[0] in cell_IDs:
                all_cells_coorinates_by_cluster[unique_clusters.index(i[1])][0].append(float(coordinates[cell_IDs.index(i[0])][1]))    # x-coord, keeps  sorting..
                all_cells_coorinates_by_cluster[unique_clusters.index(i[1])][1].append(float(coordinates[cell_IDs.index(i[0])][2]))    # y-coord
                all_cells_in_clusters[unique_clusters.index(i[1])].append(i[0])

    # 3) make average coordinates for each cluster, keeps sorting as i unique_clusters
    average_coordinates = [['', 'tSNE_1', 'tSNE_2']]
    for t, i in enumerate(all_cells_coorinates_by_cluster):
        temp = []
        temp.append(unique_clusters[t])
        for j in i:
            temp.append(str(np.average(j)))
        average_coordinates.append(temp)     
    
    # 4) Make new fake expression matrix
    new_fake_expression_matrix = []
    temp = copy.copy(unique_clusters)
    temp.insert(0, "")
    new_fake_expression_matrix.append(temp)  # add all clusters (new fake cells), and a empty space on top

    for x, gene_column in enumerate(matrix_expression):  # go through all genes in original matrix
        if x > 0: # first column are the cell IDs.
            current_gene = gene_column[0]
            if any(current_gene in sublist for sublist in gene_lists):
                new_gene_entry = []  # final list to be added to new_fake_gene_expression_matrix
                new_gene_entry.append(current_gene)  # fist is the gene name
                for cells_list_in_cluster in all_cells_in_clusters:  # go through each list of cells (in each cluster)
                    temp_expr_values_for_current_cluster = []
                    for num, cell_ID in enumerate(cells_list_in_cluster):
                        cell_index_in_original_expr_file = matrix_expression[0].index(cell_ID)
                        temp_expr_values_for_current_cluster.append(float(gene_column[cell_index_in_original_expr_file]))
                    new_gene_entry.append(np.average(temp_expr_values_for_current_cluster)/len(cells_list_in_cluster))
                new_fake_expression_matrix.append(new_gene_entry)

    # 5) this is where the magic happens, but first add all grey pies
    grey_pies = []  # [[color, [x, y]], [color, [x, y]] ..]
    background = Image.open(output_filename)
    size = small_pie_size # defines pie size
    for x, i in enumerate(coordinates_sans_1):
        pie_coordinates = [int((float(i[1]) - x_min) * x_correction_factor), int(plot_size_in_pixels - (float(i[2]) - y_min) * y_correction_factor)]
        if color_greys_as_clusters == 'True':
            grey_pies.append([["#CCCCCC"], pie_coordinates])
        else:
            grey_pies.append([["#CCCCCC"], pie_coordinates]) # find a way to put custom colors per pie here
    print("... plotting the grey cells as a backgroud for the big pies")
    #progress_bar = widgets.FloatProgress(value=0, min=0, max=100.0, step=0.1, description='Plotting grey pies:', bar_style='info', orientation='horizontal')
    #display(progress_bar)
    
    plt.pie([1], colors=["#CCCCCC"])
    plot_output_filename = "temp_plot_grey_" + output_filename
    plt.savefig(plot_output_filename, dpi=600, format="png", transparent=True)
    plt.clf()
    piechart = Image.open(plot_output_filename)
    for pie_nr, i in enumerate(grey_pies):
        plotpie_grey(i[1])
        #progress_bar.value = (100*pie_nr)/len(grey_pies)
    piechart.close()
    os.remove(plot_output_filename)
    background.save(output_filename)
    matrix_expression = new_fake_expression_matrix
    coordinates = average_coordinates
    coordinates_sans_1 = coordinates[1:]
    print("... prepared the large bowls for the big pies")
    


# In[3]:


# code starts here

if lighten_colors == 'True':
    expression_values = [[] for i in range(len(genes_files))]  # make empty list for each gene file
    for gene_row in matrix_expression:  # go through all genes in dataset
        for list_nr, gene_list in enumerate(gene_lists):
            if gene_row[0] in gene_list:
                expression_values[list_nr].append(gene_row[1:])  # adds list of all expression values for each gene and all cells without gene name
    list_all_genes = []  # this condenses the list of lists to a list with all expression values for all genes in each group/gene list
    for sublist in expression_values:
        temp_list = []
        for gene in sublist:
            temp_list += gene
        list_all_genes.append(temp_list)

    if global_expression_colors == 'True':
        temp_list = []
        for i in list_all_genes:
            temp_list += i
        list_all_genes = [copy.copy(temp_list)]  # makes one giant list

    cutoffs = []  # these are the cutoff values (top 10%). Everything above that value is max color for the respective list
    for gene_list in list_all_genes:
        gene_list_floats = [float(i) for i in gene_list]
        sorted_gene_list = sorted(gene_list_floats)
        while sorted_gene_list[0] == 0:  # removes 0 values
            sorted_gene_list.pop(0)
        cutoffs.append(sorted_gene_list[(int(len(sorted_gene_list)*0.9))])

grey_pies_coordinates = []  # [[coordinates]]
colored_pies = []  # [[[piezises], [colors], [x, y], ...]

expression_sans_1 = matrix_expression[1:]  # exclude first row (header)
all_genes_list = [a[0] for a in expression_sans_1]  # make list of all genes


for x, i in enumerate(coordinates_sans_1):  # iterate through all cells
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
                os.remove(plot_output_filename)
                sys.exit('''The gene "''' + current_CTA + '''" was not found in expression dataset. Aborting.''')

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
        pie_coordinates = [int((float(i[1]) - x_min) * x_correction_factor), int(plot_size_in_pixels - (float(i[2]) - y_min) * y_correction_factor)]  #large y values are high up
        grey_pies_coordinates.append(pie_coordinates)

print("... done mixing the ingredients. Baking the pies now!")

background = Image.open(output_filename)
size = small_pie_size

for i in grey_pies_coordinates:
    plotpie([1], ["#CCCCCC"], i, small_pie_size)

for x, i in enumerate(colored_pies):
    if giant_pies == 'False':
        plotpie(i[0], i[1], i[2], small_pie_size)
    elif giant_pies == 'True':
        custom = int(small_pie_size[0] * np.sqrt(len(all_cells_in_clusters[x])))/4
        custom_size = custom, custom
        plotpie(i[0], i[1], i[2], custom_size)
    else:
        os.remove(plot_output_filename)
        sys.exit("giant_pies must be True or False")
background.save(output_filename)

print("... " + str(x+1) + " pies were baked")

# Make the labels file
print("... Generating the labels file")
make_the_labels_file()

print("All done! Enjoy!")


# In[ ]:




