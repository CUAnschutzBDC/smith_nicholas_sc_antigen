1.	Guidance on which QC plots are necessary to show
nFeatures, nCount, percent.mt per sample is always good
Maybe show all celltypes captured?
Definitely need the sorting strategy

2.	Updated UMAP of all cells clustered according to mRNA for cell type
(Naïve, Transitional (intermediate), Resting Memory, Activated Memory, Plasmablast)
I'm pretty confused here
My cell types:
B. Intermediate -- Trasitional intermediate
Plasmablast -- plasmablast
Naive_3 -- Naive
Naive_1 -- Naive
Resting memory -- memory
Memory_IgA -- memory
BND2 -- memory

What gets clumped?
BND, B. Intermediate into Transitional
Memory_IgA into ___ memory?
Combine Naive

3.	List of top 10 genes contributing to cell type determination
  * Not these are not really genes contributing to cell type determination as this is not how cell type was determined (this is important for you when writing the text), these are just the top 10 marker genes of each cluster.
4.	UMAP - Antigen reactive cells cluster independently (not clustering by mRNA)
INS/GAD/IA2/Multi-islet/TET/DNA
  * We decided thta this isn't a necessary plot
5.	Stacked bar chart: x axis = B cell subtype (updated), y axis = % of total cells minus other_multi_reactive
(INS/GAD/IA2/Multi-islet/TET/DNA/Negative)
6.	INSR mRNA expression violin plot separated by B cell subtype
7.	INSR mRNA expression violin plot separated by status
8.	IGHM mRNA expression violin plot separated by status
9.	IGHD mRNA expression violin plot separated by status
10.	IGHG mRNA expression violin plot separated by status
11.	IGHA mRNA expression violin plot separated by status
*make y axis of plots 7-11 consistent to show difference in values*
12.	New DEG comparisons heatmaps (first pass plots, will narrow down later)
13.	BCR signaling pathway (KEGG) genes heatmap

Can you send me this gene list?

14.	Stacked bar chart of all cells: y axis = isotype frequency, x axis = disease status 
15.	Stacked bar chart of islet-reactive cells: y axis = isotype frequency, x axis = disease status 
16.	SHM histogram - all cells
17.	SHM histogram - islet-reactive cells

16 and 17 - Separate by status

18.	Shannon alpha diversity rarefaction curve of all cells by status
19.	Jaccard beta diversity rarefaction curve of islet-reactive cells by status

Why shannon for all and jaccard for islet-reactive?

20.	Count number of islet-reactive clones found within each individual, make a stacked bar chart delineating whether the clones were public (black) or private (white), and group the individuals on the X axis according to disease status. Y axis is the number of clones present according to category. 

Do you want total number of clones or total number of cells that are part of clones? Ignore the size of the clone

21.	Do the same as 20, but instead of raw # of clones, normalize to a percent of total clones found in each individual, so if 4 clones total and 1 was public, 0.25 of the bar is black and 0.75 of the bar is white.

Again this ignores the size of the clone?

22.	Stacked bar chart of frequency of antigen reactivity type where x axis = disease status and y axis = INS/GAD/IA2/Multi-Islet TET/DNA
23.	Same as 22 but collapse all islet-reactive to one category so y axis = All-Islet TET/DNA
24.	Same as 22 & 23 but include the negative fraction of cells
25.	Same as 22 & 23 but include the other-reactive fraction of cells
26.	Same as 22 & 23 but include the other-reactive and negative fractions of cells
27.	Replicate something similar to [Figure 6 of this paper](https://journals.aai.org/jimmunol/article/198/4/1460/109668/Dysregulation-of-B-Cell-Repertoire-Formation-in) for our samples with all cells (by status)
28.	Replicate something similar to Figure 6 of this paper for our samples with islet-reactive cells (by status)
29.	Replicate something similar to [Figure 6 of this paper](https://journals.aai.org/jimmunol/article/198/4/1460/109668/Dysregulation-of-B-Cell-Repertoire-Formation-in) for our samples with all clones (by status)
30.	Replicate something similar to [Figure 6 of this paper](https://journals.aai.org/jimmunol/article/198/4/1460/109668/Dysregulation-of-B-Cell-Repertoire-Formation-in) for our samples with public clones (by status)
31.	Replicate something similar to [Figure 6 of this paper](https://journals.aai.org/jimmunol/article/198/4/1460/109668/Dysregulation-of-B-Cell-Repertoire-Formation-in) for our samples with private clones (by status)
