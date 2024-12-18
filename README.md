# Meta-SIFT

Scripts for "Engineering bacteriophages through deep mining of metagenomic motifs".  

### Required Dependencies:  
* Numpy  
* Cython  

### Setup and Installation:  
1. `git clone https://github.com/PhilHuss/MotifMining`  
2. `cd MotifMining`  
3. ` python3 compile_cython.py build_ext –inplace`  

### Usage:  
`motif_finder_tool.py -i <input_proteins> -m <matrix_table> [options]`  

### Publication Settings:  
* `-n 6 -c 1 -e 1e-50 -t 0.8`  
* `-n 10 -c 1 -e 1e-5 -t 0.045`  

### Files and Folders
* *motif_finder_tool.py*: main motif search tool
* *motif_finder_tool_modules.pyx*: Cython modules for main tool, requires compilation
* *compile_cython.py*: used to compile Cython modules
* *IMGVR_v4_human-animal-wastewater.VOG-relevant.cdhit100.faa.gz*: curated database proteins from IMG/VR (`-i`). Some proteins contained within this database are unpublished and are credited to the original authors. Please see the [IMG/VR website](https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html) and usage policy for details. 
* *prokaryotic_virus_ncbi_july2019.lim3000.trimmed.VOG-relevant.cdhit100.faa.gz*: curated database proteins from NCBI (`-i`). 
* *Hierarchical_Cluster.R*: R scripts used for hierarchical clustering
