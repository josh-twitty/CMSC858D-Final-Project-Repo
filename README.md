This repo contains the script used in this project and the results obtained in it.

There is only one script used, contained in the topmost level folder named make_plots.py. This script is intended to be used with 
the alevin folder that salmon produces after running alevin-fry on the data, and it runs both the gene count and equivalence class count matrix representation through the data analysis routine I used and produces a folder full of several plots and a text document of clustering scores in the output folder. Usage for this is

```
python make_plots.py [ALEVIN FOLDER PATH] [OUTPUT_FOLDER] [PLOT_FILENAME_BASE]
```

Where `[ALEVIN FOLDER PATH]` is the path to the alevin folder, `[OUTPUT_FOLDER]` is the path you want to create the results in, and `[PLOT_FILENAME_BASE]` is the prefix you want each plot's filename to have.

[Salmon](https://combine-lab.github.io/salmon/) was used to align the raw sequencing data. 

In the `results` folder of this repo is the output of this script on several datasets.

Software versions I used:
Python 3.10.4
numpy 1.21.6
scanpy 1.9.1
sklearn 1.0.2
raw data processed using salmon 1.8.0

Datasets used:
500 1:1 Mixture of Human HEK293T and Mouse NIH3T3 cells, 3' LT v3.1, Chromium X (hm_500_v3 in this repo)
```
https://www.10xgenomics.com/resources/datasets/500-1-1-mixture-of-human-hek-293-t-and-mouse-nih-3-t-3-cells-3-lt-v-3-1-chromium-x-3-1-low-6-1-0
```
500 Human PBMCs, 3' LT v3.1, Chromium X (pbmc_500_v3 in this repo)
```
https://www.10xgenomics.com/resources/datasets/500-human-pbm-cs-3-lt-v-3-1-chromium-x-3-1-low-6-1-0
```
1k PBMCs from a Healthy Donor (v3 chemistry) (pbmc1k in this repo)
```
https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0
```
