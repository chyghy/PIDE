# PIDE (Prophage Island Detection using ESM-2)



### Description

A framework for identification of prophage islands from bacterial genomes and metagenome-assembled genomes. It can also be used for identification of phage sequences from metagenomes.

### Instructions

(1) Create an environment for PIDE

```bash
conda create -n pide
conda activate pide
```

(2) Install the required Python packages

a, PyTorch:

Run on CPU

```bash
conda install pytorch torchvision torchaudio cpuonly -c pytorch
```

If you want to use the GPU version, please go to https://pytorch.org/get-started and get the conda or pip install command according to your device and demand.

⚠️ The GPU version of PIDE does not currently support macOS systems with Apple Silicon (M-series) chips by default. It is recommended to use the CPU version directly, or modify the relevant code to enable GPU support on M-series Macs.

b, fair-esm

``` 
pip install fair-esm
```

c, pandas

```bash
conda install pandas
```

d, biopython

```
conda install -c bioconda biopython
```

e, prodigal

```
conda install -c bioconda prodigal
```

(3) Download the model

```bash
wget https://zenodo.org/records/12759619/files/PIDE.model.tar.gz
tar xzvf PIDE.model.tar.gz
```

(4) Download the source code of PIDE from github

```
git clone https://github.com/chyghy/PIDE.git
```

### Usage

To get the HELP information

```bash
python PIDE/classification.py -h
```

```
python PIDE/classification.py [-o OUTPUT] [-g GPU] [-b BATCHSIZE] [-n MINPNUM] [-d DISTANCE] [-s PISCORE] [-m] input model
```

Explanation

```
positional arguments:
  input                 Path of the input fasta file
  model                 Path of the model parameter file

optional arguments:
  -h, --help            Show this help message
  -o OUTPUT, --output OUTPUT
                        Path of the output directory
  -g GPU, --GPU GPU     Determine which GPU(s) to use. If this parameter is not used, the GPU is used by default. Multi-GPU is also supported, IDs of different GPUs are separated by commas
  -b BATCHSIZE, --BatchSize BATCHSIZE
                        Define the batch size used in the prediction(default is 2). Note that the batch size cannot be negative and should not be smaller than the number of GPUs used
  -n MINPNUM, --MinPNum MINPNUM
                        The min prophage ORF number of a PI (default is 5)
  -d DISTANCE, --Distance DISTANCE
                        The clustering distance (bp) to use (default is 3000)
  -s PISCORE, --PIScore PISCORE
                        The threshold of PI score (default is 0.7)
  -m, --meta            Use meta mode during ORF prediction
```



### Output

1. ###### cluster.csv

   This csv file lists all the prophage islands that PIDE found in the input fasta file. Here is each column represents:

   **Contig**: The contig where the PIs is located.

   **Start**: The start site of the PI.

   **End**: The end site of the PI.

   **Score**: The PI score.

   **B**: The locations of the PI-carried bacteria genes.

   **Total_ORFs**: The total number of the PI-carried genes.

   **B_count**: The total number of the PI-carried bacteria genes.

   **B_ratio**: The ratio of the PI-carried bacteria genes.

2. ###### virus_contig.txt

   This csv file lists all the phage contigs that PIDE found in the input fasta file. Here is each column represents:

   **Contig_name**: The phage contig name.

   **Length**: The length of this contig.
   
   **Proportion**: The phage ORFs proportion of this contig.
