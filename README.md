# TF-miRNA Feed Forward Loop

This script will calculate and export Transcription factor (TF) feed forward loops, miRNA 
feed forward loops, TFs-miRNAs composite feed forward loops and random permutation results.

### Prerequisites:
This script has been tested in Python 2.7 and require the following python modules: 

    1.optparse (version 2.3 or later)
    2.random (version 2.4 or later)
    3.numpy (version 1.10 or later)
    4.collections (version 2.5 or later)
### Usage

    Usage: FFL.py [options]

    Options:
     -h, --help       show this help message and exit
     -u, --usage      show more info on how to use this script
     -n NETWORK_FILE  regulatory network file with TFs/miRNAs at first column and
                      target gene at second column
     -m M             number of genes for random selection
     -g G             gene list for random selection
     -l L             file containing TFs and miRNAs
     -t TF_MI_GENE    TFs and miRNAs target all coding gene file with TFs/miRNAs
                      at first column target gene at second column
     -i TF_MI         file containing fixed TF-miRNA interaction
     -s SIMS          number of random simulations [default: 1000]
  
### Required input:

All files must be a tab-separated table with sources at the first column and targets at the second column.

For calculation of actural feed forward loop only, one network file is required:
    
    1. network file with TFs/miRNAs at the first column, target genes at the second column.
  
For calculation of actural feed forward loop and random permutation, you MUST provide all following files:

    1. network file with TFs/miRNAs at first column target and gene at second column (-n)
  
    2. gene list for random selection (-g)
  
    3. fixed TFs and miRNAs list (-l)
  
    4. TFs/miRNAs and all human coding gene regulation information with TFs/miRNAs at first column target the gene at second column (-t)
  
    5. fixed TFs and miRNAs regulation information with TFs/miRNAs at first column target the miRNAs/TF at second column(-i)

### Output files:

  Results will be saved as:

    1. TFs feed forward loop: TF_FFL_output_core.txt

    2. miRNAs feed forward loop: miRNA_FFL_output_core.txt

    3. TFs miRNAs composite feed forward loop: TF_miRNAs_composite_FFL_output_core.txt

    4. Random permutation loop counts: FFL_random_count.txt

    5. Z score: z_score_out.txt
    
### Example

Here's an example, provided the files are in the same directory as this python script:

    1. For calculation of actural feed forward loop only
    ./python FFL.py -n network.txt 

    2. Include random permutation
    ./python FFL.py -n network.txt -m 10 -g gene.txt -l fixed_TFs_miRNAs.txt -t TFs_miRNAs_to_gene.txt -i TFs_miRNAs_regulation.txt -s 1000

### Note
miRNAs must be provided as standard miRNAs ID, for example, hsa-miR-122a-5p or hsa-let-7a-5p.
