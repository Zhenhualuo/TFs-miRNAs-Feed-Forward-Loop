# TFs-miRNAs-Feed-Forward-Loop

This program will calculate and export Transcription factors (TFs) feed forward loop, miRNAs 
feed forward loop and TFs-miRNAs composite feed forward loop, random permutation.

Required input:

    All files must be a tab-separated table.

For calculation of actural feed forward loop only:
    
    1. network file with TFs/miRNAs at first column, target the gene at second column
  
For calculation of actural feed forward loop and random permutation, you MUST provide all following files:

    1. network file with TFs/miRNAs at first column target and gene at second column (-n)
  
    2. gene list for random selection (-G)
  
    3. fixed TFs and miRNAs list (--tm)
  
    4. TFs/miRNAs and all human coding gene regulation information with TFs/miRNAs at first column target the gene at second column (-t)
  
    5. fixed TFs and miRNAs regulation information with TFs/miRNAs at first column target the miRNAs/TF at second column(--F)

Output files:

  Results will be saved as:

    TFs feed forward loop: TF_FFL_output_core.txt

    miRNAs feed forward loop: miRNA_FFL_output_core.txt

    TFs miRNAs composite feed forward loop: TF_miRNAs_composite_FFL_output_core.txt

    Random permutation loop counts: FFL_random_count.txt

    Z score: z_score_out.txt

Here's an example, provided the files are in the same directory as this python script:

    1. For calculation of actural feed forward loop only
    ./python FFL.py -n network.txt 

    2. Include random permutation
    ./python FFL.py -n network.txt -g 10 -G gene.txt --tm fixed_TFs_miRNAs.txt -t TFs_miRNAs_to_gene.txt -F TFs_miRNAs_regulation.txt -s 1000
