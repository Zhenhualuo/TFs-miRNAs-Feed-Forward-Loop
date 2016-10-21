import optparse
import random
import sys
import numpy as np
from collections import defaultdict

# ====================================================================================================================================
def print_usage(option, opt, value, parser):


    usage_message = """

# ----------------------------------------------------------------------------------------------------------------------------

This program will calculate and export Transcription factors (TFs) feed forward loop, miRNAs 
feed forward loop and TFs-miRNAs composite feed forward loop, random permutation.

* Required input:
All files must be a tab-separated table.

For calculation of actural feed forward loop only:
  1.network file with TFs/miRNAs at first column, target the gene at second column
  
For calculation of actural feed forward loop and random permutation, you MUST provide all following files:
  1.network file with TFs/miRNAs at first column target and gene at second column (-n)
  2.gene list for random selection (-g)
  3.fixed TFs and miRNAs list (-l)
  4.TFs/miRNAs and all human coding gene regulation information with TFs/miRNAs at first column, targeted gene at second column (-t)
  5.fixed TFs and miRNAs regulation information with TFs/miRNAs at first column target the miRNAs/TF at second column(-i)
  
Here's an example, provided the files are in the same directory as this python script:

1. For calculation of actural feed forward loop only
./python FFL.py -n network.txt 

2. Include random permutation
./python FFL.py -n network.txt -m 10 -g gene.txt -l fixed_TFs_miRNAs.txt -t TFs_miRNAs_to_gene.txt -i TFs_miRNAs_regulation.txt -s 1000

# -----------------------------------------------------------------------------------------------------------------------------

    """

    print usage_message

    sys.exit()

# ====================================================================================================================================  
    
if __name__ == '__main__':

    parser = optparse.OptionParser()

    parser.add_option('-u', '--usage',
                      help    ='show more info on how to use this script',
                      action="callback", callback=print_usage)

    parser.add_option('-n', 
                      help    ='regulatory network file with TFs/miRNAs at first column and target gene at second column',
                      dest    = 'network_file',
                      default ='none',
                      type    = "string")
    
    parser.add_option('-m', 
                      help    ='number of genes for random selection',
                      dest    ='m',
                      default ='0',
                      type    = "int")
    
    parser.add_option('-g', 
                      help    ='gene list for random selection',
                      dest    ='g',
                      default ='none',
                      type    = "string")    

    parser.add_option('-l', 
                      help    ='file containing TFs and miRNAs',
                      dest    ='l',
                      default ='none',
                      type    = "string")

    parser.add_option('-t', 
                      help    ='TFs and miRNAs target all coding gene file with TFs/miRNAs at first column target gene at second column',
                      dest    ='TF_MI_gene',
                      default ='none',
                      type    = "string")
    
    parser.add_option('-i', 
                      help    ='file containing fixed TF-miRNA interaction',
                      dest    ='TF_MI',
                      default ='none',
                      type    = "string")     

    parser.add_option('-s', 
                      help    ='number of random simulations [default: 1000]',
                      dest    ='sims',
                      default ='1000',
                      type    = "int")

    
    (opts, args) = parser.parse_args()
    network_file = opts.network_file
    sims         = opts.sims
    l            = opts.l
    m            = opts.m
    g            = opts.g
    TF_MI_gene   = opts.TF_MI_gene
    TF_MI        = opts.TF_MI
    
    # checking for input:
    if network_file == 'none':
        print """
        ERROR: you must specify an network file, for example:
        ./python FFL.py -n network.txt
    
        For more information, type
        ./python FFL.py -u
            
        """
        sys.exit(0)    
    
print "> Calculating TF feed forward loop (TF_FFL), miRNAs feed forward loop (miRNAs_FFL) and TF_miRNAs composited loop (TF_miRNA_FFL)..."
print "> Actual network results:"

TF_miRNAs_core_network = network_file #network file

# load the data as { "A": set(["B", "C", "D", ... ]) }
data_TF_miRNAs_core = defaultdict(set)
with open(TF_miRNAs_core_network) as TF_miRNAs_core:
    for line_core in TF_miRNAs_core:
        a,b = line_core.rstrip().split("\t") #a and b are genes
        if a != b:          # no self-loops
            data_TF_miRNAs_core[a].add(b)  


# find all triplets such that A -> B -> C and  A -> C, for TF_FFL and miRNA_FFL only
found_TF_FFL = []
found_miRNA_FFL = []
for a,bs in data_TF_miRNAs_core.items():
    bint = bs.intersection
    for b in bs:
        if a not in data_TF_miRNAs_core[b]:
            for c in bint(data_TF_miRNAs_core[b]):
                if a not in data_TF_miRNAs_core[c] and b not in data_TF_miRNAs_core[c]:
                    if '-miR-' in a or '-let-' in a:
                        found_miRNA_FFL.append("{},{},{}".format(a, b, c)) #this return miR_FFL
                    if '-miR-' in b or '-let-' in b:
                        found_TF_FFL.append("{},{},{}".format(a, b, c)) #this return TF_FFL
                    
print '> Number of TF_FFL is ' + str(len(found_TF_FFL)) # this is the number of TF_FFL 
print '> Number of miRNA_FFL is ' + str(len(found_miRNA_FFL)) # this is the number of miR_FFL

#Export TF_FFL and miRNA_FFL
found_TF_out_core=open(r'TF_FFL_output_core.txt','w')
found_TF_out_core.write('A are TFs, B are miRNAs, C are core genes. A target B, and A,B target C'+'\n')
for item in range(len(found_TF_FFL)):
    found_TF_out_core.write(found_TF_FFL[item]+'\n')

found_miRNA_out_core=open(r'miRNA_FFL_output_core.txt','w')
found_miRNA_out_core.write('A are miRNAs, B are TFs, C are core genes. A target B, and A,B target C'+'\n')
for item in range(len(found_miRNA_FFL)):
    found_miRNA_out_core.write(found_miRNA_FFL[item]+'\n')

# find all triplets such that A -> B -> C and  A -> C, for TF_miRNA_composite_FFL only
found_TF_miRNAs_FFL = []
for a,bs in data_TF_miRNAs_core.items():
    bint = bs.intersection
    for b in bs:
        if a in data_TF_miRNAs_core[b]:
            for c in bint(data_TF_miRNAs_core[b]):
                if a not in data_TF_miRNAs_core[c] and b not in data_TF_miRNAs_core[c]:
                    if '-miR-' not in c and '-let-' not in c and '-miR-' in b or '-let-' in b:
                        found_TF_miRNAs_FFL.append("{},{},{}".format(a, b, c))

print '> Number of TF_miRNA_composite_FFL is ' + str(len(found_TF_miRNAs_FFL)) 

#Export TFs_miRNAs_FFL
found_TF_miRNAs_out_core=open(r'TF_miRNAs_composite_FFL_output_core.txt','w')
found_TF_miRNAs_out_core.write('A are TFs, B are miRNAs, C are core genes. A target B, B target A and A,B target C'+'\n')
for item in range(len(found_TF_miRNAs_FFL)):
    found_TF_miRNAs_out_core.write(found_TF_miRNAs_FFL[item]+'\n')


found_TF_out_core.close()
found_miRNA_out_core.close()
found_TF_miRNAs_out_core.close()

FFL_actural=open(r'FFL_actural_count.txt','w')
FFL_actural.write('Number of TF_FFL is ' + str(len(found_TF_FFL))+'\n'+'Number of miRNA_FFL is '+str(len(found_miRNA_FFL))+'\n'+'Number of TF_miRNA_composite_FFL is '+str(len(found_TF_miRNAs_FFL)))
FFL_actural.close()

if sims == 'none' or m == '0' or g == 'none' or l == 'none' or TF_MI_gene == 'none' or TF_MI == 'none':
    print '> Results has been saved as following files:'
    print '> TFs feed forward loop: TF_FFL_output_core.txt'
    print '> miRNAs feed forward loop: miRNA_FFL_output_core.txt'
    print '> TFs miRNAs composite feed forward loop: TF_miRNAs_composite_FFL_output_core.txt'
    print '> Finished!'
    pass
else:
#random permutation
    print "> Start Random Permutation...."
    gene_list=open(g,'r').read().split('\r') #Biliary diease genes
    Fix_TF_miRNAs=open(l,'r').read().split('\r') # fixed Tf and miRNAs list
    TF_miRNA_gene_regulation=open(TF_MI_gene,'r').read().split('\r')
    TF_miRNA_regulation=open(TF_MI,'r').read().split('\r') #This file contains TF-miRNA and miRNA-TF regulation

    count_TF_FFL_random=[]
    count_miRNA_FFL_random=[]
    count_TF_miRNAs_FFL_random=[]
    count_random=[]
    q=0
    while q<sims:
        print '> ', q+1, 'results:'
        gene_random=random.sample(gene_list, g) #change N to number of genes, randomly select genes from human coding genes
        TF_miRNA_gene_regulation_select=[]
        for j in range(len(TF_miRNA_gene_regulation)):
            k=TF_miRNA_gene_regulation[j].strip().split("\t") #split TF_gene regulation table
            if k[0] in Fix_TF_miRNAs:     #fix TF
                if k[1] in gene_random: # select gene 
                    TF_miRNA_gene_regulation_select.append(TF_miRNA_gene_regulation[j]) #select regulation base on new gene and fix TF
            
#Merge TF_gene_regulation_select, miRNA_gene_regulation_select and TF_miRNA_regulation(FIXed)
        TF_miRNA_gene_loop=TF_miRNA_gene_regulation_select+TF_miRNA_regulation #combine TF-gene, miRNAs-gene, TF-miRNA and miRNA-TF regulation
        data = defaultdict(set)
        for line in TF_miRNA_gene_loop:
            if line!="":
                a,b=line.split("\t")
                if a != b:          # no self-loops
                    data[a].add(b)


# find all triplets such that A -> B -> C and  A -> C, for TF_FFL and miRNA_FFL only
        found_TF_FFL_random = []
        found_miRNA_FFL_random = []
        for a,bs in data.items():
            bint = bs.intersection
            for b in bs:
                if a not in data[b]:
                    for c in bint(data[b]):
                        if a not in data[c] and b not in data[c]:
                            if '-miR-' in a or '-let-' in a:
                                found_miRNA_FFL_random.append("{} {} {}".format(a, b, c)) #this return miR_FFL
                            if '-miR-' in b or '-let-' in b:
                                found_TF_FFL_random.append("{} {} {}".format(a, b, c)) #this return TF_FFL

        print '> Number of TF_FFL_random is ' + str(len(found_TF_FFL_random)) # this is the number of TF_FFL 
        print '> Number of miRNA_FFL_random is ' + str(len(found_miRNA_FFL_random)) # this is the number of miR_FFL

# find all triplets such that A -> B -> C and  A -> C, for TF_miRNA_composite_FFL only
        found_TF_miRNAs_FFL_random = []
        for a,bs in data.items():
            bint = bs.intersection
            for b in bs:
                if a in data[b]:
                    for c in bint(data[b]):
                        if a not in data[c] and b not in data[c]:
                            if '-miR-' not in c and '-let-' not in c and '-miR-' in b or '-let-' in b:
                                found_TF_miRNAs_FFL_random.append("{} {} {}".format(a, b, c))

        print '> Number of TF_miRNA_composite_FFL_random is ' + str(len(found_TF_miRNAs_FFL_random)) 
        count_TF_FFL_random.append(int(len(found_TF_FFL_random)))
        count_miRNA_FFL_random.append(int(len(found_miRNA_FFL_random)))
        count_TF_miRNAs_FFL_random.append(int(len(found_TF_miRNAs_FFL_random)))
        count_random.append(str(q+1)+'\n'+'TF_FFL'+','+str(len(found_TF_FFL_random))+'\n'+'miRNA_FFL'+','+str(len(found_miRNA_FFL_random))+'\n'+'TF_miRNA_composite_FFL'+','+ str(len(found_TF_miRNAs_FFL_random)))
        q=q+1
        
    count_random_out=open(r'FFL_random_count.txt','w')
    for n in range(len(count_random)):
        count_random_out.write(count_random[n]+'\n')
    count_random_out.close()
    
    #Caculate Z score
    count_TF_FFL_random_mean=np.mean(count_TF_FFL_random)
    count_TF_FFL_random_std=np.std(count_TF_FFL_random)
    count_TF_FFL_random_z_score = (1.*len(found_TF_FFL) - count_TF_FFL_random_mean)/count_TF_FFL_random_std
 
    count_miRNA_FFL_random_mean=np.mean(count_miRNA_FFL_random)
    count_miRNA_FFL_random_std=np.std(count_miRNA_FFL_random)
    count_miRNA_FFL_random_z_score = (1.*len(found_miRNA_FFL) - count_miRNA_FFL_random_mean)/count_miRNA_FFL_random_std

    count_TF_miRNAs_FFL_random_mean=np.mean(count_TF_miRNAs_FFL_random)
    count_TF_miRNAs_FFL_random_std=np.std(count_TF_miRNAs_FFL_random)
    count_TF_miRNAs_FFL_random_z_score = (1.*len(found_TF_miRNAs_FFL) - count_TF_miRNAs_FFL_random_mean)/count_TF_miRNAs_FFL_random_std

    z_score=open(r'z_score_out.txt','w')
    z_score.write('Z score of TF_FFL is ' + str(count_TF_FFL_random_z_score)+'\n'+'Z score of miRNA_FFL is '+str(count_miRNA_FFL_random_z_score)+'\n'+'Z score of TF_miRNA_composite_FFL is '+str(count_TF_miRNAs_FFL_random_z_score))
    z_score.close()
    
    print '> Results have been saved as following files:'
    print '> TFs feed forward loop: TF_FFL_output_core.txt'
    print '> miRNAs feed forward loop: miRNA_FFL_output_core.txt'
    print '> TFs miRNAs composite feed forward loop: TF_miRNAs_composite_FFL_output_core.txt'
    print '> Actural loop counts: FFL_actural_count.txt'
    print '> Random permutation loop counts: FFL_random_count.txt'
    print '> z score: z_score_out.txt'
    print '> Finished!'
