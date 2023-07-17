#!/usr/bin/env python
# coding: utf-8

# # Whole-Genome Data Analysis
# 
# ## Unraveling the Power of Genetic Sequences from Whole Genomes: A Pathway to Understanding Insect Detoxification
# 
# In this project, our main objectives are:
# 
# 1. Develop a function to efficiently extract genes of interest from one or more genomes.
# 2. Implement a function to filter sequences to the desired base pair (bp) length for targeted analysis.
# 3. Utilize MAFFT to align the retrieved sequences for accurate comparisons.
# 4. Generate and visualize gene trees to gain valuable insights.
# 
# Focus Gene: ATPase Alpha Subunit
# The ATPase alpha subunit plays a crucial role in the detoxification process employed by insects against toxic substances present in the plants they consume. Our analysis will encompass eighteen genomes of milkweed longhorn beetles, as well as other insect orders such as butterflies and true bugs. 
# 
# Through this Python-based genomic data analysis, we aim to enhance our understanding of these essential genetic sequences and their implications in insect detoxification mechanisms.

# ## Installing required software
# 
# BioPython is already installed so let's confirm we are in the bioinfo environment:

# In[9]:


get_ipython().system('echo $CONDA_DEFAULT_ENV')


# Load libraries and install iqtree and ete.

# In[271]:


from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from ete3 import Tree, TreeStyle, NodeStyle

# conda install -c bioconda iqtree ete3


# See content of current directory.

# In[33]:


get_ipython().system('ls -lh #$WD')


# Looks like we are in the right directory, now let's create a variable with the directory where the genome files are located:

# In[34]:


genomes_dir = 'functionally_annotated_genomes/'
get_ipython().system('ls -l $genomes_dir')


# ## 1. Extract gene secuences from genome files
# 
# We have 18 genome annotation files from where I want to extract only the sequences that correspond to ATPase subunit alpha. To avoid doing it one by one, let's create a function that uses as argument the name of the gene we are looking for and the directory where the genomes are located. The function will:
# 
# **a)** Extract sequences that match the desired gene.
# 
# **b)** Append sequences to a file.
# 
# **c)** Rename sequences to indicate what species they correspond to.

# In[94]:


def genefind(gene_name, genomes_dir, split_name=True, ignore_case=True):
  """
  Looks for genes one or several genomes. The sequences found are renamed using the file name where they were 
  found and are appended to a file. 
  param gene_name is name of the gene of interest.
  param genomes_dir is directory where the genome(s) are.
  param split_name looks for the name of the gene with its letters arranged in all order combinations.
  param ignore_case makes gene_name lowercase.
  """
  # Split the gene name because it can be in a different order (default)
  if split_name:
    terms = gene_name.split() # "ATPase subunit Alpha" --> ["ATPase","subunit", "Alpha"]
  else:
    terms = [gene_name]
  # Make all terms lowercase if ignoring case (default)
  if ignore_case:
    terms = [term.lower() for term in terms] # --> ["atpase","subunit", "alpha"]
  
  # Define variables
  genomes_files = 'genomes_files.list'
  get_ipython().system('ls $genomes_dir > $genomes_files')
  get_ipython().system('cat $genomes_files')

  gene_name_file = gene_name.replace(" ", "_") + ".fasta"

  # Clear the file so that it starts empty
  with open(gene_name_file, 'w') as f:
    pass

  with open(genomes_files, 'r') as f:
    for filename in f:
      # Remove space between names
      filename = filename.strip()
      print("\n" + filename)
      seqrecs = SeqIO.parse(genomes_dir + filename, 'fasta') 

      # Go over each record in this file
      for rec in seqrecs:
        # Create variable for the description of each record
        desc = rec.description
        if ignore_case:
          desc = desc.lower()
        # Extract sequences that match the gene
        if all([term in desc for term in terms]):
          # Append sequences to a file where all sequences will be stored
          print(rec.description)
          with open(gene_name_file, 'a') as g:   
            # Rename sequences to indicate what species they correspond to
            rec.id = filename.replace('.fasta', '_') + rec.id
            rec.description = ''
            g.write(rec.format('fasta'))
  # Count the number of sequences found                      
  number_sequences = get_ipython().getoutput('grep -c "^>" $gene_name_file')
  # Redefine data type of the variable (from string to integer)
  number_sequences = int(number_sequences[0])
  print(f"\n{number_sequences} sequences matching \"{gene_name}\" have been retrieved in the file \"{gene_name_file}\"")


# Let's use genefind() to look for the sequences that match the name of the gene ATPase subunit alpha. The function is designed to find sequences that match that name even when the words are not in the right order. For example, it will find sequences named "ATPase subunit alpha" as well as those named "ATPase alpha subunit" and "subunit alpha ATPase", etc.: 

# In[95]:


genefind("ATPase subunit alpha", "functionally_annotated_genomes/")


# Result: 123 sequences matching "ATPase subunit alpha" have been retrieved in the file "ATPase_subunit_alpha.fasta"

# ## 2. Filter the sequences to a desired bp length

# Using genefind(), 123 sequences were retrieved. The length of the sequences varies from 3460 to 50 bp. Some of the sequences are very short to be included in further steps of the pipeline. Before moving forward with sequence alignment, we will use the seqfilter() function to fliter sequences shorter than 1000 bp (the minimum length chosen is an arbitrary):

# In[96]:


def seqfilter(fasta_file, min_length):
  """
  Filter sequences with a user-defined threshold and creates two files, one containing sequences that have 
  passed the threshold and one with sequences that did not. 
  param fasta_file is name of the file.
  param min_length is the number of bp of the threshold.
  """
  # Create fasta files where to storage sequences 
  long_seqs_file = fasta_file.replace("fasta", "filtered.fasta")
  short_seqs_file = fasta_file.replace("fasta", "shortseqs.fasta")
    
  # Define seqrecs using SeqIO
  seqrecs = SeqIO.parse(fasta_file, "fasta") 
    
  # Clear the files so that they start empty
  with open(long_seqs_file, 'w') as f:
    pass
  with open(short_seqs_file, 'w') as f:
    pass     
  
  # Go over each record in the file
  for rec in seqrecs:
    # If the record is longer than the min_length
    if len(rec) > min_length:
      # Append sequences to a file where all long sequences will be stored
      print(rec.description)
      with open(long_seqs_file, 'a') as g:   
        g.write(rec.format('fasta'))
    else:
      # Append sequences to a file where all short sequences will be stored
      print(rec.description)
      with open(short_seqs_file, 'a') as g:   
        g.write(rec.format('fasta'))
  # Count the number of sequences found           
  num_long_sequences = get_ipython().getoutput('grep -c "^>" $long_seqs_file')
  # Redefine data type of the variable (from string to integer)
  num_long_sequences = int(num_long_sequences[0])
  print(f"\n{num_long_sequences} sequences passed the {min_length} bp threshold. They have been retrieved in the file \"{long_seqs_file}\"")          


# Let's use seqfilter() to filter the sequences in the "ATPase_subunit_alpha.fasta" file: 

# In[97]:


seqfilter("ATPase_subunit_alpha.fasta", 1000)


# Result: only 70 of the 123 sequences were longer than 1000 bp. Sequences of less than 1000 bp were stored in a different file in case they need to be inspected later.

# ## 3. Align the sequences retrieved using mafft

# Now that we have only sequences longer than 1000 bp, we'll proceed with aligning them.

# In[14]:


# Create variable for aligned sequences
filtered_atpases_file = "ATPase_subunit_alpha.filtered.fasta"
# Create a file to storage aligned sequences
aligned_file = filtered_atpases_file.replace(".fasta", ".mafft.fasta")

# Use mafft to align the sequences
get_ipython().system('mafft $filtered_atpases_file > $aligned_file')


# ## 4. Generate a gene tree
# 
# The next step is to use the aligned sequences to generate a gene tree. For this step we will assume that all recovered sequences are descendants of a common ancestor and we will use a maximum likelihood algorithm to obtain a gene tree.

# In[ ]:


# Use ModelFinderPlus (MFP) in IQTREE to determine the best-fit model:
get_ipython().system('iqtree -s $aligned_file -m MFP')


# What model best fitted our data?

# In[17]:


get_ipython().system('cat ATPase_subunit_alpha.filtered.mafft.fasta.iqtree')


# Looks like GTR+F+G4 is the best-fit model. Now let's look at the tree:

# In[26]:


get_ipython().system('cat ATPase_subunit_alpha.filtered.mafft.fasta.treefile.nwk')


# Looks good but not very useful to see the topology. The next thing to do is use ete to visualize the tree.

# ### Visualize the gene tree
# 
# We'll use ete to see the tree.

# In[31]:


# Create a variable for the gene tree generated with IQTREE
tree_file = "ATPase_subunit_alpha.filtered.mafft.fasta.treefile.nwk"

# Assign the tree to the variable iqtree
iqtree = Tree(tree_file)


# Let's take a look at the tree:

# In[32]:


print(iqtree)


# So far the tree is unrooted. Let's root the tree to one of the sequences of *Danaus plexippus*, the monarch butterfly:

# In[39]:


# Define a variable called Dplexippus

Dplexippus = "Dplexippus_DPOGS213623-TA"

# Define Dplexippus as the root of the tree
iqtree.set_outgroup(Dplexippus)
print(iqtree)


# ### Create a nice-looking figure

# We can get an idea of the topology by looking at the visualization created in the previous step. So far we see that some of the sequences are grouped by species, but there are others - especially in the genus *Tetraopes* (identifiable by a "T" before the species name) - where it would be interesting to take a closer look. To better visualize the topology we will color the tips (or leaves) according to their genus.
# 
# First let's print the leaf names:

# In[43]:


for leaf in iqtree:
    print(leaf.name)


# Now we need to assign each species to a subgroup. We will create a dictionary where each leaf.name will be a key and the genus name will be its value:

# In[254]:


# Create an empty dictionary called insect_group
insect_group = {}

# For each leave in the tree
for leaf in iqtree:
    # Evaluate its name
    if "Aglabripennis" in leaf.name:
    # and create a new set of key:value pair if it coitains the name of the genus
      insect_group[leaf.name] = 'Aglabripennis'
    elif "Tcastaneum" in leaf.name:
      insect_group[leaf.name] = 'Tcastaneum'
    elif "Dplexippus" in leaf.name:
      insect_group[leaf.name] = 'Dplexippus'
    elif "Ofasciatus" in leaf.name:
      insect_group[leaf.name] = 'Ofasciatus'
    elif "Tannulatus" in leaf.name:
      insect_group[leaf.name] = 'Tannulatus'
    elif "Tpilosus" in leaf.name:
      insect_group[leaf.name] = 'Tpilosus'
    elif "Tdiscoideus" in leaf.name:
      insect_group[leaf.name] = 'Tdiscoideus'
    elif "Tskillmani" in leaf.name:
      insect_group[leaf.name] = 'Tskillmani'
    elif "Tfemoratus" in leaf.name:
      insect_group[leaf.name] = 'Tfemoratus'
    elif "Ttetrophthalmus" in leaf.name:
      insect_group[leaf.name] = 'Ttetrophthalmus'
    elif "Clectularius" in leaf.name:
      insect_group[leaf.name] = 'Clectularius'
    else:
      break
# Print the key:value pairs  
for key, value in insect_group.items():
    print(key, ' : ', value)


# Result: It correctly assigned each leaf.name to the genus it corresponded.
# 
# BUT after having errors when trying to do:
# 
#     for leaf in iqtree:
#         leaf.insect_gropu = insect_group[leaf.name]
# 
# I tried a different approach: use the output to create a dictionary by hand.

# In[323]:


order = {
  'Dplexippus_DPOGS213623-TA'                                 :  "Dplexippus",
  'Aglabripennis_lcl|NW_019416398.1_cds_XP_018562049.1_6847'  :  "Aglabripennis",
  'Aglabripennis_lcl|NW_019416398.1_cds_XP_018562054.1_6850'  :  'Aglabripennis',
  'Aglabripennis_lcl|NW_019416398.1_cds_XP_018562052.1_6851'  :  "Aglabripennis",
  'Aglabripennis_lcl|NW_019416398.1_cds_XP_018562055.1_6853'  :  "Aglabripennis",
  'Aglabripennis_lcl|NW_019416398.1_cds_XP_018562050.1_6848'  :  "Aglabripennis",
  'Aglabripennis_lcl|NW_019416398.1_cds_XP_018562053.1_6852'  :  "Aglabripennis",
  'Aglabripennis_lcl|NW_019416398.1_cds_XP_018562056.1_6854'  :  "Aglabripennis",
  'Aglabripennis_lcl|NW_019416398.1_cds_XP_018562051.1_6849'  :  "Aglabripennis",
  'Aglabripennis_lcl|NW_019416398.1_cds_XP_018562057.1_6855'  :  "Aglabripennis",
  'Tfemoratus2_4327769.g89'                                   :  "Tfemoratus",
  'Tannulatus1_2465216.g23'                                   :  "Tannulatus",
  'Tpilosus2_4819112.g7'                                      :  "Tpilosus",
  'Tannulatus2_5311409.g40'                                   :  "Tannulatus",
  'Tpilosus1_4271307.g61'                                     :  "Tpilosus",
  'Ttetrophthalmus_maker-163393-exonerate_protein2genome-gene-0.5-mRNA-1'  :  "Ttetrophthalmus",
  'Tfemoratus1_3077924.g10'                                   :  "Tfemoratus",
  'Tfemoratus2_3114834.g47'                                   :  "Tfemoratus",
  'Tdiscoideus1_1854031.g58'                                  :  "Tdiscoideus",
  'Tdiscoideus2_2890578.g52'                                  :  "Tdiscoideus",
  'Tskillmani1_13332957.g20'                                  :  "Tskillmani",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_015837163.1_15311'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_008196417.1_15312'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_008196422.1_15315'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_015837164.1_15316'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_008196420.1_15317'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_008196423.1_15319'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_008196418.1_15313'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_008196421.1_15318'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_008196424.1_15320'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_008196419.1_15314'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007422.5_cds_XP_008196425.1_15321'       :  "Tcastaneum",
  'Dplexippus_DPOGS214640-TA'                                 :  "Dplexippus",
  'Clectularius_CLEC006776-RA'                                :  "Clectularius",
  'Ofasciatus_lcl|KK854249.1_cds_PTY13707.1_4820'             :  "Ofasciatus",
  'Ofasciatus_lcl|KK854243.1_cds_PTY13634.1_4747'             :  "Ofasciatus",
  'Ofasciatus_lcl|KK854243.1_cds_PTY13635.1_4748'             :  "Ofasciatus",
  'Ofasciatus_lcl|KK854243.1_cds_PTY13626.1_4739'             :  "Ofasciatus",
  'Ofasciatus_lcl|KK854697.1_cds_PTY17974.1_9087'             :  "Ofasciatus",
  'Ofasciatus_lcl|KK854787.1_cds_PTY18716.1_9829'             :  "Ofasciatus",
  'Aglabripennis_lcl|NW_019417389.1_cds_XP_023312013.1_17588' :  "Aglabripennis",
  'Tannulatus1_5852914.g76'                                   :  "Tannulatus",
  'Tannulatus2_6685946.g8'                                    :  "Tannulatus",
  'Tpilosus1_5418222.g21'                                     :  "Tpilosus",
  'Tpilosus2_7357420.g51'                                     :  "Tpilosus",
  'Tpilosus1_5497471.g120'                                    :  "Tpilosus",
  'Tdiscoideus1_6770603.g14'                                  :  "Tdiscoideus",
  'Tdiscoideus2_3746555.g55'                                  :  "Tdiscoideus",
  'Tskillmani1_12903966.g5'                                   :  "Tskillmani",
  'Tfemoratus1_9193525.g41'                                   : "Tfemoratus",
  'Tfemoratus2_584367.g107'                                   :  "Tfemoratus",
  'Ttetrophthalmus_maker-158804-exonerate_protein2genome-gene-0.0-mRNA-1'  :  "Ttetrophthalmus",
  'Tdiscoideus2_3689124.g70'                                  :  "Tdiscoideus",
  'Tcastaneum_lcl|NC_007423.3_cds_XP_971478.1_16258'          :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007423.3_cds_XP_015837921.1_16259'       :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007423.3_cds_XP_015837922.1_16260'       :  "Tcastaneum",
  'Dplexippus_DPOGS211186-TA'                                 :  "Dplexippus",
  'Tcastaneum_lcl|NC_007417.3_cds_XP_015840522.1_3417'        :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007417.3_cds_XP_008190275.2_3418'        :  "Tcastaneum",
  'Tcastaneum_lcl|NC_007417.3_cds_XP_972369.1_3419'           :  "Tcastaneum",
  'Aglabripennis_lcl|NW_019416354.1_cds_XP_018565491.1_5686'  :  "Aglabripennis",
  'Tannulatus2_6776066.g34'                                   :  "Tannulatus",
  'Tpilosus1_2065760.g33'                                     :  "Tpilosus",
  'Tpilosus2_4168272.g56'                                     :  "Tpilosus",
  'Tfemoratus1_9243440.g114'                                  :  "Tfemoratus",
  'Tfemoratus2_8651013.g277'                                  :  "Tfemoratus",
  'Ttetrophthalmus_maker-167874-exonerate_protein2genome-gene-0.0-mRNA-1'  :  "Ttetrophthalmus",
  'Tdiscoideus1_5565729.g65'                                  :  "Tdiscoideus",
  'Tskillmani1_13290908.g54'                                  :  "Tskillmani",
}


# I could not resolve the unhashable type: 'dict_keys' error. I convert it to a tuple but got a different error.

# In[324]:


for leaf in iqtree:
    leaf.order = order[leaf.name]


# Part of this script works fine but I do not know how to get ride of the last part:
#   
#   The genera of dict_keys(['Dplexippus_DPOGS213623-TA'...

# In[326]:


for leaf in iqtree:
    print(f"The genera of {leaf.name} is {leaf.order}")


# In[ ]:




