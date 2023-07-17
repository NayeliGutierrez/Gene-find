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
  !ls $genomes_dir > $genomes_files
  !cat $genomes_files

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
  number_sequences = !grep -c "^>" $gene_name_file
  # Redefine data type of the variable (from string to integer)
  number_sequences = int(number_sequences[0])
  print(f"\n{number_sequences} sequences matching \"{gene_name}\" have been retrieved in the file \"{gene_name_file}\"")
