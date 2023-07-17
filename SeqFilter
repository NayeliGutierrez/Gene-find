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
  num_long_sequences = !grep -c "^>" $long_seqs_file
  # Redefine data type of the variable (from string to integer)
  num_long_sequences = int(num_long_sequences[0])
  print(f"\n{num_long_sequences} sequences passed the {min_length} bp threshold. They have been retrieved in the file \"{long_seqs_file}\"")          
