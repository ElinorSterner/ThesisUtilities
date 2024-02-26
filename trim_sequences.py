#elinor used chatGPT to write this, and slightly modified it

import os
from Bio import SeqIO
from pathlib import Path

def trim_sequence(sequence, start, end):
    """
    Trim a nucleotide sequence from start to end positions.
    
    Args:
        sequence (str): The input nucleotide sequence.
        start (int): The start position for trimming (1-based indexing).
        end (int): The end position for trimming (1-based indexing).
    
    Returns:
        str: The trimmed nucleotide sequence.
    """
    return sequence[start-1:end]  # Adjusting for 0-based indexing

def trim_sequences_in_fasta(fasta_file, start, end):
    """
    Trim sequences in a FASTA file from start to end positions.
    
    Args:
        fasta_file (str): Path to the input FASTA file.
        start (int): The start position for trimming (1-based indexing).
        end (int): The end position for trimming (1-based indexing).
    
    Returns:
        list: List of tuples containing sequence ID and trimmed sequence.
    """
    trimmed_sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
       
        #Added this block 2/26/24 to skip trimming non- At sequences
        if record.id.startswith("At_"):
            trimmed_seq = trim_sequence(str(record.seq), start, end)
            trimmed_sequences.append((record.id, trimmed_seq))
        else:
            trimmed_seq = record.seq
            trimmed_sequences.append((record.id, trimmed_seq))
    return trimmed_sequences

# Folder containing FASTA files
folder_path = "/Users/elinorsterner/Documents/katzlab/allo/TrimmingGuidance/Renamed_Guidance"
start_position = 11
end_position = -10
Path(f'Trimmed_Renamed_Guidance').mkdir(parents=True, exist_ok=True)#makes output folder
outpath = "Trimmed_Renamed_Guidance"


# Iterate through each file in the folder

for filename in os.listdir(folder_path):
    if filename.endswith(".fasta") or filename.endswith(".fa"):
        file_path = os.path.join(folder_path, filename)
        trimmed_sequences = trim_sequences_in_fasta(file_path, start_position, end_position)
        print("Trimmed Sequences in", filename)

        with open(f'{outpath}/{filename}_trimmed.fasta', "w") as o:
            for seq_id, trimmed_seq in trimmed_sequences:
                o.write(f'>{seq_id}\n{trimmed_seq}\n')

