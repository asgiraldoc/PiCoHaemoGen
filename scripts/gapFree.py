from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

def remove_gap_columns(input_alignment_file, output_alignment_file, trim_at_position):
    # Parse the input alignment file
    alignments = list(SeqIO.parse(input_alignment_file, "fasta"))

    # Get the number of sequences and the alignment length
    num_sequences = len(alignments)
    alignment_length = len(alignments[0].seq)

    # Set to store the positions with gaps
    gap_positions = set()

    # Identify the positions with gaps in any sequence
    for i in range(alignment_length):
        has_gap = False
        for j in range(num_sequences):
            # Check if the current position has a gap in any sequence
            if alignments[j].seq[i] == "-":
                has_gap = True
                break
        if has_gap:
            # Add the position to the gap positions set
            gap_positions.add(i)

    # Create a new alignment without the gap columns
    new_alignments = []
    for alignment in alignments:
        # Generate the sequence without the gap positions
        seq = "".join(
            [alignment.seq[i] for i in range(alignment_length) if i not in gap_positions]
        )
        if trim_at_position == "yes":
            # Trim the sequence at position 5866 (adjust as needed)
            trimmed_seq = seq[:5866]

            # Create a SeqRecord for the new alignment sequence
            new_alignment = SeqRecord(Seq(trimmed_seq), id=alignment.id, description="")

            # Add the new alignment sequence to the list
            new_alignments.append(new_alignment)
        else:
            new_alignment = SeqRecord(Seq(seq), id=alignment.id, description="")
            new_alignments.append(new_alignment)

    # Write the new alignment to the output file in FASTA format
    SeqIO.write(new_alignments, output_alignment_file, "fasta")




# filesG0  = [f for f in os.listdir() if f.endswith('_aliHap.fasta')]
# for file in filesG0:
#     outf = str(file).split("-")[0] + "_aliHap-ng.fasta"
#     remove_gap_columns(file, outf, "no")