from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from collections import Counter

def countHaplotypes(alignment_file, output_file):
    alignment = AlignIO.read(alignment_file, "fasta")

    # We store the first sequence and exclude it from the analysis
    first_seq = alignment[0]
    alignment = alignment[1:]
    
    # We trim each sequence at position 500 before comparing
    seq_dict = Counter(str(record.seq) for record in alignment) # trim each seq at this pos due to primers [:5855]

    # We calculate the threshold that corresponds to 1% of the total number of sequences
    threshold = len(alignment) * 0.05

    # We create a new alignment with the first sequence
    new_alignment = MultipleSeqAlignment([SeqRecord(seq=first_seq.seq, id=first_seq.id, description=first_seq.description)])
    hap = 0
    for seq, count in seq_dict.items():
        # If the sequence appears less times than the threshold, we exclude it
        if count < threshold:
            continue
        hap += 1
        # We print the sequence and the number of times it appears
        #print(f'Sequence: {seq}, Number of times it appears: {count}')
        new_record = SeqRecord(seq=Seq(seq), id=f'{first_seq.id}_Hap-{hap}_({count})', description='')
        new_alignment.append(new_record)

    # We write the new alignment to the output file
    AlignIO.write(new_alignment, output_file, "fasta")
