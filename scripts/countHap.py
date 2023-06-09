from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from collections import Counter

def countHaps(alignment_file, output_file):
    alignment = AlignIO.read(alignment_file, "fasta")

    # Almacenamos la primera secuencia y la excluimos del análisis
    first_seq = alignment[0]
    alignment = alignment[1:]
    # Comparamos las secuencias y contamos las secuencias idénticas
    seq_dict = Counter(str(record.seq) for record in alignment)

    # Calculamos el umbral que corresponde al 1% del total de secuencias
    threshold = len(alignment) * 0.01

    # Creamos un nuevo alineamiento con la primera secuencia
    new_alignment = MultipleSeqAlignment([first_seq])
    hap = 0
    for seq, count in seq_dict.items():
        # Si la secuencia aparece menos veces que el umbral, la excluimos
        if count < threshold:
            continue
        hap += 1
        # Imprimimos la secuencia y el número de veces que aparece
        #print(f'Secuencia: {seq}, Número de veces que aparece: {count}')
        new_record = SeqRecord(seq=Seq(seq), id=f'{first_seq.id}_Hap-{hap}_({count})', description='')
        new_alignment.append(new_record)

    # Escribimos el nuevo alineamiento en el archivo de salida
    AlignIO.write(new_alignment, output_file, "fasta")

