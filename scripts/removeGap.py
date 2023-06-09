from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def remove_gaps_from_multifasta(input_file, output_file):
    # Leer el archivo multifasta de entrada
    records = SeqIO.parse(input_file, "fasta")

    # Crear una lista para almacenar los registros de secuencia sin gaps
    gap_removed_records = []

    # Recorrer cada registro en el archivo multifasta
    for record in records:
        # Obtener la secuencia del registro y eliminar los gaps
        sequence = record.seq
        sequence = sequence.replace("-", "")

        # Crear un nuevo registro sin gaps
        gap_removed_record = SeqRecord(Seq(sequence), id=record.id, description="")

        # Agregar el registro sin gaps a la lista
        gap_removed_records.append(gap_removed_record)

    # Escribir los registros sin gaps en un nuevo archivo multifasta
    with open(output_file, "w") as handle:
        SeqIO.write(gap_removed_records, handle, "fasta")
