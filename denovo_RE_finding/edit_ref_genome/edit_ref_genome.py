from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

ref="/public_data_resources/reference/GRCh38/GRCh38Decoy_no_alt.fa"

for seq_record in SeqIO.parse(ref,"fasta"):
    if seq_record.name=='chr17':
        print(seq_record.id)
        refSeq = str(seq_record.seq)
        old_chrom_length = str(len(str(refSeq)))
        refSeq=refSeq[:50199703]+'GGGGCCT'*int(sys.argv[1])+refSeq[50199703:]
        new_length = str(len(refSeq))
        new_seq=Seq(str(refSeq))
        new_record = SeqRecord(
            new_seq, 
            id=seq_record.id, 
            name=seq_record.name,
            description=seq_record.description.replace(old_chrom_length, new_length)
        )
    else:
        new_record=seq_record
    with open('mod.fa', "a") as f:
        SeqIO.write([new_record], f, "fasta")




# samtools fastq -1 test_1.fastq -2 test_2.fastq -s test_singletons.fastq /path/to/bam/DNA_01.bam
# bwa index mod.fa
# bwa mem mod.fa output_1.fastq output_2.fastq | samtools view -h -bS - > aligned_reads.bam
# samtools sort -o aligned_reads_std.bam aligned_reads.bam
# samtools index aligned_reads_std.bam

