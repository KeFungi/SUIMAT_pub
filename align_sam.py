from itertools import takewhile, dropwhile
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from random import randint
import pysam
import argparse
import sys

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", metavar='input.sam', help="input .sam/.bam file", type=str)
parser.add_argument("-r", metavar='reali.fasta', help="alignment of reference sequences in fasta format", type=str)
parser.add_argument("-o", metavar='output sam', help="output .sam file", type=str)

args = parser.parse_args()

def make_longCIGAR(CIGAR_list):
    """decompress CIGAR string to one character per posistion"""
    longCIGAR = [str(type) * count for (type, count) in CIGAR_list]
    flat_list = [item for sublist in longCIGAR for item in sublist]
    return flat_list


def compress(s):
    """compress longCIGAR string to numberXcharacter"""
    if not s:
        return []
    return [(int(s[0]), (len(list(takewhile(lambda c: c == s[0], s)))))] + compress(
        "".join(list(dropwhile(lambda c: c == s[0], s))))


def get_oldgap_pos(seq):
    """get gap positions in a sequence progressively"""
    pos=0
    gap_list = []
    for i in range(len(seq)):
        if seq[i] != '-':
            pos = pos + 1
        elif seq[i] == '-':
            gap_list = gap_list + [pos+1]
    return gap_list


def amb_convert(base, random=False):
    ambiguities = {"AA": "A", "TT": "T", "CC": "C", "GG": "G",
                   "AT": "W", "TA": "W", "AC": "M", "CA": "M",
                   "AG": "R", "GA": "R", "TC": "Y", "CT": "Y",
                   "TG": "K", "GT": "K", "CG": "S", "GC": "S",
                   "..": "N"}

    if base[0] == '-' and base[1] == '-':
        return '-'
    elif base[0] == '-':
        return base[1]
    elif base[1] == '-':
        return base[0]
    if random:
        return base[randint(0, 1)]
    else:
        return ambiguities[base]


def consensus2(seq1, seq2):
    if len(seq1) != len(seq2):
        print('different length in alignment')
        exit(1)
    else:
        convseq = []
        for i in range(len(seq1)):
            base = amb_convert(seq1[i]+seq2[i], random=True)
            convseq.append(base)
    return ''.join(convseq)

# read reference sequences
refseq = SeqIO.parse(args.r, 'fasta')
refali = AlignIO.read(args.r, 'fasta')

# consensus sequence
refcon = consensus2(refali[0], refali[1])

# read sam header
with open(args.i, 'r') as insam:
    sam = pysam.AlignmentFile(insam, 'r')
    refnames = list(sam.header.references) + ['consensus']
    lengths = list(sam.header.lengths) + [len(refcon)]

# create output files
outsam = args.o
outfasta = args.o.rstrip('sam|bam') + 'consensus.fasta'
merge_sam = pysam.AlignmentFile(outsam, "w", reference_names=refnames, reference_lengths=lengths)
new_ref_file = open(outfasta, 'w')

if args.o == args.i:
    print('cannot have the same output path as input path')
    sys.exit(1)

# all positions are 0-base
for ref in refseq:
    seq = ref.seq
    raw_gap_pos = get_oldgap_pos(seq)  # gap position in old_ref positions
    sam = pysam.AlignmentFile(args.i, 'r')  # read sam

    for seg in sam:
        if seg.reference_name == ref.name:
            longCIGAR = make_longCIGAR(seg.cigar)
            seg_start = seg.reference_start
            seg_end = seg.reference_end

            # remove insertions (cigar I) in sequence and cigar
            noD_CIGAR = [base for base in longCIGAR if base != '2']  # not consider cigar-D; it does not have accordant sequence base
            keep_pos = [i for i in (range(len(noD_CIGAR))) if noD_CIGAR[i]!='1']
            new_vseq = [seg.query_sequence[i] for i in keep_pos]
            new_vqua = [seg.qual[i] for i in keep_pos]
            noI_CIGAR = [base for base in longCIGAR if base != '1']

            # adjust positions according to gaps in the reference
            passed_gap = len([pos for pos in raw_gap_pos if pos < seg_start])  # gaps before seg; add to query_start
            considered_gap_pos = [pos for pos in raw_gap_pos if seg_start <= pos < seg_end]  # gaps to add in seg (in old_ref positions)
            seq_gap_pos = [pos - seg_start for pos in considered_gap_pos]  # gap in seg positions (indices)

            # add gap to seg sequence progressively
            added_gap = 0
            for pos in seq_gap_pos:
                noI_CIGAR.insert(pos+added_gap, '2')
                added_gap = added_gap + 1

            # create new seg object
            new_seg = pysam.AlignedSegment(merge_sam.header)

            new_seg.reference_start = seg.reference_start + passed_gap  # adjust to new start position
            new_seg.seq = ''.join(new_vseq)  # new sequence
            new_seg.qual = ''.join(new_vqua)  # new quality
            new_seg.cigar = compress(noI_CIGAR)  # new CIGAR
            new_seg.reference_name = "consensus"  # new reference name
            new_seg.query_name = seg.query_name
            new_seg.flag = seg.flag
            new_seg.mapping_quality = seg.mapping_quality
            new_seg.tags = seg.tags

            merge_sam.write(new_seg)

# write and close
SeqIO.write(SeqIO.SeqRecord(seq=Seq(refcon, IUPAC.ambiguous_dna), id='consensus', description=''), new_ref_file, 'fasta')
new_ref_file.close()
merge_sam.close()
