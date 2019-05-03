from make_prg_from_msa import  *
from Bio import AlignIO
import sys

def test_get_interval_seqs(inputFile):
    print(f"Testing {inputFile}")
    seqs = AlignIO.read(inputFile, "fasta")
    print(f"Input:{seqs}\nOutput:{get_interval_seqs(seqs)}")
    print()


for i in range(6):
    try:
        test_get_interval_seqs(f"test/get_interval_seqs_test{i}.fa")
    except:
        print(sys.exc_info()[0])