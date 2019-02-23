# A first try at C++ make_prg

Graph-based approach, I'd say many things are wrong, but maybe can be built upon this (or not!)

## Build and compile:
`mkdir build && cd build && cmake .. && make`


## Running (will run the test cases):
`./make_prg`


## Output:
Last 5 tests fail:
```
===================================================
[TEST 1]: OK
File: ../test/nonmatch.fa
COMPUTED PRG:  5 AAACGTGGTT 6 CCCCCCCCCC 5
CORRECT  PRG:  5 AAACGTGGTT 6 CCCCCCCCCC 5
===================================================
[TEST 2]: OK
File: ../test/match.nonmatch.fa
COMPUTED PRG: AAACG 5 TGGTT 6 CCCCC 5
CORRECT  PRG: AAACG 5 TGGTT 6 CCCCC 5
===================================================
[TEST 3]: OK
File: ../test/nonmatch.match.fa
COMPUTED PRG:  5 AAACGT 6 CCCCCC 5 GGTT
CORRECT  PRG:  5 AAACGT 6 CCCCCC 5 GGTT
===================================================
[TEST 4]: OK
File: ../test/match.nonmatch.match.fa
COMPUTED PRG: AAACG 5 T 6 C 5 GGTT
CORRECT  PRG: AAACG 5 T 6 C 5 GGTT
===================================================
[TEST 5]: OK
File: ../test/shortmatch.nonmatch.match.fa
COMPUTED PRG:  5 AAACGT 6 ATTTTC 5 GGTT
CORRECT  PRG:  5 AAACGT 6 ATTTTC 5 GGTT
===================================================
[TEST 6]: OK
File: ../test/match.nonmatch.shortmatch.fa
COMPUTED PRG: AAAC 5 GTGGTT 6 CCCCCT 5
CORRECT  PRG: AAAC 5 GTGGTT 6 CCCCCT 5
===================================================
[TEST 7]: ***FAILED***
File: ../test/match.staggereddash.fa
COMPUTED PRG: AAACG 5 TG-- 6 --TG 5 GTT
CORRECT  PRG: AAACGTGGTT
===================================================
[TEST 8]: ***FAILED***
File: ../test/contains_n.fa
COMPUTED PRG: AAA 5 CGT 6 NGC 5 GGTT
CORRECT  PRG: AAACG 5 T 6 C 5 GGTT
===================================================
[TEST 9]: ***FAILED***
File: ../test/contains_RYKMSW.fa
COMPUTED PRG: AAACG 5 Y 6 C 5 GGTT
CORRECT  PRG: AAACG 5 T 6 C 5 GGTT
===================================================
[TEST 10]: ***FAILED***
File: ../test/contains_n_and_RYKMSW.fa
COMPUTED PRG:  5 AANCGY 6 AAACGC 5 GGTT
CORRECT  PRG: AAACG 5 T 6 C 5 GGTT
===================================================
[TEST 11]: ***FAILED***
File: ../test/contains_n_and_RYKMSW_no_variants.fa
COMPUTED PRG: AAA 5 CGT 6 TNC 5 GGTT
CORRECT  PRG: AAACGTGGTT
===================================================
```



# ORIGINAL README:



# make_prg
Code to create a PRG for input to Pandora (https://github.com/rmcolq/pandora) from a Multiple Sequence Alignment file.

__Requirements__
Expects you to have python3 and nextflow installed and in your path and a config file for nextflow set up if you are working on a cluster. Can run on python2.7+ if the command in nextflow file is edited.
Note that nextflow does not play nicely when files are in mounted or shared folders.

__Usage__

    Usage: nextflow run make_prg_nexflow.nf <arguments>
  
    Required arguments:
      --tsv_in  FILENAME  An index file of MSA to build PRGs of
      --pipeline_root DIRECTORY Absolute path to make_prg
    
    Optional arguments:
  
__Download__
```
git clone https://github.com/rmcolq/make_prg.git
cd make_prg
pip3 install -r requirements.txt
pytest 
```

__Input__
Multiple Sequence Alignment files for genes/dna sequences for which we want PRGs, and an tab-separated index of these in the form:
```
sample_id       infile
GC0000001   /absolute/path/to/GC0000001_na_aln.fa.gz
GC0000002   /absolute/path/to/GC0000002_na_aln.fa
```

__Changing parameters__

There are some parameters at the top of the nextflow file which could be changed but which I have not made command line parameters:
```
max_nesting             This is the maximum number depth of bubbles in PRG, setting to 1 will allow variants, \\ 
                        but no nesting
min_match_length        Controls graph complexity 
alignment_format        Any format accepted by biopython's AlignIO
max_forks_make_prg      If working on a cluster which allows unlimited parallel jobs per user, this will be \\
                        used by nextflow to control maximum number of processes of this type that can run in \\
                        parallel. 
max_forks_make_fasta   
```
