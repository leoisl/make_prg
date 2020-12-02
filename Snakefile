from snakemake.utils import min_version
from glob import glob
from pathlib import Path
import fileinput

min_version("5.30")

configfile: "config.yaml"


def get_fasta_files_in_dir(directory):
    fasta_files = []
    for fasta_file in glob(f"{directory}/*.fa"):
        fasta_files.append(str(Path(fasta_file).resolve()))
    return fasta_files


rule all:
    input:
        aggregated_prg = config["output_dir"] + "/prg.fa"


rule run_make_prg:
    input:
        MSAs = lambda wildcards: get_fasta_files_in_dir(directory = config["msas_dir"])
    output:
        prgs = directory(config["output_dir"] + "/prgs/")
    threads: int(config["threads"])
    resources:
        mem_mb = lambda wildcards, attempt: {1: 20000, 2: 50000, 3: 100000}.get(attempt, 150000)
    params:
        max_nesting_lvl = config.get("max_nesting_lvl", 5),
        min_match_length = config.get("min_match_length", 7),
    singularity: config["containers"]["make_prg"]
    shadow: "shallow"
    log:
        "logs/run_make_prg.log"
    script:
        "scripts/run_make_prg.py"



def is_header(line):
    return line.startswith(">")


def get_PRG_sequence(line):
    prg_sequence = line.rstrip()
    line_ends_digit = prg_sequence[-1].isdigit()
    if line_ends_digit:
        prg_sequence += " "
    return prg_sequence


def concatenate_several_prgs_into_one(input_prgs, output_prg):
    with open(output_prg, "w") as fout, fileinput.input(input_prgs) as fin:
        for line in fin:
            if is_header(line):
                fout.write(line)
            else:
                prg_sequence = get_PRG_sequence(line)
                fout.write(prg_sequence + "\n")


rule aggregate_prgs:
    input:
        prgs = rules.run_make_prg.output.prgs
    output:
        aggregated_prg = config["output_dir"] + "/prg.fa"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/aggregate_prgs.log"
    run:
        prgs = get_fasta_files_in_dir(input.prgs)
        concatenate_several_prgs_into_one(prgs, output.aggregated_prg)

