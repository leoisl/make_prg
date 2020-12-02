from pathlib import Path
from multiprocessing import Pool
import os
import time
import logging
logging.basicConfig(
    format="%(asctime)s [%(levelname)s]: %(message)s",
    level=logging.INFO,
    filename=str(snakemake.log),
    filemode="w",
)


def infer_gene_from_msa(msa):
    MSA_path = Path(msa)
    MSA_path = MSA_path.with_suffix("")
    return MSA_path.name


def run_make_prg(msa):
    gene = infer_gene_from_msa(msa)

    start = time.time()

    logging.info(f"Running make_prg for {gene}")
    os.system(f"make_prg from_msa --max_nesting {max_nesting_lvl} --min_match_length {min_match_length} "
              f"--prefix {gene} {msa}")
    os.makedirs(output_folder, exist_ok=True)
    os.rename(f"{gene}.max_nest{max_nesting_lvl}.min_match{min_match_length}.prg", f"{output_folder}/{gene}.prg.fa")

    stop = time.time()
    runtime = stop - start
    logging.info(f"Finished running make_prg for {gene}")
    logging.info(f"make_prg runtime for {gene} in seconds: {runtime:.3f}")


MSAs = snakemake.input.MSAs
MSAs = [str(msa) for msa in MSAs]
processes = snakemake.threads
output_folder = snakemake.output.prgs
max_nesting_lvl = snakemake.params.max_nesting_lvl
min_match_length = snakemake.params.min_match_length

with Pool(processes=processes) as pool:
    pool.map(run_make_prg, MSAs)
