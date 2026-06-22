#!/usr/bin/env python3
"""随机从基因组中截取片段作为查询序列。"""

import random
import sys
from pathlib import Path
from blastmini.io import parse_fasta

def main():
    genome_file = snakemake.input.genome   # noqa: F821
    out_file = snakemake.output.queries    # noqa: F821
    n = snakemake.params.n                 # noqa: F821
    length = snakemake.params.length       # noqa: F821

    # 读取基因组序列
    records = list(parse_fasta(genome_file))
    if not records:
        raise ValueError("No genome sequences found")
    seq = records[0].sequence  # 假设只有一条序列
    seq_len = len(seq)

    if seq_len < length:
        raise ValueError(f"Genome length {seq_len} < query length {length}")

    # 随机截取
    with open(out_file, 'w') as f:
        for i in range(n):
            start = random.randint(0, seq_len - length)
            frag = seq[start:start+length]
            f.write(f">query_{i+1}\n{frag}\n")

if __name__ == "__main__":
    main()