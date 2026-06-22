#!/usr/bin/env python3
"""对比 blastmini 和 blastn 的结果，生成性能与一致性图表。"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_tsv(filepath):
    """读取 TSV，假设有 header 且与 blastmini 输出兼容"""
    return pd.read_csv(filepath, sep='\t')

def main():
    blastmini_file = snakemake.input.blastmini_results   # noqa: F821
    blastn_file = snakemake.input.blastn_results         # noqa: F821
    plot_file = snakemake.output.plot                    # noqa: F821

    # 加载数据
    mini = load_tsv(blastmini_file)
    n = load_tsv(blastn_file)

    # 打印前几行查看
    print("Blastmini results:")
    print(mini.head())
    print("\nBlastn results:")
    print(n.head())

    # 简单对比：按 query_id 合并，比较得分或 identity
    # 这里假设有 query_id 和 subject_id 列，且 blastn 输出列名为 qseqid, sseqid, score, pident 等
    # 统一列名
    mini_cols = {'query_id': 'qseqid', 'subject_id': 'sseqid', 'score': 'score', 'identity_percent': 'pident'}
    n_cols = {'qseqid': 'qseqid', 'sseqid': 'sseqid', 'score': 'score', 'pident': 'pident'}

    mini_renamed = mini.rename(columns=mini_cols)[['qseqid', 'sseqid', 'score', 'pident']]
    n_renamed = n.rename(columns=n_cols)[['qseqid', 'sseqid', 'score', 'pident']]

    # 只保留每个 query 的最佳 hit（按 score 最高）
    mini_best = mini_renamed.loc[mini_renamed.groupby('qseqid')['score'].idxmax()]
    n_best = n_renamed.loc[n_renamed.groupby('qseqid')['score'].idxmax()]

    # 合并
    merged = pd.merge(mini_best, n_best, on='qseqid', suffixes=('_mini', '_n'))

    # 绘图：散点图比较得分
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax1 = axes[0]
    ax1.scatter(merged['score_mini'], merged['score_n'], alpha=0.7)
    ax1.plot([0, max(merged['score_mini'].max(), merged['score_n'].max())],
             [0, max(merged['score_mini'].max(), merged['score_n'].max())], 'r--')
    ax1.set_xlabel('Blastmini Score')
    ax1.set_ylabel('Blastn Score')
    ax1.set_title('Score Comparison')

    ax2 = axes[1]
    ax2.scatter(merged['pident_mini'], merged['pident_n'], alpha=0.7)
    ax2.plot([0, 100], [0, 100], 'r--')
    ax2.set_xlabel('Blastmini Identity (%)')
    ax2.set_ylabel('Blastn Identity (%)')
    ax2.set_title('Identity Comparison')

    plt.tight_layout()
    plt.savefig(plot_file, dpi=150)
    print(f"Benchmark plot saved to {plot_file}")

if __name__ == "__main__":
    main()