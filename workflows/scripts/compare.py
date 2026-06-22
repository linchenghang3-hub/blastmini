#!/usr/bin/env python3
"""
对比 blastmini 和 blastn 的搜索结果。

输入：
    - blastmini_results.tsv  (由 blastmini search --format tsv 生成)
    - blastn_results.tsv     (由 blastn -outfmt 6 生成，使用与 blastmini 类似的列)

输出：
    - 统计报告 (stdout 或文本文件)
    - 对比图 (PNG/PDF)
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
print("snakemake" in globals())

def load_blastmini_tsv(filepath):
    """加载 blastmini 的 TSV 输出（有 header）。"""
    df = pd.read_csv(filepath, sep='\t')
    # 确保列名统一为小写，便于处理
    df.columns = df.columns.str.lower()
    return df


def load_blastn_tsv(filepath):
    columns = [
        'qseqid', 'sseqid', 'bitscore', 'pident', 'length',
        'qstart', 'qend', 'sstart', 'send', 'evalue'
    ]
    df = pd.read_csv(filepath, sep='\t', header=None, names=columns)
    # 添加 score 列（使用 bitscore 的值）以便统一处理
    df['score'] = df['bitscore']
    return df


def get_best_hit_per_query(df, score_col='score'):
    """
    对于每个查询，选取得分最高的 hit（如果得分相同，则选 E-value 最小的）。
    适用于 blastmini 和 blastn。
    """
    # 如果 df 中有 evalue 列，优先用 evalue 作为次要排序
    if 'evalue' in df.columns:
        df_sorted = df.sort_values([score_col, 'evalue'], ascending=[False, True])
    else:
        df_sorted = df.sort_values(score_col, ascending=False)
    return df_sorted.drop_duplicates(subset=['qseqid'], keep='first')


def compare_results(mini_df, blastn_df):
    """
    合并两个数据框，按 qseqid 进行内连接，比较得分和 identity。
    返回合并后的 DataFrame 和统计字典。
    """
    # 重命名列以避免冲突
    mini_renamed = mini_df.rename(columns={
        'score': 'score_mini',
        'pident': 'pident_mini',
        'evalue': 'evalue_mini',
        'qseqid': 'qseqid',
        'sseqid': 'sseqid_mini'
    })
    blastn_renamed = blastn_df.rename(columns={
        'score': 'score_blastn',
        'pident': 'pident_blastn',
        'evalue': 'evalue_blastn',
        'sseqid': 'sseqid_blastn'
    })

    # 合并（只保留两者都匹配到的查询）
    merged = pd.merge(mini_renamed, blastn_renamed, on='qseqid', how='inner')

    # 计算差异
    merged['score_diff'] = merged['score_blastn'] - merged['score_mini']
    merged['pident_diff'] = merged['pident_blastn'] - merged['pident_mini']

    # 统计信息
    stats = {
        'total_queries': len(merged),
        'score_correlation': merged['score_mini'].corr(merged['score_blastn']),
        'pident_correlation': merged['pident_mini'].corr(merged['pident_blastn']),
        'mean_score_diff': merged['score_diff'].mean(),
        'std_score_diff': merged['score_diff'].std(),
        'mean_pident_diff': merged['pident_diff'].mean(),
        'std_pident_diff': merged['pident_diff'].std(),
        'same_subject_ratio': (merged['sseqid_mini'] == merged['sseqid_blastn']).mean()
    }
    return merged, stats


def plot_comparison(merged, output_plot):
    """生成对比散点图。"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # 得分散点图
    ax = axes[0]
    ax.scatter(merged['score_mini'], merged['score_blastn'], alpha=0.6, edgecolors='w', s=60)
    max_score = max(merged['score_mini'].max(), merged['score_blastn'].max())
    ax.plot([0, max_score], [0, max_score], 'r--', lw=2, label='Identity Line')
    ax.set_xlabel('Blastmini Score')
    ax.set_ylabel('Blastn Score')
    ax.set_title('Score Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 身份百分比散点图
    ax = axes[1]
    ax.scatter(merged['pident_mini'], merged['pident_blastn'], alpha=0.6, edgecolors='w', s=60)
    ax.plot([0, 100], [0, 100], 'r--', lw=2, label='Identity Line')
    ax.set_xlabel('Blastmini Identity (%)')
    ax.set_ylabel('Blastn Identity (%)')
    ax.set_title('Identity Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_plot, dpi=150)
    print(f"Plot saved to {output_plot}")


def generate_report(stats, output_txt=None):
    """生成文本报告。"""
    lines = []
    lines.append("=" * 50)
    lines.append("COMPARISON REPORT: blastmini vs blastn")
    lines.append("=" * 50)
    lines.append(f"Number of common queries: {stats['total_queries']}")
    lines.append(f"Score correlation: {stats['score_correlation']:.4f}")
    lines.append(f"Identity correlation: {stats['pident_correlation']:.4f}")
    lines.append(f"Mean score difference (blastn - mini): {stats['mean_score_diff']:.2f}")
    lines.append(f"Std score difference: {stats['std_score_diff']:.2f}")
    lines.append(f"Mean identity difference (blastn - mini): {stats['mean_pident_diff']:.2f}%")
    lines.append(f"Std identity difference: {stats['std_pident_diff']:.2f}%")
    lines.append(f"Same subject ID ratio: {stats['same_subject_ratio']:.2%}")
    lines.append("=" * 50)

    report = "\n".join(lines)
    if output_txt:
        with open(output_txt, 'w') as f:
            f.write(report)
        print(f"Report saved to {output_txt}")
    else:
        print(report)


def main():
    # 如果直接在命令行运行，可以使用参数
    import argparse
    parser = argparse.ArgumentParser(description='Compare blastmini and blastn results')
    parser.add_argument('--mini', required=True, help='blastmini results TSV')
    parser.add_argument('--blastn', required=True, help='blastn results TSV')
    parser.add_argument('--plot', default='comparison.png', help='Output plot file')
    parser.add_argument('--report', help='Output report text file')
    args = parser.parse_args()

    # 如果通过 Snakemake 调用，会使用 snakemake.input / output
    # 这里兼容两种方式
    if 'snakemake' in globals():
        mini_file = snakemake.input.mini
        blastn_file = snakemake.input.blastn
        plot_file = snakemake.output.plot
        report_file = snakemake.output.get('report', None)
    else:
        mini_file = args.mini
        blastn_file = args.blastn
        plot_file = args.plot
        report_file = args.report

    # 加载数据
    mini_df = load_blastmini_tsv(mini_file)
    blastn_df = load_blastn_tsv(blastn_file)

    # 统一列名：blastmini 可能叫 query_id, subject_id, score, identity_percent
    # 转换为标准：qseqid, sseqid, score, pident
    mini_std = mini_df.rename(columns={
    'query_id': 'qseqid',
    'subject_id': 'sseqid',
    'score': 'score',
    'identity_percent': 'pident'
    })
    # 如果还有其他列，保留
    blastn_std = blastn_df.rename(columns={
        'pident': 'pident'  # 已有
    })

    # 选择最佳 hit (score 最高)
    mini_best = get_best_hit_per_query(mini_std, score_col='score')
    blastn_best = get_best_hit_per_query(blastn_std, score_col='score')

    # 合并对比
    merged, stats = compare_results(mini_best, blastn_best)

    # 生成报告
    if report_file:
        generate_report(stats, report_file)
    else:
        generate_report(stats)

    # 生成图表
    if merged.shape[0] > 1:
        plot_comparison(merged, plot_file)
    else:
        print("Not enough data points for plot.")


if __name__ == "__main__":
    main()