#!/usr/bin/env python3
"""
blastmini 性能基准测试

循环测试不同 k 值（8, 11, 15）和不同数据库大小（400kb, 1Mb, 4.6Mb），
记录索引构建耗时、内存峰值、搜索耗时和命中数，并生成可视化图表。
"""

import os
import sys
import time
import random
import tempfile
import gc
from pathlib import Path

# 添加项目根目录到 sys.path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import matplotlib.pyplot as plt
import psutil

from blastmini.api import BlastMini
from blastmini.io import parse_fasta

# 配置
GENOME_FILE = Path("data/NC_000913.3.fa")
QUERY_LENGTH = 200
NUM_QUERIES = 30
K_VALUES = [8, 11, 15]
DB_SIZES = [400_000, 1_000_000, 4_641_652]
DB_LABELS = ['400kb', '1Mb', '4.6Mb']

random.seed(42)


def get_memory_usage():
    """返回当前进程的 RSS 内存占用（MB）"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)


def run_benchmark():
    """执行基准测试并生成结果"""
    if not GENOME_FILE.exists():
        raise FileNotFoundError(f"基因组文件不存在: {GENOME_FILE}")
    records = list(parse_fasta(GENOME_FILE))
    if not records:
        raise ValueError("基因组文件为空或无法解析")
    genome_seq = records[0].sequence
    genome_len = len(genome_seq)
    print(f"基因组长度: {genome_len:,} bp")

    # 生成固定查询集（从基因组随机截取）
    queries = []
    for _ in range(NUM_QUERIES):
        start = random.randint(0, genome_len - QUERY_LENGTH - 1)
        queries.append(genome_seq[start:start + QUERY_LENGTH])

    results = []

    for db_size, db_label in zip(DB_SIZES, DB_LABELS):
        # 创建临时数据库文件
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
            tmp.write(f">db_{db_label}\n{genome_seq[:db_size]}\n")
            db_file = tmp.name

        for k in K_VALUES:
            print(f"正在测试: 数据库={db_label}, k={k} ...")

            # 垃圾回收，减少干扰
            gc.collect()
            mem_before = get_memory_usage()

            # ----- 构建索引 -----
            build_start = time.time()
            blast = BlastMini.from_fasta(db_file, k=k, verbose=False)
            build_time = time.time() - build_start

            mem_after = get_memory_usage()
            memory_peak = max(mem_before, mem_after)

            # ----- 搜索查询 -----
            search_start = time.time()
            total_hits = 0
            for query in queries:
                result = blast.search(query, top_n=1)
                total_hits += result.num_hits
            search_time = time.time() - search_start

            # 记录结果
            results.append({
                'db_size': db_label,
                'db_bytes': db_size,
                'k': k,
                'build_time_sec': build_time,
                'memory_peak_mb': memory_peak,
                'search_time_sec': search_time,
                'total_hits': total_hits
            })

            del blast
            gc.collect()

        # 删除临时数据库文件
        os.unlink(db_file)

    # 转换为 DataFrame 并保存
    df = pd.DataFrame(results)
    df.to_csv('benchmark_results.csv', index=False)
    print("\n结果已保存到 benchmark_results.csv")
    print(df.to_string(index=False))

    # ---------- 绘图 ----------
    # 1. 搜索时间 vs k
    fig1, axes1 = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
    for i, db_label in enumerate(DB_LABELS):
        sub = df[df['db_size'] == db_label]
        ax = axes1[i]
        ax.plot(sub['k'], sub['search_time_sec'], marker='o', linestyle='-', linewidth=2)
        ax.set_title(f'数据库大小: {db_label}')
        ax.set_xlabel('k-mer 大小')
        ax.set_ylabel('搜索总耗时 (秒)')
        ax.grid(True, linestyle='--', alpha=0.7)
        for _, row in sub.iterrows():
            ax.annotate(f'{row["search_time_sec"]:.2f}s',
                        (row['k'], row['search_time_sec']),
                        textcoords="offset points", xytext=(0, 10), ha='center')
    plt.tight_layout()
    plt.savefig('benchmark_search_time.png', dpi=150)
    print("搜索耗时图已保存为 benchmark_search_time.png")

    # 2. 构建时间 vs k
    fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))
    for i, db_label in enumerate(DB_LABELS):
        sub = df[df['db_size'] == db_label]
        ax = axes2[i]
        ax.plot(sub['k'], sub['build_time_sec'], marker='s', linestyle='--', color='red')
        ax.set_title(f'数据库大小: {db_label}')
        ax.set_xlabel('k-mer 大小')
        ax.set_ylabel('索引构建耗时 (秒)')
        ax.grid(True, linestyle='--', alpha=0.7)
        for _, row in sub.iterrows():
            ax.annotate(f'{row["build_time_sec"]:.2f}s',
                        (row['k'], row['build_time_sec']),
                        textcoords="offset points", xytext=(0, 10), ha='center')
    plt.tight_layout()
    plt.savefig('benchmark_build_time.png', dpi=150)
    print("构建耗时图已保存为 benchmark_build_time.png")

    # 3. 内存峰值 vs k
    fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))
    for i, db_label in enumerate(DB_LABELS):
        sub = df[df['db_size'] == db_label]
        ax = axes3[i]
        ax.plot(sub['k'], sub['memory_peak_mb'], marker='^', linestyle=':', color='green')
        ax.set_title(f'数据库大小: {db_label}')
        ax.set_xlabel('k-mer 大小')
        ax.set_ylabel('内存峰值 (MB)')
        ax.grid(True, linestyle='--', alpha=0.7)
        for _, row in sub.iterrows():
            ax.annotate(f'{row["memory_peak_mb"]:.1f}MB',
                        (row['k'], row['memory_peak_mb']),
                        textcoords="offset points", xytext=(0, 10), ha='center')
    plt.tight_layout()
    plt.savefig('benchmark_memory.png', dpi=150)
    print("内存峰值图已保存为 benchmark_memory.png")


if __name__ == "__main__":
    run_benchmark()