
---

```markdown
# 🧬 blastmini — 轻量级 BLAST 教学实现

[![CI](https://github.com/linchenghang3-hub/blastmini/actions/workflows/ci.yml/badge.svg)](https://github.com/linchenghang3-hub/blastmini/actions/workflows/ci.yml)
[![Docs](https://github.com/linchenghang3-hub/blastmini/actions/workflows/docs.yml/badge.svg)](https://github.com/linchenghang3-hub/blastmini/actions/workflows/docs.yml)
[![Docker](https://img.shields.io/badge/docker-available-blue)](https://hub.docker.com/)

**blastmini** 是一个面向教学与科研原型的轻量级序列搜索工具，从零实现了 BLAST 算法的核心启发式策略——**seed‑and‑extend**。本项目以 Python 包形式开发，旨在帮助开发者深入理解大规模序列比对中“速度‑灵敏度”权衡的工程本质。

---

## ✨ 主要功能

-  **种子搜索**：基于哈希表的 k‑mer 索引，快速定位候选匹配区域。
-  **X‑dropoff 双向延伸**：模拟 BLAST 的 ungapped 延伸，支持正负链。
-  **统计显著性估计**：通过置换检验估算 E‑value 和 p‑value。
-  **命令行接口**：统一 CLI，支持 `build`、`search`、`stats`、`view` 等子命令。
-  **Python API**：提供 `BlastMini` 高级接口，方便在脚本或 Jupyter 中调用。
-  **可复现工作流**：集成 Snakemake，支持与真实 BLASTn 对比并生成报告。
-  **容器化**：提供 Dockerfile，一键构建环境，消除依赖烦恼。
-  **自动文档**：Sphinx 构建，API 文档自动从代码注释生成。
-  **持续集成**：GitHub Actions 自动运行 pytest + flake8，保证代码质量。

---

## 📦 安装

### 方式一：使用 pip 安装（开发模式）

```bash
git clone https://github.com/linchenghang3-hub/blastmini.git
cd blastmini
pip install -e .
```

### 方式二：使用 Conda 创建完整开发环境

```bash
conda env create -f environment.yml
conda activate blastmini-dev
```

### 方式三：使用 Docker（无需安装 Python 环境）

```bash
docker build -t blastmini:latest .
docker run --rm blastmini:latest blastmini --help
```

---

## 🚀 快速上手

### 1. 准备数据

将你的基因组序列（FASTA 格式）放在 `data/` 目录下，例如 `data/NC_000913.3.fa`。

### 2. 构建索引

```bash
blastmini build -i data/NC_000913.3.fa -o data/index.json -k 11
```

### 3. 搜索查询

```bash
blastmini search -q query.fa -d data/index.json --database-fasta data/NC_000913.3.fa -o results.tsv --format tsv --top 10
```

### 4. 查看结果

```bash
blastmini view -r results.tsv --show-alignment
```

### 5. Python API 示例

```python
from blastmini.api import BlastMini

blast = BlastMini.from_fasta("data/NC_000913.3.fa", k=11)
result = blast.search("AGCTTTTCATTCTGACTGCA", top_n=5)

for hit in result.hits:
    print(f"{hit.rank}: {hit.hit.subject_id}, score={hit.raw_score}")
```

---

## 📂 项目结构

```
blastmini/
├── .github/workflows/     # CI/CD 工作流（pytest + flake8 + docs），已经上传github
├── benchmarks/            # 性能基准测试脚本 + 结果图表
├── docs/                  # Sphinx 文档源码
├── workflows/             # Snakemake 可复现工作流
├── src/blastmini/         # 核心源码
│   ├── models.py          # 数据模型（SequenceRecord, Hit, AlignmentConfig）
│   ├── io.py              # FASTA 解析与 I/O
│   ├── index.py           # k‑mer 索引构建
│   ├── seeding.py         # 种子搜索
│   ├── extension.py       # X‑dropoff 双向延伸
│   ├── scoring.py         # 打分、排序与格式化
│   ├── stats.py           # 统计显著性估计（E‑value / p‑value）
│   ├── api.py             # 高级 Python API
│   └── cli.py             # 命令行接口
├── tests/                 # pytest 单元测试
├── notebooks/             # 演示 PPT
├── data/                  # 测试数据（大肠杆菌基因组等）
├── Dockerfile             # 容器化配置
├── environment.yml        # Conda 完整开发环境
├── pyproject.toml         # 项目元数据与依赖声明
└── README.md              # 项目说明（本文件）
```

---

## 🧪 运行测试

```bash
pytest tests/ -v
```

---

## 📖 文档

在线文档：[https://linchenghang3-hub.github.io/blastmini/](https://linchenghang3-hub.github.io/blastmini/)

并在本地构建文档：

```bash
cd docs
make html
# 打开 build/html/index.html
```

---

## 🐳 Docker 用法

```bash
# 构建镜像
docker build -t blastmini:latest .

# 构建索引（挂载数据目录）
docker run --rm -v $(pwd)/data:/data blastmini:latest \
    blastmini build -i /data/NC_000913.3.fa -o /data/index.json -k 11

# 搜索
docker run --rm -v $(pwd)/data:/data blastmini:latest \
    blastmini search -q /data/queries.fa -d /data/index.json \
    --database-fasta /data/NC_000913.3.fa -o /data/results.tsv

# 运行 Snakemake 工作流
docker run --rm -v $(pwd):/workspace -w /workspace blastmini:latest \
    snakemake -s workflows/Snakefile --cores 2
```

---

## 📊 与 BLASTn 对比结果

运行 Snakemake 工作流后，对比图保存在 `results/comparison.png`，报告在 `results/comparison_report.txt`。

示例报告（部分）：

```
==================================================
COMPARISON REPORT: blastmini vs blastn
==================================================
Number of common queries: 50
Score correlation: 0.985
Identity correlation: 0.997
Same subject ID ratio: 100.00%
```

---

## 🤝 贡献

欢迎提交 Issue 和 Pull Request。请确保代码通过 `flake8` 检查且 `pytest` 全部通过。

---

## 📄 分工

在本项目中，吴宇轩负责基础功能和scr文件夹下的实现，PPT的制作，林成航负责所有代码的测试和benchmarks、docs、Snakemeke、github上传、data数据库构建、Docker容器化、pyproject配置、environment环境编写，汤羽翔负责后期整理、介绍和视频剪辑。

---

## 🙏 致谢

本项目为 BIO2502 课程项目。感谢老师一学期以来生物编程的教学。

本项目部分使用deepseek辅助编程，为我们提供了很大帮助。

---

## 🚀 提交并推送

```bash
cd /home/linchenghang/biancheng/blastmini

# 添加 README
git add README.md

# 如果有其他未提交的更改，一并添加
git add .

# 提交
git commit -m "Add comprehensive README.md"

# 推送到远程 main 分支
git push origin HEAD:main
```

