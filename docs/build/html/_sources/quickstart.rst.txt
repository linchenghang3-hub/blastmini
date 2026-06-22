Quickstart
==========

本文档将带你快速体验 blastmini 的完整工作流：构建索引、搜索序列，并查看结果。

准备测试数据
------------

我们使用大肠杆菌（*E. coli* K-12）的基因组作为数据库，并随机截取一条短序列作为查询。

.. code-block:: bash

    # 确保项目 data 目录下存在基因组文件
    # 如果还没有，可以从 NCBI 下载或使用项目中已有的 data/NC_000913.3.fa

    # 查看基因组文件
    head -n 2 data/NC_000913.3.fa

    # 手动创建一条简单的查询序列（也可从基因组随机截取）
    echo -e ">query1\nAGCTTTTCATTCTGACTGCA" > my_query.fa

构建索引
--------

使用 `blastmini build` 命令构建 k-mer 索引。k 值越大，搜索越快，但内存占用也越高。

.. code-block:: bash

    blastmini build -i data/NC_000913.3.fa -o data/index.json -k 11

执行搜索
--------

使用 `blastmini search` 对查询序列进行搜索，输出格式可选 `text`、`tsv`、`json` 或 `bed`。

.. code-block:: bash

    blastmini search -q my_query.fa -d data/index.json --database-fasta data/NC_000913.3.fa -o results.tsv --format tsv --top 5

查看结果
--------

使用 `blastmini view` 查看格式化后的比对结果：

.. code-block:: bash

    blastmini view -r results.tsv --show-alignment

你应该会看到类似下面的输出，包含得分、比对长度和一致性百分比：

::

    ================================================================================
    BLAST Search Results (5 hits found)
    ================================================================================

    Rank: 1
    Subject: NC_000913.3
    Score: 1000 (bitscore: 144.3)
    E-value: 1.00e-250
    Identity: 100.0%
    ...

运行完整工作流（可选）
----------------------

如果你想一次性完成数据生成、索引构建、搜索和对比，可以直接运行项目自带的 Snakemake 工作流：

.. code-block:: bash

    snakemake -s workflows/Snakefile --cores 2

该流程会：
1. 从基因组随机截取 50 条查询序列
2. 分别用 blastmini 和 NCBI BLASTn 进行搜索
3. 生成对比图表和统计报告

Python API 示例
---------------

如果你更喜欢在 Python 脚本中使用，可以这样调用：

.. code-block:: python

    from blastmini.api import BlastMini

    # 加载数据库
    blast = BlastMini.from_fasta("data/NC_000913.3.fa", k=11)

    # 搜索
    result = blast.search("AGCTTTTCATTCTGACTGCA", top_n=5)

    # 打印结果
    for hit in result.hits:
        print(f"Hit {hit.rank}: {hit.hit.subject_id}, score={hit.raw_score}")