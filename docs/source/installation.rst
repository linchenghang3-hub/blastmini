Installation
============

blastmini 支持多种安装方式。推荐使用 **pip 开发模式安装** 以便于测试和修改源码。

系统要求
--------

- Python 3.9 或更高版本
- 建议使用 Linux 或 WSL 2 环境
- 可选：Docker（用于容器化运行）

安装方法
--------

1. 从源码安装（开发模式）
~~~~~~~~~~~~~~~~~~~~~~~~~~

将包链接到你的 Python 环境中，代码修改后立即生效。

.. code-block:: bash

    # 克隆或进入项目根目录
    cd /path/to/blastmini

    # 安装依赖并注册命令行工具
    pip install -e .

    # 验证安装
    blastmini --help

2. 安装可选依赖（用于工作流和基准测试）
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

如果你想运行 Snakemake 工作流或生成性能图表，建议安装以下依赖：

.. code-block:: bash

    pip install snakemake pandas matplotlib seaborn psutil memory-profiler

或者直接使用项目提供的 `requirements.txt`。

3. 使用 Docker（无需配置环境）
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

如果你不想在本地安装 Python 环境，可以使用已构建的 Docker 镜像：

.. code-block:: bash

    # 构建镜像
    docker build -t blastmini:latest .

    # 运行搜索（需挂载数据目录）
    docker run --rm -v $(pwd)/data:/data blastmini:latest blastmini search -q /data/queries.fa -d /data/index.json --database-fasta /data/NC_000913.3.fa

验证安装
--------

打开终端，输入以下命令，如果显示帮助信息则说明安装成功：

.. code-block:: bash

    blastmini --version
    # 输出示例: blastmini 0.1.0