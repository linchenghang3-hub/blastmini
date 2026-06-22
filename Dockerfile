# 使用 Python 3.10 官方轻量级镜像
FROM python:3.10-slim as builder

# 设置工作目录
WORKDIR /app

# 更换 apt 源为清华镜像
RUN sed -i 's|deb.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list || true \
    && sed -i 's|deb.debian.org|mirrors.tuna.tsinghua.edu.cn|g' /etc/apt/sources.list.d/debian.sources || true

# 安装系统依赖（编译工具和库）
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# 复制 pyproject.toml 并安装 Python 依赖（利用 Docker 缓存）
COPY pyproject.toml ./
COPY src/ ./src/

# 设置 pip 源为清华镜像（加速所有 pip 操作）
RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple

# 一次性安装 blastmini 包及所有额外工具
RUN pip install --no-cache-dir -e . snakemake pandas matplotlib seaborn psutil memory-profiler

# 复制其余代码（如果之前未复制）
# 注意：由于我们已复制 src/，但工作流和测试等可能还需要
COPY workflows/ ./workflows/
COPY benchmarks/ ./benchmarks/
COPY tests/ ./tests/
COPY notebooks/ ./notebooks/

# 设置环境变量（可选）
ENV PYTHONUNBUFFERED=1

# 默认命令：显示帮助信息
CMD ["blastmini", "--help"]