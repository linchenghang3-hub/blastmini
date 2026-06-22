"""
blastmini - A lightweight BLAST implementation for educational purposes.

blastmini provides a minimal but functional implementation of the BLAST
algorithm, including k-mer indexing, seed-and-extend search, X-dropoff
extension, scoring, and statistical significance estimation.
"""

__version__ = "0.1.0"

# Core data models
from .models import SequenceRecord, Hit, AlignmentConfig

# I/O utilities
from .io import (
    parse_fasta,
    parse_fasta_multi,
    save_hits_to_tsv,
    load_hits_from_tsv,
    save_index,
    load_index,
    get_sequence_stats,
)

# Indexing
from .index import KmerIndex, build_index_from_fasta, load_index_from_file

# Seeding
from .seeding import Seed, SeedCluster, SeedFinder

# Extension
from .extension import ExtensionResult, ExtensionStats, SeedExtender

# Scoring
from .scoring import (
    HitScorer,
    ScoredHit,
    format_hits_as_text,
    format_hits_as_tsv,
    format_hits_as_json,
    format_hits_as_bed,
    visualize_alignment,
)

# Statistics
from .stats import (
    ScoreDistribution,
    SignificanceResult,
    SignificanceEstimator,
    format_significance_results,
)

# High-level API
from .api import (
    BlastMini,
    SearchResult,
    BatchSearchResult,
    quick_search,
    quick_batch_search,
)

# CLI entry point (optional, can be used programmatically)
from .cli import main as cli_main

# Public API
__all__ = [
    # version
    "__version__",
    # models
    "SequenceRecord",
    "Hit",
    "AlignmentConfig",
    # io
    "parse_fasta",
    "parse_fasta_multi",
    "save_hits_to_tsv",
    "load_hits_from_tsv",
    "save_index",
    "load_index",
    "get_sequence_stats",
    # index
    "KmerIndex",
    "build_index_from_fasta",
    "load_index_from_file",
    # seeding
    "Seed",
    "SeedCluster",
    "SeedFinder",
    # extension
    "ExtensionResult",
    "ExtensionStats",
    "SeedExtender",
    # scoring
    "HitScorer",
    "ScoredHit",
    "format_hits_as_text",
    "format_hits_as_tsv",
    "format_hits_as_json",
    "format_hits_as_bed",
    "visualize_alignment",
    # stats
    "ScoreDistribution",
    "SignificanceResult",
    "SignificanceEstimator",
    "format_significance_results",
    # api
    "BlastMini",
    "SearchResult",
    "BatchSearchResult",
    "quick_search",
    "quick_batch_search",
    # cli
    "cli_main",
]