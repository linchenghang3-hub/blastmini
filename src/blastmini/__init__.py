"""
blastmini - A lightweight BLAST implementation for educational purposes.

blastmini provides a minimal but functional implementation of the BLAST
algorithm, including k-mer indexing, seed-and-extend search, X-dropoff
extension, scoring, and statistical significance estimation.
"""

__version__ = "0.1.0"

# High-level API
from .api import (BatchSearchResult, BlastMini, SearchResult,
                  quick_batch_search, quick_search)
# CLI entry point (optional, can be used programmatically)
from .cli import main as cli_main
# Extension
from .extension import ExtensionResult, ExtensionStats, SeedExtender
# Indexing
from .index import KmerIndex, build_index_from_fasta, load_index_from_file
# I/O utilities
from .io import (get_sequence_stats, load_hits_from_tsv, load_index,
                 parse_fasta, parse_fasta_multi, save_hits_to_tsv, save_index)
# Core data models
from .models import AlignmentConfig, Hit, SequenceRecord
# Scoring
from .scoring import (HitScorer, ScoredHit, format_hits_as_bed,
                      format_hits_as_json, format_hits_as_text,
                      format_hits_as_tsv, visualize_alignment)
# Seeding
from .seeding import Seed, SeedCluster, SeedFinder
# Statistics
from .stats import (ScoreDistribution, SignificanceEstimator,
                    SignificanceResult, format_significance_results)

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
