"""Command-line interface for blastmini.

This module provides a unified CLI for all blastmini functionality,
including index building, sequence searching, and result formatting.

Usage:
    blastmini build -i database.fa -k 11 -o index.json
    blastmini search -q query.fa -d index.json --top 10 -o results.tsv
    blastmini stats -i index.json
    blastmini view -r results.tsv --format text
"""

import argparse
import os
import sys
import time

from .extension import SeedExtender
from .index import build_index_from_fasta, load_index_from_file
from .io import load_hits_from_tsv, parse_fasta
from .models import AlignmentConfig
from .scoring import (HitScorer, ScoredHit, format_hits_as_bed,
                      format_hits_as_json, format_hits_as_text,
                      format_hits_as_tsv, visualize_alignment)
from .seeding import SeedFinder
from .stats import SignificanceEstimator, format_significance_results

# ============================================================================
# Utility Functions
# ============================================================================


def print_error(message: str) -> None:
    """Print error message to stderr."""
    print(f"ERROR: {message}", file=sys.stderr)


def print_info(message: str) -> None:
    """Print info message to stderr."""
    print(f"INFO: {message}", file=sys.stderr)


def print_success(message: str) -> None:
    """Print success message to stderr."""
    print(f"✓ {message}", file=sys.stderr)


def format_duration(seconds: float) -> str:
    """Format duration in human-readable format."""
    if seconds < 60:
        return f"{seconds:.2f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m {seconds % 60:.0f}s"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}h"


def validate_fasta_file(filepath: str) -> bool:
    """Validate that a FASTA file exists and is readable."""
    if not os.path.exists(filepath):
        print_error(f"File not found: {filepath}")
        return False

    try:
        # Try to parse at least one record
        for _ in parse_fasta(filepath):
            return True
    except Exception as e:
        print_error(f"Failed to parse FASTA file {filepath}: {e}")
        return False

    print_error(f"No sequences found in {filepath}")
    return False


# ============================================================================
# Command Handlers
# ============================================================================

def cmd_build(args: argparse.Namespace) -> int:
    """Build k-mer index from FASTA file."""
    start_time = time.time()

    print_info(f"Building index from {args.input}")
    print_info(f"  k-mer size: {args.kmer}")
    print_info(f"  Canonical: {not args.no_canonical}")

    if args.min_occurrences > 1 or args.max_occurrences is not None:
        print_info(
            f"  Filter: min={args.min_occurrences}, max={args.max_occurrences or 'unlimited'}")

    try:
        # Build index
        idx = build_index_from_fasta(
            fasta_file=args.input,
            k=args.kmer,
            canonical=not args.no_canonical,
            min_occurrences=args.min_occurrences,
            max_occurrences=args.max_occurrences,
            output_file=args.output,
            progress=not args.quiet
        )

        # Print statistics
        stats = idx.stats
        print_success("Index built successfully")
        print(f"  Sequences: {stats.total_sequences:,}")
        print(f"  Total length: {stats.total_length:,} bp")
        print(f"  Unique k-mers: {stats.unique_kmers:,}")
        print(f"  Total k-mer positions: {stats.total_kmers:,}")
        print(f"  Memory estimate: {stats.memory_estimate_mb:.2f} MB")
        print(f"  Time: {format_duration(time.time() - start_time)}")

        return 0

    except Exception as e:
        print_error(f"Failed to build index: {e}")
        return 1


def cmd_search(args: argparse.Namespace) -> int:
    """Search query sequences against indexed database."""
    start_time = time.time()

    # Load index
    print_info(f"Loading index from {args.database}")
    try:
        idx = load_index_from_file(args.database)
    except Exception as e:
        print_error(f"Failed to load index: {e}")
        return 1

    print_info(
        f"  k-mer size: {idx.k}, sequences: {idx.stats.total_sequences:,}")

    # Load queries
    print_info(f"Loading queries from {args.query}")
    try:
        queries = list(parse_fasta(args.query))
        if args.max_queries and len(queries) > args.max_queries:
            queries = queries[:args.max_queries]
            print_info(f"  Using first {len(queries)} queries")
    except Exception as e:
        print_error(f"Failed to load queries: {e}")
        return 1

    if not queries:
        print_error("No queries found")
        return 1

    print_info(f"  Found {len(queries)} query sequences")

    # Load subject sequences (for extension)
    if args.database_fasta:
        print_info(f"Loading subject sequences from {args.database_fasta}")
        try:
            subject_sequences = {}
            for record in parse_fasta(args.database_fasta):
                subject_sequences[record.id] = record.sequence
            print_info(f"  Loaded {len(subject_sequences)} subject sequences")
        except Exception as e:
            print_error(f"Failed to load subject sequences: {e}")
            return 1
    else:
        print_error(
            "Database FASTA file required for extension (--database-fasta)")
        return 1

    # Create configuration
    config = AlignmentConfig(
        kmer_size=idx.k,
        match_score=args.match_score,
        mismatch_penalty=args.mismatch_penalty,
        x_dropoff=args.dropoff,
        top_n=args.top
    )

    # Create components
    finder = SeedFinder(idx, config=config)
    extender = SeedExtender(config=config, min_extension_score=args.min_score)
    scorer = HitScorer(
        config=config,
        min_score=args.min_score,
        min_identity=args.min_identity,
        top_n=args.top
    )

    # Process each query
    all_results = []
    total_seeds = 0
    total_extensions = 0

    for query_idx, query in enumerate(queries, 1):
        if not args.quiet:
            print_info(
                f"Processing query {query_idx}/{len(queries)}: {query.id}")

        # Find seeds
        seeds = finder.find_best_seeds(
            query=query,
            top_n=args.max_seeds,
            subject_ids=list(subject_sequences.keys())
        )

        total_seeds += len(seeds)

        if not seeds:
            if not args.quiet:
                print_info(f"  No seeds found for {query.id}")
            continue

        if not args.quiet:
            print_info(f"  Found {len(seeds)} seeds")

        # Extend seeds
        results = extender.extend_seeds(
            query=query,
            seeds=seeds,
            subject_sequences=subject_sequences,
            max_results=args.max_results,
            progress=not args.quiet
        )

        total_extensions += len(results)

        # Score and rank
        scored_hits = scorer.score_hits(
            hits=[r.to_hit() for r in results],
            query_length=len(query),
            db_total_length=idx.stats.total_length,
            progress=False
        )

        all_results.extend(scored_hits)

    # Sort all results
    all_results.sort(key=lambda s: (
        s.raw_score, s.identity_percent), reverse=True)

    # Assign global ranks
    for i, scored in enumerate(all_results):
        scored.rank = i + 1

    # Keep top N
    if args.top and len(all_results) > args.top:
        all_results = all_results[:args.top]

    # Output results
    if args.output:
        if args.output_format == 'tsv':
            output = format_hits_as_tsv(all_results)
            with open(args.output, 'w') as f:
                f.write(output)
        elif args.output_format == 'json':
            output = format_hits_as_json(all_results)
            with open(args.output, 'w') as f:
                f.write(output)
        elif args.output_format == 'bed':
            output = format_hits_as_bed(all_results)
            with open(args.output, 'w') as f:
                f.write(output)
        else:
            output = format_hits_as_text(all_results, top_n=args.top)
            with open(args.output, 'w') as f:
                f.write(output)
        print_success(f"Results saved to {args.output}")
    else:
        # Print to stdout
        if args.output_format == 'tsv':
            print(format_hits_as_tsv(all_results))
        elif args.output_format == 'json':
            print(format_hits_as_json(all_results))
        elif args.output_format == 'bed':
            print(format_hits_as_bed(all_results))
        else:
            print(format_hits_as_text(all_results, top_n=args.top))

    # Summary
    print_info("Summary:")
    print(f"  Queries processed: {len(queries)}")
    print(f"  Total seeds found: {total_seeds}")
    print(f"  Total extensions: {total_extensions}")
    print(f"  Results saved: {len(all_results)}")
    print(f"  Time: {format_duration(time.time() - start_time)}")

    return 0


def cmd_stats(args: argparse.Namespace) -> int:
    """Display index statistics."""
    print_info(f"Loading index from {args.index}")

    try:
        idx = load_index_from_file(args.index)
    except Exception as e:
        print_error(f"Failed to load index: {e}")
        return 1

    stats = idx.stats

    print("Index Statistics")
    print("=" * 40)
    print(f"  k-mer size: {stats.kmer_size}")
    print(f"  Sequences: {stats.total_sequences:,}")
    print(f"  Total length: {stats.total_length:,} bp")
    print(f"  Unique k-mers: {stats.unique_kmers:,}")
    print(f"  Total positions: {stats.total_kmers:,}")
    print(f"  Memory estimate: {stats.memory_estimate_mb:.2f} MB")

    if args.verbose:
        print(
            f"  Average positions per k-mer: {stats.total_kmers / stats.unique_kmers:.2f}")

        # Get k-mer frequency distribution
        kmers_with_counts = idx.get_kmers_with_counts()

        if kmers_with_counts:
            print("  Most frequent k-mers:")
            for kmer, count in kmers_with_counts[:10]:
                print(f"    {kmer}: {count} occurrences")

    return 0


def cmd_view(args: argparse.Namespace) -> int:
    """View and visualize search results."""
    print_info(f"Loading results from {args.results}")

    try:
        hits = load_hits_from_tsv(args.results)
    except Exception as e:
        print_error(f"Failed to load results: {e}")
        return 1

    if not hits:
        print_error("No hits found in results file")
        return 1

    print_info(f"Loaded {len(hits)} hits")

    # 将 Hit 对象转换为 ScoredHit 以便格式化
    scored_hits = []
    for i, hit in enumerate(hits, 1):
        scored = ScoredHit(
            hit=hit,
            raw_score=hit.score,
            identity_percent=hit.identity_percent,
            coverage_percent=0.0,          # 未提供查询长度，设为0
            bit_score=hit.bit_score or 0.0,
            evalue=hit.evalue or 1.0,
            rank=i
        )
        scored_hits.append(scored)

    if args.format == 'text':
        output = format_hits_as_text(scored_hits, top_n=args.top)
    elif args.format == 'tsv':
        output = format_hits_as_tsv(scored_hits)
    elif args.format == 'json':
        output = format_hits_as_json(scored_hits)
    elif args.format == 'bed':
        output = format_hits_as_bed(scored_hits)
    else:
        output = format_hits_as_text(scored_hits, top_n=args.top)

    print(output)

    # 如果需要显示比对可视化
    if args.show_alignment and scored_hits:
        print("\n" + "=" * 80)
        print("Alignment Visualization (Top Hit)")
        print("=" * 80)

        top_hit = scored_hits[0].hit  # 从 ScoredHit 中提取原始 Hit
        if top_hit.query_alignment and top_hit.subject_alignment:
            print(visualize_alignment(
                top_hit.query_alignment,
                top_hit.subject_alignment,
                line_width=args.line_width
            ))
        else:
            print("No alignment data available")

    return 0


def cmd_significance(args: argparse.Namespace) -> int:
    """Estimate statistical significance of hits."""
    start_time = time.time()

    print_info(f"Loading results from {args.results}")

    try:
        hits = load_hits_from_tsv(args.results)
    except Exception as e:
        print_error(f"Failed to load results: {e}")
        return 1

    if not hits:
        print_error("No hits found in results file")
        return 1

    print_info(f"Loaded {len(hits)} hits")

    # Load database for background estimation
    if args.database:
        print_info(f"Loading database from {args.database}")
        try:
            subject_sequences = {}
            for record in parse_fasta(args.database):
                subject_sequences[record.id] = record.sequence
            print_info(f"  Loaded {len(subject_sequences)} sequences")
        except Exception as e:
            print_error(f"Failed to load database: {e}")
            return 1
    else:
        print_error(
            "Database FASTA file required for significance estimation (--database)")
        return 1

    # Load query
    if args.query:
        print_info(f"Loading query from {args.query}")
        try:
            queries = list(parse_fasta(args.query))
            if not queries:
                print_error("No queries found in file")
                return 1
            query = queries[0]
        except Exception as e:
            print_error(f"Failed to load query: {e}")
            return 1
    else:
        print_error("Query file required for significance estimation (--query)")
        return 1

    # Calculate database size
    db_size = sum(len(seq) for seq in subject_sequences.values())

    # Create estimator
    estimator = SignificanceEstimator(
        n_permutations=args.permutations,
        random_seed=args.seed,
        min_score_threshold=args.min_score
    )

    # Estimate background distribution
    if not args.quiet:
        print_info(
            f"Estimating background distribution with {args.permutations} permutations...")

    # Estimate significance
    if not args.quiet:
        print_info(f"Estimating significance for {len(hits)} hits...")

    results = estimator.estimate_significance_batch(
        hits=hits,
        query_length=len(query.sequence),
        db_size=db_size,
        use_extreme_distribution=args.use_extreme,
        progress=not args.quiet
    )

    # Apply multiple testing correction
    if args.correction:
        results = estimator.adjust_for_multiple_testing(
            results, method=args.correction)

    # Output results
    if args.output:
        output = format_significance_results(results, format=args.format)
        with open(args.output, 'w') as f:
            f.write(output)
        print_success(f"Significance results saved to {args.output}")
    else:
        output = format_significance_results(results, format=args.format)
        print(output)

    # Summary
    significant = sum(1 for r in results if r.is_significant)
    print_info("Summary:")
    print(f"  Hits analyzed: {len(results)}")
    print(f"  Significant hits: {significant}")
    print(f"  Time: {format_duration(time.time() - start_time)}")

    return 0


def cmd_info(args: argparse.Namespace) -> int:
    """Display package information."""
    from . import __version__

    print(f"blastmini v{__version__}")
    print("A lightweight BLAST implementation for educational purposes")
    print("")
    print("Features:")
    print("  - k-mer indexing with canonical support")
    print("  - Seed-and-extend search")
    print("  - X-dropoff extension")
    print("  - Statistical significance estimation")
    print("  - Multiple output formats (TSV, JSON, BED)")
    print("")
    print("Usage:")
    print("  blastmini build -i database.fa -o index.json")
    print("  blastmini search -q query.fa -d index.json --database-fasta database.fa")
    print("  blastmini stats -i index.json")
    print("  blastmini view -r results.tsv")
    print("  blastmini significance -r results.tsv -d database.fa -q query.fa")

    return 0


# ============================================================================
# Command Line Parser
# ============================================================================

def create_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        prog="blastmini",
        description="Lightweight BLAST implementation for educational purposes"
    )

    parser.add_argument(
        "--version",
        action="store_true",
        help="Show version and exit"
    )

    # Subparsers for commands
    subparsers = parser.add_subparsers(
        dest="command",
        help="Command to execute",
        required=False
    )

    # ===== build command =====
    build_parser = subparsers.add_parser(
        "build",
        help="Build k-mer index from FASTA file"
    )
    build_parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file"
    )
    build_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output index file (JSON)"
    )
    build_parser.add_argument(
        "-k", "--kmer",
        type=int,
        default=11,
        help="k-mer size (default: 11)"
    )
    build_parser.add_argument(
        "--no-canonical",
        action="store_true",
        help="Disable canonical k-mers (default: enabled)"
    )
    build_parser.add_argument(
        "--min-occurrences",
        type=int,
        default=1,
        help="Minimum occurrences for a k-mer (default: 1)"
    )
    build_parser.add_argument(
        "--max-occurrences",
        type=int,
        help="Maximum occurrences for a k-mer (prune repetitive)"
    )
    build_parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress progress messages"
    )
    build_parser.set_defaults(func=cmd_build)

    # ===== search command =====
    search_parser = subparsers.add_parser(
        "search",
        help="Search query sequences against indexed database"
    )
    search_parser.add_argument(
        "-q", "--query",
        required=True,
        help="Query FASTA file"
    )
    search_parser.add_argument(
        "-d", "--database",
        required=True,
        help="Index file (JSON)"
    )
    search_parser.add_argument(
        "--database-fasta",
        required=True,
        help="Database FASTA file (for extension)"
    )
    search_parser.add_argument(
        "-o", "--output",
        help="Output file (if not specified, prints to stdout)"
    )
    search_parser.add_argument(
        "-k", "--kmer",
        type=int,
        help="k-mer size (overrides index value)"
    )
    search_parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of top hits to report (default: 10)"
    )
    search_parser.add_argument(
        "--max-seeds",
        type=int,
        default=50,
        help="Maximum seeds per query (default: 50)"
    )
    search_parser.add_argument(
        "--max-results",
        type=int,
        default=100,
        help="Maximum results per query (default: 100)"
    )
    search_parser.add_argument(
        "--max-queries",
        type=int,
        help="Maximum queries to process"
    )
    search_parser.add_argument(
        "--match-score",
        type=int,
        default=5,
        help="Match score (default: 5)"
    )
    search_parser.add_argument(
        "--mismatch-penalty",
        type=int,
        default=-4,
        help="Mismatch penalty (default: -4)"
    )
    search_parser.add_argument(
        "--dropoff",
        type=int,
        default=5,
        help="X-dropoff threshold (default: 5)"
    )
    search_parser.add_argument(
        "--min-score",
        type=int,
        default=10,
        help="Minimum score for a hit (default: 10)"
    )
    search_parser.add_argument(
        "--min-identity",
        type=float,
        default=70.0,
        help="Minimum identity percentage (default: 70.0)"
    )
    search_parser.add_argument(
        "--format",
        dest="output_format",
        choices=["text", "tsv", "json", "bed"],
        default="text",
        help="Output format (default: text)"
    )
    search_parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress messages"
    )
    search_parser.set_defaults(func=cmd_search)

    # ===== stats command =====
    stats_parser = subparsers.add_parser(
        "stats",
        help="Display index statistics"
    )
    stats_parser.add_argument(
        "-i", "--index",
        required=True,
        help="Index file (JSON)"
    )
    stats_parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Show detailed statistics"
    )
    stats_parser.set_defaults(func=cmd_stats)

    # ===== view command =====
    view_parser = subparsers.add_parser(
        "view",
        help="View and visualize search results"
    )
    view_parser.add_argument(
        "-r", "--results",
        required=True,
        help="Results file (TSV)"
    )
    view_parser.add_argument(
        "--format",
        choices=["text", "tsv", "json", "bed"],
        default="text",
        help="Output format (default: text)"
    )
    view_parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Number of top hits to show (default: 20)"
    )
    view_parser.add_argument(
        "--show-alignment",
        action="store_true",
        help="Show alignment visualization"
    )
    view_parser.add_argument(
        "--line-width",
        type=int,
        default=60,
        help="Alignment line width (default: 60)"
    )
    view_parser.set_defaults(func=cmd_view)

    # ===== significance command =====
    significance_parser = subparsers.add_parser(
        "significance",
        help="Estimate statistical significance of hits"
    )
    significance_parser.add_argument(
        "-r", "--results",
        required=True,
        help="Results file (TSV)"
    )
    significance_parser.add_argument(
        "-d", "--database",
        required=True,
        help="Database FASTA file"
    )
    significance_parser.add_argument(
        "-q", "--query",
        required=True,
        help="Query FASTA file"
    )
    significance_parser.add_argument(
        "-o", "--output",
        help="Output file"
    )
    significance_parser.add_argument(
        "-n", "--permutations",
        type=int,
        default=100,
        help="Number of permutations (default: 100)"
    )
    significance_parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for reproducibility"
    )
    significance_parser.add_argument(
        "--format",
        choices=["text", "tsv", "json"],
        default="text",
        help="Output format (default: text)"
    )
    significance_parser.add_argument(
        "--correction",
        choices=["bonferroni", "fdr"],
        help="Multiple testing correction method"
    )
    significance_parser.add_argument(
        "--use-extreme",
        action="store_true",
        help="Use extreme value distribution (faster)"
    )
    significance_parser.add_argument(
        "--min-score",
        type=int,
        default=10,
        help="Minimum score threshold (default: 10)"
    )
    significance_parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress progress messages"
    )
    significance_parser.set_defaults(func=cmd_significance)

    # ===== info command =====
    info_parser = subparsers.add_parser(
        "info",
        help="Display package information"
    )
    info_parser.set_defaults(func=cmd_info)

    return parser


# ============================================================================
# Main Entry Point
# ============================================================================

def main() -> int:
    parser = create_parser()
    args = parser.parse_args()

    # Handle version flag
    if args.version:
        from . import __version__
        print(f"blastmini {__version__}")
        return 0

    # If no command is given, show help
    if args.command is None:
        parser.print_help()
        return 1

    # Execute command
    if hasattr(args, 'func'):
        try:
            return args.func(args)
        except KeyboardInterrupt:
            print("\nInterrupted by user", file=sys.stderr)
            return 130
        except Exception as e:
            print_error(f"Unexpected error: {e}")
            if args.func.__name__ in ['cmd_search', 'cmd_significance']:
                print_info(
                    "Try running with --quiet to suppress progress messages")
            return 1
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
