"""High-level API for blastmini.

This module provides a clean, user-friendly interface for all blastmini
functionality. It is designed for use in Python scripts, Jupyter notebooks,
and interactive sessions.

Key features:
    - Simple one-line search interface
    - Context managers for resource management
    - Batch processing support
    - Progress callbacks
    - Result caching
    - Easy access to all blastmini features

Examples:
    >>> from blastmini.api import BlastMini
    >>>
    >>> # Initialize with database
    >>> blast = BlastMini.from_fasta("database.fa", k=11)
    >>>
    >>> # Search a query
    >>> results = blast.search("ATCGATCGATCG", top_n=10)
    >>>
    >>> # Get top hit
    >>> top = results[0]
    >>> print(f"Best hit: {top.subject_id}, score={top.raw_score}")
"""

import os
import sys
import tempfile
from typing import Optional, List, Dict, Any, Union, Iterator, Callable, Tuple
from pathlib import Path
from dataclasses import dataclass, field
import json

from .models import SequenceRecord, Hit, AlignmentConfig
from .io import parse_fasta, save_hits_to_tsv, load_hits_from_tsv
from .index import KmerIndex, build_index_from_fasta, load_index_from_file
from .seeding import SeedFinder, Seed
from .extension import SeedExtender, ExtensionResult
from .scoring import HitScorer, ScoredHit, format_hits_as_text, format_hits_as_tsv
from .stats import SignificanceEstimator, SignificanceResult, format_significance_results

# ============================================================================
# Progress Callback Types
# ============================================================================

ProgressCallback = Callable[[str, int, int], None]


def default_progress_callback(stage: str, current: int, total: int) -> None:
    """Default progress callback that prints to stderr."""
    if total > 0:
        percent = (current / total) * 100
        sys.stderr.write(f"\r{stage}: {percent:.1f}% ({current}/{total})")
        if current >= total:
            sys.stderr.write("\n")


# ============================================================================
# Search Result Container
# ============================================================================

@dataclass
class SearchResult:
    """Container for search results.

    This class provides easy access to search results and associated
    metadata.

    Attributes:
        query_id: Query sequence identifier.
        query_sequence: Query sequence string.
        hits: List of scored hits.
        num_seeds_found: Number of seeds found.
        num_extensions: Number of extensions performed.
        search_time: Search time in seconds.
        parameters: Search parameters used.
    """
    query_id: str
    query_sequence: str
    hits: List[ScoredHit]
    num_seeds_found: int = 0
    num_extensions: int = 0
    search_time: float = 0.0
    parameters: Dict[str, Any] = field(default_factory=dict)

    @property
    def top_hit(self) -> Optional[ScoredHit]:
        """Get the top-scoring hit."""
        return self.hits[0] if self.hits else None

    @property
    def num_hits(self) -> int:
        """Get the number of hits."""
        return len(self.hits)

    @property
    def has_significant_hits(self) -> bool:
        """Check if there are any significant hits."""
        return any(h.is_significant for h in self.hits if hasattr(h, 'is_significant'))

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'query_id': self.query_id,
            'query_sequence': self.query_sequence[:100] + "..." if len(
                self.query_sequence) > 100 else self.query_sequence,
            'num_hits': self.num_hits,
            'num_seeds_found': self.num_seeds_found,
            'num_extensions': self.num_extensions,
            'search_time': self.search_time,
            'parameters': self.parameters,
            'hits': [h.to_dict() for h in self.hits]
        }

    def to_text(self) -> str:
        """Format as human-readable text."""
        return format_hits_as_text(self.hits)

    def to_tsv(self) -> str:
        """Format as TSV."""
        return format_hits_as_tsv(self.hits)

    def to_json(self) -> str:
        """Format as JSON."""
        return json.dumps(self.to_dict(), indent=2)

    def save(self, filepath: Union[str, Path], format: str = 'tsv') -> None:
        """Save results to file.

        Args:
            filepath: Output file path.
            format: Output format ('tsv', 'text', 'json').
        """
        if format == 'tsv':
            content = self.to_tsv()
        elif format == 'json':
            content = self.to_json()
        else:
            content = self.to_text()

        with open(filepath, 'w') as f:
            f.write(content)

    def __repr__(self) -> str:
        return f"SearchResult(query='{self.query_id}', hits={self.num_hits})"


# ============================================================================
# Batch Search Result
# ============================================================================

@dataclass
class BatchSearchResult:
    """Container for batch search results.

    Attributes:
        results: List of SearchResult objects.
        total_hits: Total number of hits across all queries.
        total_time: Total search time in seconds.
    """
    results: List[SearchResult] = field(default_factory=list)
    total_time: float = 0.0

    @property
    def total_hits(self) -> int:
        """Get total number of hits."""
        return sum(r.num_hits for r in self.results)

    @property
    def total_queries(self) -> int:
        """Get total number of queries."""
        return len(self.results)

    def get_best_hits(self, top_n: int = 10) -> List[ScoredHit]:
        """Get the best hits across all queries."""
        all_hits = []
        for result in self.results:
            all_hits.extend(result.hits)
        all_hits.sort(key=lambda h: h.raw_score, reverse=True)
        return all_hits[:top_n]

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'total_queries': self.total_queries,
            'total_hits': self.total_hits,
            'total_time': self.total_time,
            'results': [r.to_dict() for r in self.results]
        }

    def to_json(self) -> str:
        """Format as JSON."""
        return json.dumps(self.to_dict(), indent=2)

    def save(self, filepath: Union[str, Path], format: str = 'tsv') -> None:
        """Save all results to a single file.

        Args:
            filepath: Output file path.
            format: Output format ('tsv', 'json').
        """
        if format == 'tsv':
            # Combine all TSV lines
            lines = []
            for result in self.results:
                lines.append(result.to_tsv())
            content = "\n".join(lines)
        else:
            content = self.to_json()

        with open(filepath, 'w') as f:
            f.write(content)

    def __repr__(self) -> str:
        return f"BatchSearchResult(queries={self.total_queries}, hits={self.total_hits})"


# ============================================================================
# Main API Class
# ============================================================================

class BlastMini:
    """Main API class for blastmini.

    This class provides a high-level interface for all blastmini
    functionality. It manages the database index, configuration,
    and provides convenient methods for searching.

    Attributes:
        config: Alignment configuration.
        index: K-mer index for the database.
        subject_sequences: Dictionary of subject sequences.
        verbose: Whether to print progress messages.

    Examples:
        >>> # Load from FASTA
        >>> blast = BlastMini.from_fasta("database.fa", k=11)
        >>>
        >>> # Search a single query
        >>> results = blast.search("ATCGATCGATCG")
        >>> print(results.top_hit)
        >>>
        >>> # Search multiple queries from a file
        >>> results = blast.search_fasta("queries.fa")
        >>>
        >>> # Configure search parameters
        >>> blast = BlastMini.from_fasta("database.fa", config=AlignmentConfig(
        ...     match_score=2, mismatch_penalty=-3, x_dropoff=10
        ... ))
    """

    def __init__(
            self,
            index: KmerIndex,
            subject_sequences: Optional[Dict[str, str]] = None,
            config: Optional[AlignmentConfig] = None,
            verbose: bool = True
    ):
        """Initialize BlastMini with an index.

        Args:
            index: K-mer index for the database.
            subject_sequences: Dictionary mapping subject ID to sequence.
            config: Alignment configuration.
            verbose: Whether to print progress messages.
        """
        self.index = index
        self.subject_sequences = subject_sequences or {}
        self.config = config or AlignmentConfig(kmer_size=index.k)
        self.verbose = verbose

        # Initialize components
        self._finder: Optional[SeedFinder] = None
        self._extender: Optional[SeedExtender] = None
        self._scorer: Optional[HitScorer] = None
        self._estimator: Optional[SignificanceEstimator] = None

        self._init_components()

    @classmethod
    def from_fasta(
            cls,
            fasta_file: Union[str, Path],
            k: int = 11,
            config: Optional[AlignmentConfig] = None,
            canonical: bool = True,
            min_occurrences: int = 1,
            max_occurrences: Optional[int] = None,
            max_records: Optional[int] = None,
            verbose: bool = True
    ) -> 'BlastMini':
        """Create BlastMini from a FASTA file.

        This is the recommended way to initialize BlastMini.

        Args:
            fasta_file: Path to FASTA file.
            k: k-mer size.
            config: Alignment configuration.
            canonical: Use canonical k-mers.
            min_occurrences: Minimum occurrences for k-mers.
            max_occurrences: Maximum occurrences for k-mers.
            max_records: Maximum records to index.
            verbose: Whether to print progress messages.

        Returns:
            BlastMini instance.

        Examples:
            >>> blast = BlastMini.from_fasta("database.fa", k=11)
        """
        # Build index
        if verbose:
            print(f"Building index from {fasta_file}...", file=sys.stderr)

        records = []
        for record in parse_fasta(fasta_file):
            records.append(record)
            if max_records and len(records) >= max_records:
                break

        if not records:
            raise ValueError(f"No sequences found in {fasta_file}")

        index = KmerIndex.build(
            records=records,
            k=k,
            canonical=canonical,
            min_kmer_occurrences=min_occurrences,
            max_kmer_occurrences=max_occurrences,
            progress=verbose
        )

        # Build subject sequences dict
        subject_sequences = {rec.id: rec.sequence for rec in records}

        # Create config
        if config is None:
            config = AlignmentConfig(kmer_size=k)

        return cls(
            index=index,
            subject_sequences=subject_sequences,
            config=config,
            verbose=verbose
        )

    @classmethod
    def from_index(
            cls,
            index_file: Union[str, Path],
            fasta_file: Optional[Union[str, Path]] = None,
            config: Optional[AlignmentConfig] = None,
            verbose: bool = True
    ) -> 'BlastMini':
        """Create BlastMini from a saved index file.

        Args:
            index_file: Path to index JSON file.
            fasta_file: Path to FASTA file (for subject sequences).
            config: Alignment configuration.
            verbose: Whether to print progress messages.

        Returns:
            BlastMini instance.
        """
        if verbose:
            print(f"Loading index from {index_file}...", file=sys.stderr)

        index = load_index_from_file(index_file)

        # Load subject sequences if provided
        subject_sequences = {}
        if fasta_file:
            for record in parse_fasta(fasta_file):
                subject_sequences[record.id] = record.sequence

        if config is None:
            config = AlignmentConfig(kmer_size=index.k)

        return cls(
            index=index,
            subject_sequences=subject_sequences,
            config=config,
            verbose=verbose
        )

    def _init_components(self) -> None:
        """Initialize internal components."""
        self._finder = SeedFinder(
            index=self.index,
            config=self.config
        )
        self._extender = SeedExtender(
            config=self.config,
            track_stats=True
        )
        self._scorer = HitScorer(
            config=self.config,
            min_score=0,
            min_identity=0,
            top_n=self.config.top_n
        )

    def _get_progress_callback(self) -> Optional[ProgressCallback]:
        """Get progress callback based on verbose setting."""
        if self.verbose:
            return default_progress_callback
        return None

    def search(
            self,
            query: Union[str, SequenceRecord],
            top_n: Optional[int] = None,
            min_score: Optional[int] = None,
            min_identity: Optional[float] = None,
            max_seeds: Optional[int] = None,
            max_extensions: Optional[int] = None,
            estimate_significance: bool = False,
            significance_params: Optional[Dict[str, Any]] = None
    ) -> SearchResult:
        """Search a single query sequence.

        Args:
            query: Query sequence or SequenceRecord.
            top_n: Number of top hits to return.
            min_score: Minimum score threshold.
            min_identity: Minimum identity percentage.
            max_seeds: Maximum seeds to consider.
            max_extensions: Maximum extensions to perform.
            estimate_significance: Whether to estimate statistical significance.
            significance_params: Parameters for significance estimation.

        Returns:
            SearchResult object.

        Examples:
            >>> result = blast.search("ATCGATCGATCG")
            >>> print(result.top_hit)
            >>>
            >>> # With custom parameters
            >>> result = blast.search("ATCGATCGATCG", top_n=5, min_score=50)
        """
        import time
        start_time = time.time()

        # Parse query
        if isinstance(query, SequenceRecord):
            query_id = query.id
            query_seq = query.sequence
        else:
            query_id = "query"
            query_seq = str(query)

        # Set parameters
        top_n = top_n or self.config.top_n
        min_score = min_score or 0
        min_identity = min_identity or 0.0
        max_seeds = max_seeds or 100
        max_extensions = max_extensions or 100

        # Update components for this search
        self._scorer.top_n = top_n
        self._scorer.min_score = min_score
        self._scorer.min_identity = min_identity

        # Find seeds
        if self.verbose:
            print(f"Searching query: {query_id}", file=sys.stderr)

        seeds = self._finder.find_best_seeds(
            query=query,
            top_n=max_seeds,
            subject_ids=list(self.subject_sequences.keys())
        )

        if not seeds:
            return SearchResult(
                query_id=query_id,
                query_sequence=query_seq,
                hits=[],
                num_seeds_found=0,
                search_time=time.time() - start_time,
                parameters={'top_n': top_n, 'min_score': min_score}
            )

        # Extend seeds
        results = self._extender.extend_seeds(
            query=query,
            seeds=seeds[:max_extensions],
            subject_sequences=self.subject_sequences,
            max_results=top_n * 2,
            progress=self.verbose
        )

        # Score hits
        hits = [r.to_hit() for r in results]
        scored_hits = self._scorer.score_hits(
            hits=hits,
            query_length=len(query_seq),
            db_total_length=self.index.stats.total_length,
            progress=self.verbose
        )

        # Estimate significance if requested
        if estimate_significance and scored_hits:
            params = significance_params or {}
            estimator = SignificanceEstimator(
                n_permutations=params.get('n_permutations', 100),
                random_seed=params.get('random_seed')
            )

            # Use existing background distribution if available
            if hasattr(self, '_background_distribution'):
                estimator.background_distribution = self._background_distribution

            # Estimate significance
            sign_results = estimator.estimate_significance_batch(
                hits=[s.hit for s in scored_hits],
                query_length=len(query_seq),
                db_size=self.index.stats.total_length,
                use_extreme_distribution=params.get('use_extreme', False),
                progress=self.verbose
            )

            # Merge significance into scored hits
            for scored, sign in zip(scored_hits, sign_results):
                scored.evalue = sign.evalue
                # Add significance attributes
                scored.is_significant = sign.is_significant
                scored.pvalue = sign.pvalue

        # Create result
        result = SearchResult(
            query_id=query_id,
            query_sequence=query_seq,
            hits=scored_hits,
            num_seeds_found=len(seeds),
            num_extensions=len(results),
            search_time=time.time() - start_time,
            parameters={
                'top_n': top_n,
                'min_score': min_score,
                'min_identity': min_identity,
                'max_seeds': max_seeds,
                'max_extensions': max_extensions,
                'estimate_significance': estimate_significance
            }
        )

        return result

    def search_fasta(
            self,
            fasta_file: Union[str, Path],
            top_n: Optional[int] = None,
            min_score: Optional[int] = None,
            min_identity: Optional[float] = None,
            max_queries: Optional[int] = None,
            estimate_significance: bool = False,
            significance_params: Optional[Dict[str, Any]] = None,
            progress: bool = True
    ) -> BatchSearchResult:
        """Search multiple queries from a FASTA file.

        Args:
            fasta_file: Path to FASTA file with queries.
            top_n: Number of top hits per query.
            min_score: Minimum score threshold.
            min_identity: Minimum identity percentage.
            max_queries: Maximum number of queries to process.
            estimate_significance: Whether to estimate significance.
            significance_params: Parameters for significance estimation.
            progress: Whether to show progress.

        Returns:
            BatchSearchResult object.

        Examples:
            >>> results = blast.search_fasta("queries.fa")
            >>> for result in results.results:
            ...     print(f"{result.query_id}: {result.num_hits} hits")
        """
        import time
        start_time = time.time()

        # Load queries
        queries = list(parse_fasta(fasta_file))
        if max_queries and len(queries) > max_queries:
            queries = queries[:max_queries]

        if not queries:
            raise ValueError(f"No queries found in {fasta_file}")

        if self.verbose:
            print(f"Processing {len(queries)} queries...", file=sys.stderr)

        results = []
        for i, query in enumerate(queries):
            if progress and self.verbose:
                print(f"  {i + 1}/{len(queries)}: {query.id}", file=sys.stderr)

            result = self.search(
                query=query,
                top_n=top_n,
                min_score=min_score,
                min_identity=min_identity,
                estimate_significance=estimate_significance,
                significance_params=significance_params
            )
            results.append(result)

        return BatchSearchResult(
            results=results,
            total_time=time.time() - start_time
        )

    def estimate_background(
            self,
            query: Union[str, SequenceRecord],
            n_permutations: int = 100,
            random_seed: Optional[int] = None
    ) -> 'BlastMini':
        """Estimate background distribution for significance testing.

        This caches the background distribution for use in future searches.

        Args:
            query: Query sequence or SequenceRecord.
            n_permutations: Number of permutations.
            random_seed: Random seed for reproducibility.

        Returns:
            Self (for chaining).

        Examples:
            >>> blast.estimate_background("ATCGATCGATCG", n_permutations=200)
            >>> results = blast.search("ATCGATCGATCG", estimate_significance=True)
        """
        estimator = SignificanceEstimator(
            n_permutations=n_permutations,
            random_seed=random_seed
        )

        self._background_distribution = estimator.estimate_background_distribution(
            query=query,
            subject_sequences=self.subject_sequences,
            n_permutations=n_permutations,
            progress=self.verbose
        )

        return self

    def get_stats(self) -> Dict[str, Any]:
        """Get statistics about the database and index.

        Returns:
            Dictionary of statistics.
        """
        stats = self.index.stats

        return {
            'kmer_size': stats.kmer_size,
            'total_sequences': stats.total_sequences,
            'total_length': stats.total_length,
            'unique_kmers': stats.unique_kmers,
            'total_kmers': stats.total_kmers,
            'memory_estimate_mb': stats.memory_estimate_mb,
            'subject_sequences_loaded': len(self.subject_sequences) > 0
        }

    def save_index(self, filepath: Union[str, Path]) -> None:
        """Save the index to a file.

        Args:
            filepath: Output file path.
        """
        self.index.save(filepath)
        if self.verbose:
            print(f"Index saved to {filepath}", file=sys.stderr)

    def set_verbose(self, verbose: bool) -> 'BlastMini':
        """Set verbose mode.

        Args:
            verbose: Whether to print progress messages.

        Returns:
            Self (for chaining).
        """
        self.verbose = verbose
        return self

    def set_config(self, config: AlignmentConfig) -> 'BlastMini':
        """Update configuration.

        Args:
            config: New configuration.

        Returns:
            Self (for chaining).
        """
        self.config = config
        self._init_components()
        return self

    def __repr__(self) -> str:
        return (f"BlastMini(k={self.index.k}, "
                f"sequences={self.index.stats.total_sequences}, "
                f"kmers={self.index.stats.unique_kmers:,})")


# ============================================================================
# Convenience Functions
# ============================================================================

def quick_search(
        query: Union[str, SequenceRecord],
        database: Union[str, Path],
        k: int = 11,
        top_n: int = 10,
        **kwargs
) -> SearchResult:
    """Quick one-line search.

    This is a convenience function for simple searches without creating
    a BlastMini instance.

    Args:
        query: Query sequence or SequenceRecord.
        database: Database FASTA file.
        k: k-mer size.
        top_n: Number of top hits.
        **kwargs: Additional arguments passed to BlastMini.search.

    Returns:
        SearchResult object.

    Examples:
        >>> result = quick_search("ATCGATCGATCG", "database.fa")
        >>> print(result.top_hit)
    """
    blast = BlastMini.from_fasta(database, k=k, verbose=False)
    return blast.search(query, top_n=top_n, **kwargs)


def quick_batch_search(
        queries_file: Union[str, Path],
        database: Union[str, Path],
        k: int = 11,
        top_n: int = 10,
        **kwargs
) -> BatchSearchResult:
    """Quick one-line batch search.

    Args:
        queries_file: FASTA file with queries.
        database: Database FASTA file.
        k: k-mer size.
        top_n: Number of top hits per query.
        **kwargs: Additional arguments passed to BlastMini.search_fasta.

    Returns:
        BatchSearchResult object.
    """
    blast = BlastMini.from_fasta(database, k=k, verbose=False)
    return blast.search_fasta(queries_file, top_n=top_n, **kwargs)


# ============================================================================
# Module Exports
# ============================================================================

__all__ = [
    'BlastMini',
    'SearchResult',
    'BatchSearchResult',
    'quick_search',
    'quick_batch_search',
    'ProgressCallback',
    'default_progress_callback'
]


# ============================================================================
# Command Line Interface for API Testing
# ============================================================================

def main():
    """Simple CLI for testing the API."""
    import argparse

    parser = argparse.ArgumentParser(description="blastmini API test")
    parser.add_argument("-d", "--database", required=True,
                        help="Database FASTA file")
    parser.add_argument("-q", "--query", required=True,
                        help="Query FASTA file or sequence string")
    parser.add_argument("-k", "--kmer", type=int, default=11,
                        help="k-mer size")
    parser.add_argument("--top", type=int, default=10,
                        help="Number of top hits")
    parser.add_argument("-o", "--output", help="Output file")

    args = parser.parse_args()

    # Create BlastMini instance
    print(f"Loading database from {args.database}...", file=sys.stderr)
    blast = BlastMini.from_fasta(args.database, k=args.kmer)

    print(f"Database stats: {blast.get_stats()}", file=sys.stderr)

    # Load or parse query
    try:
        queries = list(parse_fasta(args.query))
        if queries:
            query = queries[0]
        else:
            query = args.query
    except:
        query = args.query

    # Search
    print(f"Searching query...", file=sys.stderr)
    result = blast.search(query, top_n=args.top)

    print(f"\nFound {result.num_hits} hits in {result.search_time:.2f}s")
    print("\nTop hits:")
    print(result.to_text())

    if args.output:
        result.save(args.output)
        print(f"Results saved to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()