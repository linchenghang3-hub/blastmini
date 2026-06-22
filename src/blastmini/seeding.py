"""Seed search module for blastmini.

This module implements the seeding phase of the BLAST algorithm, which is
responsible for finding initial matches (seeds) between query sequences
and the indexed database. The seeding phase is critical for efficiency
as it reduces the search space from O(mn) to O(m * k_mer_matches).

Key features:
    - Query k-mer extraction with position tracking
    - Dictionary-based seed finding using the k-mer index
    - Seed clustering (merge nearby seeds)
    - Seed filtering based on score and position
    - Support for single and multiple query sequences
    - Progress tracking for large queries
"""

import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, Union

from .index import KmerIndex, extract_kmers
from .models import AlignmentConfig, SequenceRecord

# ============================================================================
# Seed Data Structures
# ============================================================================


@dataclass
class Seed:
    """A single seed match between query and subject sequences.

    A seed represents a k-mer match that serves as a starting point
    for the extension phase.

    Attributes:
        query_id: Identifier of the query sequence.
        subject_id: Identifier of the subject sequence.
        query_pos: Position in the query sequence (0-based).
        subject_pos: Position in the subject sequence (0-based).
        kmer: The k-mer that formed this seed.
        score: Initial seed score (default 0, can be used for filtering).
        strand: '+' for same strand, '-' for reverse complement.

    Examples:
        >>> seed = Seed("Q1", "S1", 10, 100, "ATCGAT")
        >>> print(seed)
        Seed(Q1 -> S1 at query=10, subject=100)
    """
    query_id: str
    subject_id: str
    query_pos: int
    subject_pos: int
    kmer: str = ""
    score: int = 0
    strand: str = "+"

    def __repr__(self) -> str:
        return (f"Seed(query='{self.query_id}', subject='{self.subject_id}', "
                f"q_pos={self.query_pos}, s_pos={self.subject_pos}, "
                f"strand={self.strand})")

    def distance_to(self, other: 'Seed') -> int:
        """Calculate the distance between two seeds on the query axis.

        Args:
            other: Another seed to compare with.

        Returns:
            Absolute distance in query positions.
        """
        return abs(self.query_pos - other.query_pos)

    def is_on_same_diagonal(self, other: 'Seed', max_offset: int = 10) -> bool:
        """Check if two seeds lie on the same diagonal.

        Two seeds are on the same diagonal if the difference between
        their query and subject positions is similar (allowing some offset).

        Args:
            other: Another seed to compare with.
            max_offset: Maximum allowed offset in diagonal difference.

        Returns:
            True if seeds are on the same diagonal.
        """
        diag1 = self.subject_pos - self.query_pos
        diag2 = other.subject_pos - other.query_pos
        return abs(diag1 - diag2) <= max_offset


@dataclass
class SeedCluster:
    """A cluster of nearby seeds on the same diagonal.

    Clustering seeds helps identify regions with multiple matches,
    which increases confidence and improves extension quality.

    Attributes:
        query_id: Query sequence ID.
        subject_id: Subject sequence ID.
        seeds: List of seeds in the cluster.
        strand: '+' or '-' for the strand orientation.
        query_range: (start, end) of the cluster on the query.
        subject_range: (start, end) of the cluster on the subject.
        total_score: Sum of all seed scores.
        density: Number of seeds per unit length.

    Examples:
        >>> cluster = SeedCluster("Q1", "S1")
        >>> cluster.add_seed(Seed("Q1", "S1", 10, 100))
        >>> cluster.add_seed(Seed("Q1", "S1", 15, 105))
        >>> print(cluster)
        SeedCluster(Q1 -> S1: 2 seeds, range q=10-15, s=100-105)
    """
    query_id: str
    subject_id: str
    seeds: List[Seed] = field(default_factory=list)
    strand: str = "+"

    def __post_init__(self):
        if not self.seeds:
            return
        # Sort seeds by query position
        self.seeds.sort(key=lambda s: s.query_pos)
        self._update_ranges()

    def add_seed(self, seed: Seed) -> None:
        """Add a seed to the cluster."""
        self.seeds.append(seed)
        self.seeds.sort(key=lambda s: s.query_pos)
        self._update_ranges()

    def _update_ranges(self) -> None:
        """Update the query and subject ranges."""
        if not self.seeds:
            return
        self.query_range = (self.seeds[0].query_pos, self.seeds[-1].query_pos)
        # Subject range may be non-contiguous, but we use min/max
        subject_positions = [s.subject_pos for s in self.seeds]
        self.subject_range = (min(subject_positions), max(subject_positions))

        # Calculate density: seeds per 100 bases
        q_len = self.query_range[1] - self.query_range[0] + 1
        if q_len > 0:
            self.density = len(self.seeds) / q_len * 100
        else:
            self.density = 0

        # Total score (assuming each seed has equal weight for now)
        self.total_score = len(self.seeds)

    def __len__(self) -> int:
        return len(self.seeds)

    def __repr__(self) -> str:
        if not self.seeds:
            return f"SeedCluster({self.query_id} -> {self.subject_id}: empty)"
        return (f"SeedCluster({self.query_id} -> {self.subject_id}: "
                f"{len(self.seeds)} seeds, q={self.query_range}, "
                f"s={self.subject_range}, density={self.density:.1f})")

    def max_score_seed(self) -> Optional[Seed]:
        """Return the seed with the highest score in the cluster."""
        if not self.seeds:
            return None
        return max(self.seeds, key=lambda s: s.score)


# ============================================================================
# Seed Finder
# ============================================================================

class SeedFinder:
    """Finds seeds by matching query k-mers against the index.

    This class implements the core seeding algorithm, including k-mer
    extraction, dictionary lookup, and seed clustering.

    Attributes:
        index: KmerIndex object containing the database.
        config: AlignmentConfig with parameters.
        canonical: Whether to use canonical k-mers.
        min_seed_score: Minimum seed score to keep.
        max_seed_cluster_gap: Maximum gap to cluster seeds together.

    Examples:
        >>> from blastmini.index import KmerIndex
        >>> from blastmini.io import parse_fasta
        >>> records = list(parse_fasta("database.fa"))
        >>> idx = KmerIndex.build(records, k=11)
        >>> finder = SeedFinder(idx)
        >>> seeds = finder.find_seeds(records[0])
        >>> print(f"Found {len(seeds)} seeds")
    """

    def __init__(
            self,
            index: KmerIndex,
            config: Optional[AlignmentConfig] = None,
            canonical: bool = True,
            min_seed_score: int = 1,
            max_seed_cluster_gap: int = 20,
            max_seeds_per_region: int = 100
    ):
        """Initialize the seed finder.

        Args:
            index: K-mer index to search against.
            config: Alignment configuration (uses defaults if None).
            canonical: Use canonical k-mers.
            min_seed_score: Minimum score for a seed (default 1).
            max_seed_cluster_gap: Max query gap for clustering.
            max_seeds_per_region: Max seeds to return per query region.
        """
        self.index = index
        self.config = config or AlignmentConfig()
        self.canonical = canonical
        self.min_seed_score = min_seed_score
        self.max_seed_cluster_gap = max_seed_cluster_gap
        self.max_seeds_per_region = max_seeds_per_region

        # Cache for query k-mer lookups
        self._kmer_cache: Dict[str, List[Tuple[str, int]]] = {}

    def find_seeds(
            self,
            query: Union[SequenceRecord, str],
            subject_ids: Optional[List[str]] = None,
            max_seeds: Optional[int] = None,
            progress: bool = True
    ) -> List[Seed]:
        """Find seeds for a single query sequence.

        Args:
            query: Query sequence or SequenceRecord.
            subject_ids: Optional list of subject IDs to search (all if None).
            max_seeds: Maximum number of seeds to return.
            progress: Show progress feedback.

        Returns:
            List of Seed objects.

        Examples:
            >>> seeds = finder.find_seeds("ATCGATCGATCG")
            >>> for seed in seeds[:5]:
            ...     print(seed)
        """
        # Extract query ID and sequence
        if isinstance(query, SequenceRecord):
            query_id = query.id
            query_seq = query.sequence
        else:
            query_id = "query"
            query_seq = str(query)

        # Extract k-mers from query
        query_kmers = list(extract_kmers(
            query_seq,
            self.index.k,
            canonical=self.canonical
        ))

        if not query_kmers:
            return []

        # Find seeds
        seeds: List[Seed] = []
        # subject_id, query_pos, subject_pos
        seen_positions: Set[Tuple[str, int, int]] = set()

        # Determine which subject IDs to search
        if subject_ids is None:
            subject_ids = self.index.subject_ids

        subject_id_set = set(subject_ids)

        # Progress tracking
        total_kmers = len(query_kmers)
        if progress and total_kmers > 1000:
            print(
                f"Searching {total_kmers} k-mers in index...", file=sys.stderr)

        for idx, (kmer, q_pos) in enumerate(query_kmers):
            # Look up k-mer in index
            matches = self.index.lookup(
                kmer, canonical=False)  # Already canonical

            if not matches:
                continue

            # Add seeds for each match
            for subject_id, s_pos in matches:
                if subject_id not in subject_id_set:
                    continue

                # Avoid duplicate positions
                key = (subject_id, q_pos, s_pos)
                if key in seen_positions:
                    continue
                seen_positions.add(key)

                # Determine strand (simplified: if using canonical, strand is ambiguous)
                # For now, assume same strand
                strand = "+"

                seed = Seed(
                    query_id=query_id,
                    subject_id=subject_id,
                    query_pos=q_pos,
                    subject_pos=s_pos,
                    kmer=kmer,
                    score=1,  # Each seed has weight 1
                    strand=strand
                )
                seeds.append(seed)

                # Check max seeds
                if max_seeds and len(seeds) >= max_seeds:
                    if progress:
                        print(
                            f"  Reached max seeds ({max_seeds})", file=sys.stderr)
                    return self._filter_seeds(seeds)

            # Progress update
            if progress and total_kmers > 1000 and (idx + 1) % 100 == 0:
                print(f"  Processed {idx + 1}/{total_kmers} k-mers, found {len(seeds)} seeds",
                      file=sys.stderr)

        if progress and total_kmers > 1000:
            print(f"  Done: found {len(seeds)} seeds", file=sys.stderr)

        return self._filter_seeds(seeds)

    def find_seeds_batch(
            self,
            queries: List[Union[SequenceRecord, str]],
            subject_ids: Optional[List[str]] = None,
            max_seeds_per_query: Optional[int] = None,
            progress: bool = True
    ) -> Dict[str, List[Seed]]:
        """Find seeds for multiple query sequences.

        Args:
            queries: List of query sequences or SequenceRecords.
            subject_ids: Optional list of subject IDs to search.
            max_seeds_per_query: Maximum seeds per query.
            progress: Show progress feedback.

        Returns:
            Dictionary mapping query ID to list of seeds.
        """
        results = {}

        total_queries = len(queries)
        if progress:
            print(f"Processing {total_queries} queries...", file=sys.stderr)

        for i, query in enumerate(queries):
            query_id = query.id if isinstance(
                query, SequenceRecord) else f"query_{i}"
            results[query_id] = self.find_seeds(
                query=query,
                subject_ids=subject_ids,
                max_seeds=max_seeds_per_query,
                progress=False  # Don't show per-query progress in batch
            )

            if progress and (i + 1) % 10 == 0:
                print(
                    f"  Processed {i + 1}/{total_queries} queries", file=sys.stderr)

        return results

    def _filter_seeds(self, seeds: List[Seed]) -> List[Seed]:
        """Filter seeds based on quality thresholds.

        Args:
            seeds: List of seeds to filter.

        Returns:
            Filtered list of seeds.
        """
        if not seeds:
            return seeds

        # Sort by score (higher is better)
        seeds.sort(key=lambda s: s.score, reverse=True)

        # Apply max seeds per region filter
        if len(seeds) > self.max_seeds_per_region:
            # Keep only top seeds per (query_id, subject_id)
            filtered = []
            seen_subjects = defaultdict(int)

            for seed in seeds:
                key = (seed.query_id, seed.subject_id)
                if seen_subjects[key] < self.max_seeds_per_region:
                    filtered.append(seed)
                    seen_subjects[key] += 1

            return filtered

        return seeds

    def cluster_seeds(self, seeds: List[Seed]) -> List[SeedCluster]:
        """Cluster seeds by subject and diagonal.

        Seeds on the same diagonal within the same subject are clustered
        together, which helps identify high-confidence regions.

        Args:
            seeds: List of seeds to cluster.

        Returns:
            List of SeedCluster objects.
        """
        if not seeds:
            return []

        # Group seeds by (subject_id, strand)
        groups: Dict[Tuple[str, str], List[Seed]] = defaultdict(list)
        for seed in seeds:
            groups[(seed.subject_id, seed.strand)].append(seed)

        clusters = []

        for (subject_id, strand), group_seeds in groups.items():
            # Sort by query position
            group_seeds.sort(key=lambda s: s.query_pos)

            # Cluster by diagonal offset
            current_cluster = None

            for seed in group_seeds:
                if current_cluster is None:
                    # Start new cluster
                    current_cluster = SeedCluster(
                        query_id=seed.query_id,
                        subject_id=subject_id,
                        strand=strand
                    )
                    current_cluster.add_seed(seed)
                    clusters.append(current_cluster)
                else:
                    # Check if seed belongs to current cluster
                    last_seed = current_cluster.seeds[-1]

                    # Check if on same diagonal (allow small offset)
                    diag1 = last_seed.subject_pos - last_seed.query_pos
                    diag2 = seed.subject_pos - seed.query_pos
                    diag_diff = abs(diag1 - diag2)

                    # Check gap in query position
                    q_gap = seed.query_pos - last_seed.query_pos

                    if (diag_diff <= 10 and q_gap <= self.max_seed_cluster_gap):
                        # Add to current cluster
                        current_cluster.add_seed(seed)
                    else:
                        # Start new cluster
                        current_cluster = SeedCluster(
                            query_id=seed.query_id,
                            subject_id=subject_id,
                            strand=strand
                        )
                        current_cluster.add_seed(seed)
                        clusters.append(current_cluster)

        # Sort clusters by size (largest first)
        clusters.sort(key=lambda c: len(c), reverse=True)

        return clusters

    def find_best_seeds(
            self,
            query: Union[SequenceRecord, str],
            top_n: int = 10,
            subject_ids: Optional[List[str]] = None
    ) -> List[Seed]:
        """Find the best seeds for a query, prioritizing high-density regions.

        This method finds seeds, clusters them, and returns the top seeds
        from the best clusters.

        Args:
            query: Query sequence or SequenceRecord.
            top_n: Number of top seeds to return.
            subject_ids: Optional subject IDs to search.

        Returns:
            List of top seeds.
        """
        # Find all seeds
        seeds = self.find_seeds(
            query=query,
            subject_ids=subject_ids,
            max_seeds=self.max_seeds_per_region * 10,  # Allow more for clustering
            progress=False
        )

        if not seeds:
            return []

        # Cluster seeds
        clusters = self.cluster_seeds(seeds)

        if not clusters:
            return seeds[:top_n]

        # Select top seeds from top clusters
        selected_seeds = []
        for cluster in clusters:
            # Take up to 2 seeds per cluster
            cluster_seeds = cluster.seeds[:2]
            selected_seeds.extend(cluster_seeds)

            if len(selected_seeds) >= top_n:
                break

        return selected_seeds[:top_n]


# ============================================================================
# Seed Filtering Utilities
# ============================================================================

def filter_seeds_by_density(
        seeds: List[Seed],
        min_density: float = 1.0,
        window_size: int = 100
) -> List[Seed]:
    """Filter seeds to keep only those in high-density regions.

    Args:
        seeds: List of seeds.
        min_density: Minimum density threshold (seeds per 100 bases).
        window_size: Window size for density calculation.

    Returns:
        Filtered list of seeds.
    """
    if not seeds:
        return seeds

    # Group by subject
    by_subject: Dict[str, List[Seed]] = defaultdict(list)
    for seed in seeds:
        by_subject[seed.subject_id].append(seed)

    filtered = []

    for subject_id, subject_seeds in by_subject.items():
        # Sort by query position
        subject_seeds.sort(key=lambda s: s.query_pos)

        # Calculate density in sliding window
        for i, seed in enumerate(subject_seeds):
            # Count seeds within window_size
            window_start = seed.query_pos - window_size // 2
            window_end = seed.query_pos + window_size // 2

            count = sum(
                1 for s in subject_seeds
                if window_start <= s.query_pos <= window_end
            )

            # Calculate density per 100 bases
            density = count / window_size * 100

            if density >= min_density:
                filtered.append(seed)

    return filtered


def filter_seeds_by_coverage(
        seeds: List[Seed],
        max_coverage: int = 10
) -> List[Seed]:
    """Filter seeds to avoid over-sampling highly covered regions.

    Args:
        seeds: List of seeds.
        max_coverage: Maximum number of seeds per query position.

    Returns:
        Filtered list of seeds.
    """
    if not seeds:
        return seeds

    # Count seeds per query position
    position_counts: Dict[int, int] = defaultdict(int)
    for seed in seeds:
        position_counts[seed.query_pos] += 1

    # Keep only seeds from positions with count <= max_coverage
    filtered = []
    for seed in seeds:
        if position_counts[seed.query_pos] <= max_coverage:
            filtered.append(seed)

    # If filtering removed all seeds, return original (avoid empty results)
    if not filtered:
        return seeds[:len(seeds) // 2]  # Keep half as fallback

    return filtered


# ============================================================================
# Command Line Interface (for testing)
# ============================================================================

def main():
    """Simple command line interface for testing the seeding module."""
    import argparse

    parser = argparse.ArgumentParser(description="Seed search tool")
    parser.add_argument("action", choices=["search", "cluster"],
                        help="Action to perform")
    parser.add_argument("-d", "--database", required=True,
                        help="Database FASTA file or index file")
    parser.add_argument("-q", "--query", required=True,
                        help="Query FASTA file")
    parser.add_argument("-k", "--kmer", type=int, default=11,
                        help="k-mer size")
    parser.add_argument("-o", "--output", help="Output file")
    parser.add_argument("--top", type=int, default=10,
                        help="Number of top seeds to show")

    args = parser.parse_args()

    # Load or build index
    try:
        import os

        from .index import build_index_from_fasta, load_index_from_file

        if os.path.exists(args.database) and args.database.endswith('.json'):
            idx = load_index_from_file(args.database)
        else:
            idx = build_index_from_fasta(args.database, k=args.kmer)
    except Exception as e:
        print(f"Error loading index: {e}", file=sys.stderr)
        sys.exit(1)

    # Load query
    from .io import parse_fasta
    queries = list(parse_fasta(args.query))

    if not queries:
        print("No queries found", file=sys.stderr)
        sys.exit(1)

    # Create seed finder
    finder = SeedFinder(idx)

    if args.action == "search":
        print(
            f"Searching seeds for {len(queries)} queries...", file=sys.stderr)

        for query in queries:
            seeds = finder.find_best_seeds(query, top_n=args.top)
            print(f"\nQuery: {query.id} - Found {len(seeds)} seeds")
            for i, seed in enumerate(seeds[:args.top]):
                print(
                    f"  {i + 1}: {seed.subject_id} at q={seed.query_pos}, s={seed.subject_pos}")

    elif args.action == "cluster":
        for query in queries:
            seeds = finder.find_seeds(query, max_seeds=1000)
            clusters = finder.cluster_seeds(seeds)

            print(f"\nQuery: {query.id}")
            print(f"  Total seeds: {len(seeds)}")
            print(f"  Total clusters: {len(clusters)}")

            for i, cluster in enumerate(clusters[:5]):
                print(f"  Cluster {i + 1}: {len(cluster)} seeds, "
                      f"q={cluster.query_range}, s={cluster.subject_range}")


if __name__ == "__main__":
    main()
