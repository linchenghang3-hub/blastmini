"""Seed extension module for blastmini.

This module implements the X-dropoff extension algorithm, which takes
seed matches and extends them bidirectionally to find high-scoring
segment pairs (HSPs). The extension phase is where the actual alignment
is generated from the seed positions.

Key features:
    - Bidirectional X-dropoff extension
    - Ungapped extension (simplified BLAST)
    - Support for both same-strand and reverse-complement extensions
    - Scoring with match/mismatch matrices
    - Extension statistics and quality tracking
    - Parallel extension support for multiple seeds
"""

import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple, Union

from .index import reverse_complement
from .models import AlignmentConfig, Hit, SequenceRecord
from .seeding import Seed, SeedCluster

# ============================================================================
# Extension Data Structures
# ============================================================================


@dataclass
class ExtensionResult:
    """Result of a seed extension operation.

    This class stores the complete information about an extended alignment,
    including the aligned sequences, scores, and coordinates.

    Attributes:
        query_id: Query sequence identifier.
        subject_id: Subject sequence identifier.
        query_start: Start position in query (0-based, inclusive).
        query_end: End position in query (0-based, exclusive).
        subject_start: Start position in subject (0-based, inclusive).
        subject_end: End position in subject (0-based, exclusive).
        score: Total alignment score.
        identity: Number of identical matches.
        mismatches: Number of mismatches.
        alignment_length: Total alignment length.
        query_alignment: Aligned query sequence (with gaps).
        subject_alignment: Aligned subject sequence (with gaps).
        strand: '+' for same strand, '-' for reverse complement.
        seed_origin: The seed that initiated this extension.
        max_score_during_extension: Peak score during extension.
        extension_steps: Number of steps taken during extension.

    Examples:
        >>> result = ExtensionResult("Q1", "S1", 10, 50, 100, 140)
        >>> print(f"Score: {result.score}, Identity: {result.identity}/{result.alignment_length}")
    """
    query_id: str
    subject_id: str
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    score: int = 0
    identity: int = 0
    mismatches: int = 0
    alignment_length: int = 0
    query_alignment: str = ""
    subject_alignment: str = ""
    strand: str = "+"
    seed_origin: Optional[Seed] = None
    max_score_during_extension: int = 0
    extension_steps: int = 0

    def __post_init__(self) -> None:
        """Calculate derived statistics."""
        if self.alignment_length == 0 and self.query_alignment:
            self.alignment_length = len(self.query_alignment)

    @property
    def identity_percent(self) -> float:
        """Calculate identity percentage."""
        if self.alignment_length == 0:
            return 0.0
        return (self.identity / self.alignment_length) * 100

    @property
    def query_sequence(self) -> str:
        """Get the query sequence without gaps."""
        return self.query_alignment.replace("-", "")

    @property
    def subject_sequence(self) -> str:
        """Get the subject sequence without gaps."""
        return self.subject_alignment.replace("-", "")

    def to_hit(self) -> Hit:
        """Convert extension result to a Hit object for output."""
        return Hit(
            query_id=self.query_id,
            subject_id=self.subject_id,
            score=self.score,
            identity_percent=self.identity_percent,
            alignment_length=self.alignment_length,
            query_start=self.query_start,
            query_end=self.query_end,
            subject_start=self.subject_start,
            subject_end=self.subject_end,
            query_alignment=self.query_alignment,
            subject_alignment=self.subject_alignment,
        )

    def __repr__(self) -> str:
        return (f"ExtensionResult(query='{self.query_id}', subject='{self.subject_id}', "
                f"score={self.score}, identity={self.identity_percent:.1f}%, "
                f"len={self.alignment_length}, strand={self.strand})")

    def __lt__(self, other: 'ExtensionResult') -> bool:
        """Sort by score descending."""
        return self.score > other.score


@dataclass
class ExtensionStats:
    """Statistics about the extension process.

    Tracks performance and quality metrics for extension operations.

    Attributes:
        total_extensions: Number of extension attempts.
        successful_extensions: Number of successful extensions.
        failed_extensions: Number of failed extensions.
        total_extension_steps: Total steps across all extensions.
        average_extension_length: Average length of successful extensions.
        max_extension_length: Maximum extension length achieved.
        min_score: Minimum score achieved.
        max_score: Maximum score achieved.
        total_time_ms: Total time spent in milliseconds.
    """
    total_extensions: int = 0
    successful_extensions: int = 0
    failed_extensions: int = 0
    total_extension_steps: int = 0
    average_extension_length: float = 0.0
    max_extension_length: int = 0
    min_score: int = 0
    max_score: int = 0
    total_time_ms: float = 0.0

    def record_extension(self, result: ExtensionResult) -> None:
        """Record a single extension result."""
        self.total_extensions += 1
        if result.alignment_length > 0:
            self.successful_extensions += 1
            self.total_extension_steps += result.extension_steps
            self.average_extension_length = (
                (self.average_extension_length * (self.successful_extensions - 1) +
                 result.alignment_length) / self.successful_extensions
            )
            self.max_extension_length = max(
                self.max_extension_length, result.alignment_length)
            self.min_score = min(
                self.min_score, result.score) if self.min_score != 0 else result.score
            self.max_score = max(self.max_score, result.score)
        else:
            self.failed_extensions += 1

    def __repr__(self) -> str:
        if self.total_extensions == 0:
            return "ExtensionStats(no extensions)"
        return (f"ExtensionStats(extensions={self.total_extensions}, "
                f"success={self.successful_extensions}, "
                f"avg_len={self.average_extension_length:.1f}, "
                f"max_score={self.max_score})")


# ============================================================================
# Core Extension Algorithms
# ============================================================================

class SeedExtender:
    """Extends seeds to generate high-scoring segment pairs (HSPs).

    This class implements the X-dropoff extension algorithm for both
    same-strand and reverse-complement alignments.

    Attributes:
        config: Alignment configuration parameters.
        match_score: Score for matching nucleotides.
        mismatch_penalty: Penalty for mismatching nucleotides.
        x_dropoff: X-dropoff threshold for extension termination.
        gap_open_penalty: Penalty for opening a gap (for future gapped extension).
        gap_extend_penalty: Penalty for extending a gap (for future gapped extension).
        min_extension_score: Minimum score for a valid extension.
        stats: Extension statistics tracker.

    Examples:
        >>> from blastmini.index import build_index_from_fasta
        >>> from blastmini.seeding import SeedFinder
        >>> idx = build_index_from_fasta("database.fa", k=11)
        >>> finder = SeedFinder(idx)
        >>> extender = SeedExtender()
        >>> query = SequenceRecord("Q1", "ATCGATCGATCG")
        >>> seeds = finder.find_seeds(query)
        >>> results = extender.extend_seeds(query, seeds)
    """

    def __init__(
            self,
            config: Optional[AlignmentConfig] = None,
            min_extension_score: int = 10,
            track_stats: bool = False
    ):
        """Initialize the seed extender.

        Args:
            config: Alignment configuration (uses defaults if None).
            min_extension_score: Minimum score for a valid extension.
            track_stats: Whether to track extension statistics.
        """
        self.config = config or AlignmentConfig()

        # Scoring parameters
        self.match_score = self.config.match_score
        self.mismatch_penalty = self.config.mismatch_penalty
        self.x_dropoff = self.config.x_dropoff
        self.gap_open_penalty = self.config.gap_open_penalty
        self.gap_extend_penalty = self.config.gap_extend_penalty

        self.min_extension_score = min_extension_score
        self.track_stats = track_stats
        self.stats = ExtensionStats() if track_stats else None

    def extend_seed(
            self,
            query: Union[SequenceRecord, str],
            seed: Seed,
            subject_sequences: Dict[str, str]
    ) -> Optional[ExtensionResult]:
        """Extend a single seed bidirectionally.

        Args:
            query: Query sequence or SequenceRecord.
            seed: Seed to extend.
            subject_sequences: Dictionary mapping subject ID to sequence.

        Returns:
            ExtensionResult if successful, None otherwise.
        """
        # Get query sequence
        if isinstance(query, SequenceRecord):
            query_id = query.id
            query_seq = query.sequence
        else:
            query_id = "query"
            query_seq = str(query)

        # Get subject sequence
        subject_seq = subject_sequences.get(seed.subject_id)
        if subject_seq is None:
            return None

        # Handle reverse complement if needed
        if seed.strand == "-":
            subject_seq = reverse_complement(subject_seq)
            # Adjust position for reverse complement
            subject_pos = len(subject_seq) - seed.subject_pos - 1
        else:
            subject_pos = seed.subject_pos

        # Extend bidirectionally
        result = self._extend_bidirectional(
            query_seq=query_seq,
            subject_seq=subject_seq,
            query_pos=seed.query_pos,
            subject_pos=subject_pos,
            query_id=query_id,
            subject_id=seed.subject_id,
            strand=seed.strand
        )

        if result and self.track_stats and self.stats:
            self.stats.record_extension(result)

        return result

    def extend_seeds(
            self,
            query: Union[SequenceRecord, str],
            seeds: List[Seed],
            subject_sequences: Dict[str, str],
            max_results: Optional[int] = None,
            progress: bool = True
    ) -> List[ExtensionResult]:
        """Extend multiple seeds and return the best results.

        Args:
            query: Query sequence or SequenceRecord.
            seeds: List of seeds to extend.
            subject_sequences: Dictionary mapping subject ID to sequence.
            max_results: Maximum number of results to return.
            progress: Show progress feedback.

        Returns:
            List of ExtensionResult objects sorted by score.
        """
        if not seeds:
            return []

        results: List[ExtensionResult] = []
        # subject_id, q_start, q_end, s_start, s_end
        seen_pairs: Set[Tuple[str, str, int, int]] = set()

        # Get query sequence
        if isinstance(query, SequenceRecord):
            query.sequence
        else:
            str(query)

        # Filter seeds by quality
        filtered_seeds = self._filter_seeds(seeds)

        total_seeds = len(filtered_seeds)
        if progress and total_seeds > 100:
            print(f"Extending {total_seeds} seeds...", file=sys.stderr)

        for i, seed in enumerate(filtered_seeds):
            # Extend the seed
            result = self.extend_seed(query, seed, subject_sequences)

            if result and result.score >= self.min_extension_score:
                # Avoid duplicate alignments
                key = (
                    result.subject_id,
                    result.query_start,
                    result.query_end,
                    result.subject_start,
                    result.subject_end
                )
                if key not in seen_pairs:
                    seen_pairs.add(key)
                    results.append(result)

            # Progress update
            if progress and total_seeds > 100 and (i + 1) % 10 == 0:
                print(f"  Extended {i + 1}/{total_seeds} seeds, found {len(results)} HSPs",
                      file=sys.stderr)

        # Sort by score descending
        results.sort(key=lambda r: r.score, reverse=True)

        if max_results and len(results) > max_results:
            results = results[:max_results]

        if progress and total_seeds > 100:
            print(f"  Done: found {len(results)} HSPs", file=sys.stderr)

        return results

    def extend_clusters(
            self,
            query: Union[SequenceRecord, str],
            clusters: List[SeedCluster],
            subject_sequences: Dict[str, str],
            max_results: Optional[int] = None
    ) -> List[ExtensionResult]:
        """Extend seed clusters and return the best results.

        This method extends seeds from clusters, giving priority to
        clusters with more seeds (higher confidence).

        Args:
            query: Query sequence or SequenceRecord.
            clusters: List of SeedCluster objects.
            subject_sequences: Dictionary mapping subject ID to sequence.
            max_results: Maximum number of results to return.

        Returns:
            List of ExtensionResult objects sorted by score.
        """
        if not clusters:
            return []

        results: List[ExtensionResult] = []
        seen_pairs: Set[Tuple[str, int, int, int, int]] = set()

        # Sort clusters by size (largest first)
        sorted_clusters = sorted(clusters, key=lambda c: len(c), reverse=True)

        for cluster in sorted_clusters:
            # Try best seed from cluster
            best_seed = cluster.max_score_seed()
            if best_seed is None:
                continue

            result = self.extend_seed(query, best_seed, subject_sequences)

            if result and result.score >= self.min_extension_score:
                key = (
                    result.subject_id,
                    result.query_start,
                    result.query_end,
                    result.subject_start,
                    result.subject_end
                )
                if key not in seen_pairs:
                    seen_pairs.add(key)
                    results.append(result)

            # If we have enough results, stop
            if max_results and len(results) >= max_results:
                break

        # Sort by score descending
        results.sort(key=lambda r: r.score, reverse=True)

        if max_results and len(results) > max_results:
            results = results[:max_results]

        return results

    def _extend_bidirectional(
            self,
            query_seq: str,
            subject_seq: str,
            query_pos: int,
            subject_pos: int,
            query_id: str,
            subject_id: str,
            strand: str = "+"
    ) -> Optional[ExtensionResult]:
        """Extend a seed bidirectionally using X-dropoff algorithm.

        This is the core extension algorithm that extends the seed
        to the left and right simultaneously.

        Args:
            query_seq: Full query sequence.
            subject_seq: Full subject sequence.
            query_pos: Position in query sequence.
            subject_pos: Position in subject sequence.
            query_id: Query identifier.
            subject_id: Subject identifier.
            strand: Strand orientation.

        Returns:
            ExtensionResult or None if extension failed.
        """
        # Initial seed score (all match)
        total_score = self.match_score * 1  # k-mer length is just the seed

        # Initialize alignment with seed
        query_start = query_pos
        query_end = query_pos + 1
        subject_start = subject_pos
        subject_end = subject_pos + 1

        # For tracking max score during extension
        max_score = total_score

        # Store alignment characters (for visualization)
        query_chars = [query_seq[query_pos]]
        subject_chars = [subject_seq[subject_pos]]

        # Extension steps
        extension_steps = 0

        # Extend to the right
        right_query = query_pos + 1
        right_subject = subject_pos + 1
        right_score = 0
        max_right_score = 0

        while (right_query < len(query_seq) and
               right_subject < len(subject_seq) and
               right_score - max_right_score >= -self.x_dropoff):

            q_char = query_seq[right_query]
            s_char = subject_seq[right_subject]

            # Score the current position
            if q_char == s_char:
                pos_score = self.match_score
                right_score += self.match_score
            else:
                pos_score = self.mismatch_penalty
                right_score += self.mismatch_penalty

            # Update max score for X-dropoff
            if right_score > max_right_score:
                max_right_score = right_score

            # Check if we should extend
            if right_score - max_right_score >= -self.x_dropoff:
                query_chars.append(q_char)
                subject_chars.append(s_char)
                total_score += pos_score
                extension_steps += 1

                # Update max score tracking
                if total_score > max_score:
                    max_score = total_score

                right_query += 1
                right_subject += 1
            else:
                break

        # Extend to the left
        left_query = query_pos - 1
        left_subject = subject_pos - 1
        left_score = 0
        max_left_score = 0

        left_chars = []
        left_subject_chars = []

        while (left_query >= 0 and
               left_subject >= 0 and
               left_score - max_left_score >= -self.x_dropoff):

            q_char = query_seq[left_query]
            s_char = subject_seq[left_subject]

            # Score the current position
            if q_char == s_char:
                pos_score = self.match_score
                left_score += self.match_score
            else:
                pos_score = self.mismatch_penalty
                left_score += self.mismatch_penalty

            # Update max score for X-dropoff
            if left_score > max_left_score:
                max_left_score = left_score

            # Check if we should extend
            if left_score - max_left_score >= -self.x_dropoff:
                left_chars.append(q_char)
                left_subject_chars.append(s_char)
                total_score += pos_score
                extension_steps += 1

                # Update max score tracking
                if total_score > max_score:
                    max_score = total_score

                left_query -= 1
                left_subject -= 1
            else:
                break

        # Check if extension was successful
        if extension_steps < 2:  # Need at least some extension
            return None

        # Build the final alignment
        # Left side (reverse order)
        left_query_aln = ''.join(reversed(left_chars))
        left_subject_aln = ''.join(reversed(left_subject_chars))

        # Seed (central part)
        seed_char = query_seq[query_pos]
        seed_subject_char = subject_seq[subject_pos]

        # Right side
        # Skip the seed char (included in query_chars[0])
        right_query_aln = ''.join(query_chars[1:])
        right_subject_aln = ''.join(subject_chars[1:])

        # Complete alignment
        query_alignment = left_query_aln + seed_char + right_query_aln
        subject_alignment = left_subject_aln + seed_subject_char + right_subject_aln

        # Calculate identity and mismatches
        identity = sum(1 for q, s in zip(query_alignment,
                       subject_alignment) if q == s and q != '-')
        mismatches = sum(1 for q, s in zip(query_alignment, subject_alignment)
                         if q != '-' and s != '-' and q != s)

        # Update coordinates
        query_start = left_query + 1
        query_end = right_query
        subject_start = left_subject + 1
        subject_end = right_subject

        # For reverse complement, the subject coordinates need to be adjusted
        # Since we already transformed the subject sequence, we need to map back
        if strand == "-":
            # The subject positions are in the transformed coordinate space
            # We need to map them back to the original coordinates
            subject_start = len(subject_seq) - subject_end
            subject_end = len(subject_seq) - subject_start

        # Create result
        result = ExtensionResult(
            query_id=query_id,
            subject_id=subject_id,
            query_start=query_start,
            query_end=query_end,
            subject_start=subject_start,
            subject_end=subject_end,
            score=total_score,
            identity=identity,
            mismatches=mismatches,
            alignment_length=len(query_alignment),
            query_alignment=query_alignment,
            subject_alignment=subject_alignment,
            strand=strand,
            max_score_during_extension=max_score,
            extension_steps=extension_steps
        )

        return result

    def _filter_seeds(self, seeds: List[Seed]) -> List[Seed]:
        """Filter seeds before extension to improve efficiency.

        Args:
            seeds: List of seeds to filter.

        Returns:
            Filtered list of seeds.
        """
        if not seeds:
            return seeds

        # Remove duplicate seeds (same query, subject, and positions)
        seen = set()
        filtered = []

        for seed in seeds:
            key = (seed.query_id, seed.subject_id,
                   seed.query_pos, seed.subject_pos)
            if key not in seen:
                seen.add(key)
                filtered.append(seed)

        # Sort by score descending (higher score first)
        filtered.sort(key=lambda s: s.score, reverse=True)

        # Keep only the best seeds per (query_id, subject_id)
        best_per_pair = defaultdict(list)
        for seed in filtered:
            key = (seed.query_id, seed.subject_id)
            best_per_pair[key].append(seed)

        # For each pair, keep top 5 seeds
        final_seeds = []
        for seeds_in_pair in best_per_pair.values():
            final_seeds.extend(seeds_in_pair[:5])

        return final_seeds

    def estimate_extension_time(self, num_seeds: int) -> float:
        """Estimate the time needed for extension.

        Args:
            num_seeds: Number of seeds to extend.

        Returns:
            Estimated time in seconds.
        """
        # Rough estimate: 0.001 seconds per seed extension
        return num_seeds * 0.001


# ============================================================================
# Extension Quality Assessment
# ============================================================================

def evaluate_extension_result(result: ExtensionResult) -> Dict[str, float]:
    """Evaluate the quality of an extension result.

    Args:
        result: ExtensionResult to evaluate.

    Returns:
        Dictionary with quality metrics.
    """
    metrics = {
        'score': float(result.score),
        'identity_percent': result.identity_percent,
        'length': float(result.alignment_length),
        'score_per_position': result.score / result.alignment_length if result.alignment_length > 0 else 0,
        'mismatch_rate': result.mismatches / result.alignment_length if result.alignment_length > 0 else 0,
        'extension_efficiency': result.extension_steps / result.alignment_length if result.alignment_length > 0 else 0,
        'score_stability': result.score / result.max_score_during_extension
        if result.max_score_during_extension > 0
        else 0
    }
    return metrics


def is_high_quality_extension(
        result: ExtensionResult,
        min_score: Optional[int] = None,
        min_identity: Optional[float] = None,
        min_length: Optional[int] = None
) -> bool:
    """Check if an extension meets quality thresholds.

    Args:
        result: ExtensionResult to check.
        min_score: Minimum score threshold.
        min_identity: Minimum identity percentage threshold.
        min_length: Minimum alignment length threshold.

    Returns:
        True if the extension meets all thresholds.
    """
    if min_score is not None and result.score < min_score:
        return False

    if min_identity is not None and result.identity_percent < min_identity:
        return False

    if min_length is not None and result.alignment_length < min_length:
        return False

    return True


# ============================================================================
# Command Line Interface (for testing)
# ============================================================================

def main():
    """Simple command line interface for testing the extension module."""
    import argparse

    parser = argparse.ArgumentParser(description="Seed extension tool")
    parser.add_argument("-d", "--database", required=True,
                        help="Database FASTA file")
    parser.add_argument("-q", "--query", required=True,
                        help="Query FASTA file")
    parser.add_argument("-k", "--kmer", type=int, default=11,
                        help="k-mer size")
    parser.add_argument("-s", "--seed", type=int, default=5,
                        help="Number of seeds to try")
    parser.add_argument("-o", "--output", help="Output file")
    parser.add_argument("--dropoff", type=int, default=5,
                        help="X-dropoff threshold")

    args = parser.parse_args()

    # Build index
    from .index import build_index_from_fasta
    from .io import parse_fasta
    from .seeding import SeedFinder

    print("Building index...", file=sys.stderr)
    idx = build_index_from_fasta(args.database, k=args.kmer, progress=False)

    # Create seed finder
    finder = SeedFinder(idx)

    # Load queries
    queries = list(parse_fasta(args.query))
    if not queries:
        print("No queries found", file=sys.stderr)
        sys.exit(1)

    # Load subject sequences
    subject_sequences = {}
    for record in parse_fasta(args.database):
        subject_sequences[record.id] = record.sequence

    # Create extender with custom dropoff
    config = AlignmentConfig(x_dropoff=args.dropoff)
    extender = SeedExtender(config=config, track_stats=True)

    print(f"Processing {len(queries)} queries...", file=sys.stderr)

    for query in queries:
        print(f"\nQuery: {query.id}")

        # Find seeds
        seeds = finder.find_best_seeds(query, top_n=args.seed)
        print(f"  Found {len(seeds)} seeds")

        if not seeds:
            print("  No seeds found")
            continue

        # Extend seeds
        results = extender.extend_seeds(query, seeds, subject_sequences)

        print(f"  Found {len(results)} HSPs")

        # Show top results
        for i, result in enumerate(results[:5]):
            print(f"    {i + 1}: {result.subject_id} "
                  f"score={result.score}, identity={result.identity_percent:.1f}%, "
                  f"len={result.alignment_length}")

    if extender.stats:
        print("\nExtension Statistics:", file=sys.stderr)
        print(
            f"  Total extensions: {extender.stats.total_extensions}", file=sys.stderr)
        print(
            f"  Successful: {extender.stats.successful_extensions}", file=sys.stderr)
        print(
            f"  Average length: {extender.stats.average_extension_length:.1f}", file=sys.stderr)
        print(f"  Max score: {extender.stats.max_score}", file=sys.stderr)


if __name__ == "__main__":
    main()
