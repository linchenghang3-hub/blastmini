"""k-mer index construction and management for blastmini.

This module provides hash-based k-mer indexing for nucleotide sequences,
which is the foundation of the seed-and-extend search strategy. The index
maps each k-mer to all positions where it appears in the database, enabling
O(1) lookup during seed search.

Key features:
    - Canonical k-mer handling (considers reverse complement)
    - Memory-efficient indexing with position tracking
    - Index persistence (save/load JSON)
    - Support for multiple sequences in a single index
    - k-mer frequency analysis
"""

import json
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Tuple, Union

from .io import load_index, parse_fasta, save_index
from .models import SequenceRecord

# ============================================================================
# DNA Utilities
# ============================================================================


def reverse_complement(seq: str) -> str:
    """Compute the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence string (A, T, C, G, N).

    Returns:
        Reverse complement sequence.

    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("N")
        'N'
    """
    complement_map = {
        'A': 'T', 'T': 'A',
        'C': 'G', 'G': 'C',
        'a': 'T', 't': 'A',
        'c': 'G', 'g': 'C',
        'N': 'N', 'n': 'N'
    }
    return ''.join(complement_map.get(base, 'N') for base in reversed(seq))


def canonical_kmer(kmer: str) -> str:
    """Return the canonical representation of a k-mer.

    The canonical k-mer is the lexicographically smaller of the k-mer
    and its reverse complement. This ensures that both strands of DNA
    map to the same index key, improving sensitivity.

    Args:
        kmer: DNA k-mer string.

    Returns:
        Canonical k-mer string.

    Examples:
        >>> canonical_kmer("ATCG")
        'ATCG'
        >>> canonical_kmer("CGAT")  # reverse complement of ATCG is ATCG
        'ATCG'
        >>> canonical_kmer("AAAA")
        'AAAA'
    """
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc


def extract_kmers(sequence: str, k: int, canonical: bool = True) -> Iterator[Tuple[str, int]]:
    """Extract all k-mers from a sequence with their positions.

    Args:
        sequence: DNA sequence string.
        k: k-mer length.
        canonical: If True, convert to canonical k-mers.

    Yields:
        Tuple of (kmer_string, position) for each k-mer.

    Examples:
        >>> list(extract_kmers("ATCGATCG", 3))
        [('ATC', 0), ('TCG', 1), ('CGA', 2), ('GAT', 3), ('ATC', 4)]
    """
    sequence_upper = sequence.upper()
    seq_len = len(sequence_upper)

    for i in range(seq_len - k + 1):
        kmer = sequence_upper[i:i + k]
        if canonical:
            kmer = canonical_kmer(kmer)
        yield (kmer, i)


# ============================================================================
# Index Data Structures
# ============================================================================

@dataclass
class IndexStats:
    """Statistics about a k-mer index.

    Attributes:
        total_kmers: Total number of k-mer occurrences.
        unique_kmers: Number of unique k-mers.
        total_sequences: Number of sequences in the database.
        total_length: Total sequence length.
        kmer_size: Size of k-mers used.
        memory_estimate_mb: Estimated memory usage in megabytes.
    """
    total_kmers: int = 0
    unique_kmers: int = 0
    total_sequences: int = 0
    total_length: int = 0
    kmer_size: int = 0
    memory_estimate_mb: float = 0.0

    def __repr__(self) -> str:
        return (f"IndexStats(kmers={self.total_kmers:,}, "
                f"unique={self.unique_kmers:,}, "
                f"sequences={self.total_sequences}, "
                f"k={self.kmer_size}, "
                f"memory={self.memory_estimate_mb:.2f}MB)")


class KmerIndex:
    """Hash-based k-mer index for nucleotide sequences.

    This class manages the construction, querying, and persistence of
    k-mer indices. It supports canonical k-mers for strand-independent
    searching and provides statistics for understanding index properties.

    Attributes:
        index: Dictionary mapping k-mer string to list of (subject_id, position).
        k: k-mer size used for the index.
        stats: Index statistics.
        subject_ids: List of sequence IDs in the database.

    Examples:
        >>> from blastmini.io import parse_fasta
        >>> records = list(parse_fasta("database.fa"))
        >>> index = KmerIndex.build(records, k=11)
        >>> matches = index.lookup("ATCGATCGATC")
        >>> print(f"Found {len(matches)} matches")
    """

    def __init__(self, k: int = 11):
        """Initialize an empty k-mer index.

        Args:
            k: k-mer size (default 11 for DNA BLAST).
        """
        if k < 1:
            raise ValueError(f"k must be >= 1, got {k}")
        self.k = k
        self.index: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
        self.subject_ids: List[str] = []
        self.stats = IndexStats(kmer_size=k)

    @classmethod
    def build(
            cls,
            records: List[SequenceRecord],
            k: int = 11,
            canonical: bool = True,
            min_kmer_occurrences: int = 1,
            max_kmer_occurrences: Optional[int] = None,
            progress: bool = True
    ) -> 'KmerIndex':
        """Build a k-mer index from a list of sequence records.

        This is the primary factory method for creating indices. It handles
        filtering of repetitive k-mers and provides progress feedback.

        Args:
            records: List of SequenceRecord objects to index.
            k: k-mer size.
            canonical: Use canonical k-mers (consider reverse complement).
            min_kmer_occurrences: Minimum occurrences to keep (filter rare).
            max_kmer_occurrences: Maximum occurrences to keep (filter repetitive).
            progress: Show progress bar (if tqdm is installed).

        Returns:
            KmerIndex object with populated index.

        Examples:
            >>> from blastmini.models import SequenceRecord
            >>> records = [SequenceRecord("S1", "ATCGATCG"), SequenceRecord("S2", "GATCGATC")]
            >>> idx = KmerIndex.build(records, k=3)
            >>> print(idx.lookup("ATC"))
            [('S1', 0), ('S1', 4), ('S2', 2)]
        """
        index_obj = cls(k)

        # Collect all sequence IDs
        index_obj.subject_ids = [rec.id for rec in records]

        # Build frequency dictionary first to filter if needed
        if min_kmer_occurrences > 1 or max_kmer_occurrences is not None:
            freq = {}
            for rec in records:
                seq = rec.sequence.upper()
                for i in range(len(seq) - k + 1):
                    kmer = seq[i:i + k]
                    if canonical:
                        kmer = canonical_kmer(kmer)
                    freq[kmer] = freq.get(kmer, 0) + 1

            # Build filter set
            keep_kmers = {
                kmer for kmer, count in freq.items()
                if count >= min_kmer_occurrences
                and (max_kmer_occurrences is None or count <= max_kmer_occurrences)
            }
        else:
            keep_kmers = None

        # Build the actual index
        total_kmers = 0
        total_length = 0
        seen_kmers = set()

        # Progress indicator (simple)
        total_records = len(records)
        if progress and total_records > 0:
            print(
                f"Indexing {total_records} sequences with k={k}...", file=sys.stderr)

        for idx, rec in enumerate(records):
            seq = rec.sequence.upper()
            total_length += len(seq)

            for pos in range(len(seq) - k + 1):
                kmer = seq[pos:pos + k]
                if canonical:
                    kmer = canonical_kmer(kmer)

                # Apply occurrence filter
                if keep_kmers is not None and kmer not in keep_kmers:
                    continue

                index_obj.index[kmer].append((rec.id, pos))
                total_kmers += 1
                seen_kmers.add(kmer)

            if progress and (idx + 1) % 100 == 0:
                print(
                    f"  Processed {idx + 1}/{total_records} sequences...", file=sys.stderr)

        # Update statistics
        index_obj.stats = IndexStats(
            total_kmers=total_kmers,
            unique_kmers=len(seen_kmers),
            total_sequences=len(records),
            total_length=total_length,
            kmer_size=k,
            memory_estimate_mb=index_obj._estimate_memory()
        )

        if progress:
            print(f"  Done: {index_obj.stats}", file=sys.stderr)

        return index_obj

    @classmethod
    def build_from_fasta(
            cls,
            fasta_file: Union[str, os.PathLike],
            k: int = 11,
            canonical: bool = True,
            min_kmer_occurrences: int = 1,
            max_kmer_occurrences: Optional[int] = None,
            max_records: Optional[int] = None,
            progress: bool = True
    ) -> 'KmerIndex':
        """Build a k-mer index directly from a FASTA file.

        Convenience method that parses the FASTA file and builds the index
        in one step. For large files, this uses streaming but still loads
        all records into memory.

        Args:
            fasta_file: Path to FASTA file.
            k: k-mer size.
            canonical: Use canonical k-mers.
            min_kmer_occurrences: Filter rare k-mers.
            max_kmer_occurrences: Filter repetitive k-mers.
            max_records: Maximum number of records to index (for testing).
            progress: Show progress feedback.

        Returns:
            KmerIndex object.

        Examples:
            >>> idx = KmerIndex.build_from_fasta("database.fa", k=11, max_records=100)
        """
        records = []
        for rec in parse_fasta(fasta_file):
            records.append(rec)
            if max_records and len(records) >= max_records:
                break

        return cls.build(
            records=records,
            k=k,
            canonical=canonical,
            min_kmer_occurrences=min_kmer_occurrences,
            max_kmer_occurrences=max_kmer_occurrences,
            progress=progress
        )

    def lookup(self, kmer: str, canonical: bool = True) -> List[Tuple[str, int]]:
        """Look up a k-mer in the index.

        Args:
            kmer: k-mer string to look up.
            canonical: If True, lookup the canonical version.

        Returns:
            List of (subject_id, position) tuples for the k-mer.
            Returns empty list if k-mer is not found.

        Examples:
            >>> idx = KmerIndex.build(records, k=3)
            >>> idx.lookup("ATC")
            [('S1', 0), ('S1', 4)]
        """
        if len(kmer) != self.k:
            # If k-mer length doesn't match, we can still try to lookup
            # but this should be handled by the caller
            pass

        if canonical:
            kmer = canonical_kmer(kmer)

        return self.index.get(kmer, [])

    def lookup_batch(self, kmers: List[str], canonical: bool = True) -> Dict[str, List[Tuple[str, int]]]:
        """Look up multiple k-mers efficiently.

        Args:
            kmers: List of k-mer strings.
            canonical: Use canonical k-mers.

        Returns:
            Dictionary mapping k-mer string to list of matches.
        """
        results = {}
        for kmer in kmers:
            if canonical:
                kmer_key = canonical_kmer(kmer)
            else:
                kmer_key = kmer
            results[kmer_key] = self.index.get(kmer_key, [])
        return results

    def contains(self, kmer: str, canonical: bool = True) -> bool:
        """Check if a k-mer exists in the index.

        Args:
            kmer: k-mer string.
            canonical: Use canonical k-mers.

        Returns:
            True if the k-mer is in the index.
        """
        if canonical:
            kmer = canonical_kmer(kmer)
        return kmer in self.index

    def get_kmers_with_counts(
            self,
            min_count: int = 1,
            max_count: Optional[int] = None
    ) -> List[Tuple[str, int]]:
        """Get k-mers with their occurrence counts.

        Useful for analyzing k-mer distributions and filtering.

        Args:
            min_count: Minimum count threshold.
            max_count: Maximum count threshold.

        Returns:
            List of (kmer, count) tuples sorted by count descending.
        """
        counts = [(kmer, len(positions))
                  for kmer, positions in self.index.items()]

        # Apply filters
        if min_count > 1 or max_count is not None:
            counts = [
                (kmer, count) for kmer, count in counts
                if count >= min_count and (max_count is None or count <= max_count)
            ]

        # Sort by count descending
        counts.sort(key=lambda x: x[1], reverse=True)
        return counts

    def save(self, filepath: Union[str, os.PathLike]) -> None:
        """Save index to JSON file.

        This is a wrapper around io.save_index that includes metadata.

        Args:
            filepath: Path to output JSON file.
        """
        # Include metadata in the saved index
        metadata = {
            'k': self.k,
            'subject_ids': self.subject_ids,
            'stats': {
                'total_kmers': self.stats.total_kmers,
                'unique_kmers': self.stats.unique_kmers,
                'total_sequences': self.stats.total_sequences,
                'total_length': self.stats.total_length,
            }
        }

        # Save the actual index
        save_index(self.index, filepath)

        # Save metadata separately (or append to same file)
        meta_file = str(filepath) + '.meta'
        with open(meta_file, 'w') as f:
            json.dump(metadata, f, indent=2)

    @classmethod
    def load(cls, filepath: Union[str, os.PathLike]) -> 'KmerIndex':
        """Load index from JSON file.

        Args:
            filepath: Path to JSON index file.

        Returns:
            KmerIndex object.

        Raises:
            FileNotFoundError: If file doesn't exist.
        """
        # Load the index
        index_data = load_index(filepath)

        # Load metadata
        meta_file = str(filepath) + '.meta'
        if os.path.exists(meta_file):
            with open(meta_file, 'r') as f:
                metadata = json.load(f)
            k = metadata.get('k', 11)
            subject_ids = metadata.get('subject_ids', [])
            stats_data = metadata.get('stats', {})
        else:
            # Fallback: try to determine k from first k-mer
            k = 11
            subject_ids = []
            stats_data = {}

        # Create index object
        index_obj = cls(k)
        index_obj.index = defaultdict(list, index_data)
        index_obj.subject_ids = subject_ids

        # Restore stats
        index_obj.stats = IndexStats(
            total_kmers=stats_data.get('total_kmers', 0),
            unique_kmers=stats_data.get('unique_kmers', len(index_data)),
            total_sequences=stats_data.get('total_sequences', 0),
            total_length=stats_data.get('total_length', 0),
            kmer_size=k,
            memory_estimate_mb=index_obj._estimate_memory()
        )

        return index_obj

    def _estimate_memory(self) -> float:
        """Estimate memory usage in megabytes.

        This is a rough estimate based on the number of entries and
        typical Python object overhead.

        Returns:
            Estimated memory usage in MB.
        """
        # Rough estimate: each k-mer key ~50 bytes, each position tuple ~56 bytes
        num_positions = self.stats.total_kmers
        num_kmers = self.stats.unique_kmers

        # Key overhead: ~50 bytes each
        key_memory = num_kmers * 50
        # Position overhead: ~56 bytes each
        pos_memory = num_positions * 56
        # List overhead: ~56 bytes per k-mer
        list_memory = num_kmers * 56
        # Dict overhead: ~72 bytes per entry
        dict_memory = num_kmers * 72

        total_bytes = key_memory + pos_memory + list_memory + dict_memory
        return total_bytes / (1024 * 1024)

    def prune(self, max_occurrences: int) -> 'KmerIndex':
        """Remove highly repetitive k-mers from the index.

        This is useful for reducing memory usage and improving search
        speed by excluding common repetitive elements.

        Args:
            max_occurrences: Maximum occurrences to keep.

        Returns:
            New KmerIndex with pruned index.

        Examples:
            >>> idx = KmerIndex.build(records, k=3)
            >>> pruned = idx.prune(max_occurrences=10)
        """
        new_idx = KmerIndex(self.k)
        new_idx.subject_ids = self.subject_ids.copy()

        kept_kmers = 0
        for kmer, positions in self.index.items():
            if len(positions) <= max_occurrences:
                new_idx.index[kmer] = positions
                kept_kmers += 1

        # Update statistics
        new_idx.stats = IndexStats(
            total_kmers=sum(len(positions)
                            for positions in new_idx.index.values()),
            unique_kmers=kept_kmers,
            total_sequences=self.stats.total_sequences,
            total_length=self.stats.total_length,
            kmer_size=self.k,
            memory_estimate_mb=new_idx._estimate_memory()
        )

        return new_idx

    def coverage_stats(self) -> Dict[str, float]:
        """Calculate coverage statistics for the indexed sequences.

        Returns:
            Dictionary with coverage metrics:
            - 'mean_coverage': Average number of k-mers covering each position.
            - 'coverage_std': Standard deviation of coverage.
            - 'zero_coverage_positions': Number of positions with no coverage.

        Note:
            This is a simplified estimate and may be slow for large indices.
        """
        # This is a placeholder - full implementation would require
        # loading sequences to calculate true coverage
        return {
            'mean_coverage': self.stats.total_kmers / self.stats.total_length if self.stats.total_length > 0 else 0,
            'coverage_std': 0.0,
            'zero_coverage_positions': 0
        }


# ============================================================================
# Convenience Functions
# ============================================================================

def build_index_from_fasta(
        fasta_file: Union[str, os.PathLike],
        k: int = 11,
        canonical: bool = True,
        min_occurrences: int = 1,
        max_occurrences: Optional[int] = None,
        output_file: Optional[Union[str, os.PathLike]] = None,
        progress: bool = True
) -> KmerIndex:
    """Build a k-mer index from a FASTA file.

    This is a high-level convenience function that handles the complete
    index building workflow from a FASTA file.

    Args:
        fasta_file: Path to FASTA file.
        k: k-mer size.
        canonical: Use canonical k-mers.
        min_occurrences: Minimum occurrences to keep.
        max_occurrences: Maximum occurrences to keep (prune repetitive).
        output_file: If provided, save index to this file.
        progress: Show progress feedback.

    Returns:
        KmerIndex object.

    Examples:
        >>> idx = build_index_from_fasta("database.fa", k=11, output_file="index.json")
    """
    if progress:
        print(
            f"Building index from {fasta_file} with k={k}...", file=sys.stderr)

    idx = KmerIndex.build_from_fasta(
        fasta_file=fasta_file,
        k=k,
        canonical=canonical,
        min_kmer_occurrences=min_occurrences,
        max_kmer_occurrences=max_occurrences,
        progress=progress
    )

    if output_file:
        idx.save(output_file)
        if progress:
            print(f"Index saved to {output_file}", file=sys.stderr)

    return idx


def load_index_from_file(filepath: Union[str, os.PathLike]) -> KmerIndex:
    """Load a k-mer index from a file.

    Wrapper around KmerIndex.load.

    Args:
        filepath: Path to index file.

    Returns:
        KmerIndex object.
    """
    return KmerIndex.load(filepath)


# ============================================================================
# Command Line Interface (for testing)
# ============================================================================

def main():
    """Simple command line interface for testing the index module.

    This allows building and inspecting indices from the command line.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="k-mer index construction tool")
    parser.add_argument("action", choices=["build", "stats", "lookup"],
                        help="Action to perform")
    parser.add_argument("-f", "--fasta", help="Input FASTA file")
    parser.add_argument("-i", "--index", help="Index file")
    parser.add_argument("-k", "--kmer", type=int,
                        default=11, help="k-mer size")
    parser.add_argument("-o", "--output", help="Output file")
    parser.add_argument("--lookup", help="k-mer to lookup")
    parser.add_argument("--canonical", action="store_true", default=True,
                        help="Use canonical k-mers")

    args = parser.parse_args()

    if args.action == "build":
        if not args.fasta:
            print("Error: --fasta required for build action")
            sys.exit(1)

        idx = build_index_from_fasta(
            fasta_file=args.fasta,
            k=args.kmer,
            canonical=args.canonical,
            output_file=args.output
        )
        print(f"Index built: {idx.stats}")

    elif args.action == "stats":
        if not args.index:
            print("Error: --index required for stats action")
            sys.exit(1)

        idx = load_index_from_file(args.index)
        print(f"Index stats: {idx.stats}")
        print(f"Subject IDs: {idx.subject_ids[:5]}...")

    elif args.action == "lookup":
        if not args.index or not args.lookup:
            print("Error: --index and --lookup required for lookup action")
            sys.exit(1)

        idx = load_index_from_file(args.index)
        matches = idx.lookup(args.lookup, canonical=args.canonical)
        print(f"Found {len(matches)} matches for '{args.lookup}':")
        for subject_id, pos in matches[:10]:
            print(f"  {subject_id} at position {pos}")


if __name__ == "__main__":
    main()
