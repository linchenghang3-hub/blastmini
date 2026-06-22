"""Data models for blastmini.

This module defines the core data structures used throughout the blastmini
package, including sequence records, alignment hits, and configuration
parameters.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class SequenceRecord:
    """A biological sequence record with identifier and sequence data.

    This class represents a single sequence entry from a FASTA file or
    other sequence data source. It is designed to be immutable and
    hashable for use in indices and caches.

    Attributes:
        id: Unique sequence identifier (without the '>' character).
        sequence: The biological sequence string (DNA/RNA/protein).
        description: Optional description line from the FASTA header.

    Examples:
        >>> record = SequenceRecord("NC_000001", "ATCGATCGATCG")
        >>> print(record.id)
        NC_000001
    """

    id: str
    sequence: str
    description: Optional[str] = ""

    def __post_init__(self) -> None:
        """Normalize sequence to uppercase and strip whitespace."""
        self.sequence = self.sequence.strip().upper()
        if self.description is None:
            self.description = ""

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self.sequence)

    def __repr__(self) -> str:
        """Compact string representation."""
        seq_preview = self.sequence[:30] + "..." if len(self.sequence) > 30 else self.sequence
        return f"SequenceRecord(id='{self.id}', length={len(self)}, sequence='{seq_preview}')"


@dataclass
class Hit:
    """A single alignment hit between a query and a subject sequence.

    This class stores all information about a seed-and-extend alignment
    result, including coordinates, scores, and alignment details.

    Attributes:
        query_id: Identifier of the query sequence.
        subject_id: Identifier of the subject (database) sequence.
        score: Raw alignment score from the extension process.
        identity_percent: Percentage of identical residues in the alignment.
        alignment_length: Total length of the alignment region.
        query_start: Start position in the query sequence (0-based).
        query_end: End position in the query sequence (0-based, exclusive).
        subject_start: Start position in the subject sequence (0-based).
        subject_end: End position in the subject sequence (0-based, exclusive).
        query_alignment: The query sequence segment with gaps (if any).
        subject_alignment: The subject sequence segment with gaps (if any).
        evalue: Estimated E-value for this hit (if computed).
        bit_score: Bit score for this hit (if computed).

    Examples:
        >>> hit = Hit(query_id="Q1", subject_id="S1", score=45, identity_percent=95.0)
        >>> print(hit)
        Hit(query='Q1', subject='S1', score=45, identity=95.0%)
    """

    query_id: str
    subject_id: str
    score: int = 0
    identity_percent: float = 0.0
    alignment_length: int = 0
    query_start: int = -1
    query_end: int = -1
    subject_start: int = -1
    subject_end: int = -1
    query_alignment: str = ""
    subject_alignment: str = ""
    evalue: Optional[float] = None
    bit_score: Optional[float] = None

    def __post_init__(self) -> None:
        """Validate that coordinates are consistent."""
        if self.query_start >= 0 and self.query_end >= 0:
            if self.query_end <= self.query_start:
                raise ValueError(
                    f"query_end ({self.query_end}) must be > query_start ({self.query_start})"
                )
        if self.subject_start >= 0 and self.subject_end >= 0:
            if self.subject_end <= self.subject_start:
                raise ValueError(
                    f"subject_end ({self.subject_end}) must be > subject_start ({self.subject_start})"
                )

    def __repr__(self) -> str:
        """Compact string representation."""
        return (f"Hit(query='{self.query_id}', subject='{self.subject_id}', "
                f"score={self.score}, identity={self.identity_percent:.1f}%)")

    def to_tsv(self) -> str:
        """Convert hit to tab-separated values for output.

        Returns:
            A tab-separated string suitable for TSV output files.

        Examples:
            >>> hit = Hit(query_id="Q1", subject_id="S1", score=100)
            >>> print(hit.to_tsv())
            Q1\tS1\t100\t...
        """
        return "\t".join([
            self.query_id,
            self.subject_id,
            str(self.score),
            f"{self.identity_percent:.2f}",
            str(self.alignment_length),
            str(self.query_start),
            str(self.query_end),
            str(self.subject_start),
            str(self.subject_end),
            str(self.evalue if self.evalue is not None else ""),
            str(self.bit_score if self.bit_score is not None else ""),
        ])

    @classmethod
    def tsv_header(cls) -> str:
        """Return the TSV header line for hit output.

        Returns:
            A tab-separated header string.
        """
        return "\t".join([
            "query_id", "subject_id", "score", "identity_percent",
            "alignment_length", "query_start", "query_end",
            "subject_start", "subject_end", "evalue", "bit_score"
        ])


@dataclass
class AlignmentConfig:
    """Configuration parameters for alignment algorithms.

    This class centralizes all scoring and alignment parameters to avoid
    passing many individual arguments to functions.

    Attributes:
        kmer_size: Size of k-mers used for seeding (default: 11).
        match_score: Score for matching nucleotides (default: 5).
        mismatch_penalty: Score for mismatching nucleotides (negative, default: -4).
        gap_open_penalty: Penalty for opening a gap (negative, default: -5).
        gap_extend_penalty: Penalty for extending a gap (negative, default: -2).
        x_dropoff: X-dropoff threshold for extension termination (default: 5).
        top_n: Number of top hits to report (default: 10).
        min_score: Minimum score to report a hit (default: 0).
        min_identity: Minimum identity percentage to report (default: 0.0).

    Examples:
        >>> config = AlignmentConfig(kmer_size=15, match_score=1, mismatch_penalty=-1)
        >>> print(config.kmer_size)
        15
    """

    kmer_size: int = 11
    match_score: int = 5
    mismatch_penalty: int = -4
    gap_open_penalty: int = -5
    gap_extend_penalty: int = -2
    x_dropoff: int = 5
    top_n: int = 10
    min_score: int = 0
    min_identity: float = 0.0

    def __post_init__(self) -> None:
        """Validate configuration parameters."""
        if self.kmer_size < 1:
            raise ValueError(f"kmer_size must be >= 1, got {self.kmer_size}")
        if self.match_score < 1:
            raise ValueError(f"match_score must be positive, got {self.match_score}")
        if self.mismatch_penalty > 0:
            raise ValueError(f"mismatch_penalty must be <= 0, got {self.mismatch_penalty}")
        if self.x_dropoff < 1:
            raise ValueError(f"x_dropoff must be >= 1, got {self.x_dropoff}")
        if self.top_n < 1:
            raise ValueError(f"top_n must be >= 1, got {self.top_n}")