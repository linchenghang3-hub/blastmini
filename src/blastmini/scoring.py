"""Scoring and ranking module for blastmini.

This module provides scoring functions for evaluating alignment quality,
ranking hits by significance, calculating identity percentages, and
formatting results for output.

Key features:
    - Hit scoring and ranking
    - Identity and similarity calculations
    - Score normalization (bitscore calculation)
    - Statistical score estimation (empirical)
    - Hit filtering by score thresholds
    - TSV/CSV output formatting
    - Alignment visualization
"""

import math
import json
from typing import List, Dict, Tuple, Optional, Union, Any, Set
from dataclasses import dataclass, field
from collections import defaultdict
import sys

from .models import SequenceRecord, Hit, AlignmentConfig
from .extension import ExtensionResult

# ============================================================================
# Scoring Matrices
# ============================================================================

# Simple nucleotide scoring matrix (match/mismatch)
NUCLEOTIDE_SCORES = {
    'A': {'A': 5, 'T': -4, 'C': -4, 'G': -4, 'N': -1},
    'T': {'A': -4, 'T': 5, 'C': -4, 'G': -4, 'N': -1},
    'C': {'A': -4, 'T': -4, 'C': 5, 'G': -4, 'N': -1},
    'G': {'A': -4, 'T': -4, 'C': -4, 'G': 5, 'N': -1},
    'N': {'A': -1, 'T': -1, 'C': -1, 'G': -1, 'N': 0},
}

# Alternative scoring matrix (BLASTN default-like)
BLASTN_SCORES = {
    'A': {'A': 1, 'T': -3, 'C': -3, 'G': -3, 'N': -1},
    'T': {'A': -3, 'T': 1, 'C': -3, 'G': -3, 'N': -1},
    'C': {'A': -3, 'T': -3, 'C': 1, 'G': -3, 'N': -1},
    'G': {'A': -3, 'T': -3, 'C': -3, 'G': 1, 'N': -1},
    'N': {'A': -1, 'T': -1, 'C': -1, 'G': -1, 'N': 0},
}

# Match/mismatch scoring (simplest)
MATCH_SCORES = {
    'match': 1,
    'mismatch': -1,
}


def get_nucleotide_score(char1: str, char2: str, matrix: Dict[str, Dict[str, int]] = NUCLEOTIDE_SCORES) -> int:
    """Get score for a pair of nucleotides.

    Args:
        char1: First nucleotide character.
        char2: Second nucleotide character.
        matrix: Scoring matrix to use.

    Returns:
        Score for the pair.

    Examples:
        >>> get_nucleotide_score('A', 'A')
        5
        >>> get_nucleotide_score('A', 'T')
        -4
    """
    char1_upper = char1.upper()
    char2_upper = char2.upper()

    # Handle gaps
    if char1_upper == '-' or char2_upper == '-':
        return 0

    # Default to the matrix
    if char1_upper in matrix and char2_upper in matrix[char1_upper]:
        return matrix[char1_upper][char2_upper]

    # Fallback: use match/mismatch
    if char1_upper == char2_upper:
        return MATCH_SCORES['match']
    else:
        return MATCH_SCORES['mismatch']


def calculate_alignment_score(
        query_alignment: str,
        subject_alignment: str,
        matrix: Dict[str, Dict[str, int]] = NUCLEOTIDE_SCORES
) -> Tuple[int, int, int]:
    """Calculate score, identity, and mismatches from an alignment.

    Args:
        query_alignment: Aligned query sequence (may contain gaps).
        subject_alignment: Aligned subject sequence (may contain gaps).
        matrix: Scoring matrix to use.

    Returns:
        Tuple of (score, identity_count, mismatches_count).

    Examples:
        >>> score, identity, mismatches = calculate_alignment_score("A-T", "A-C")
        >>> print(score, identity, mismatches)
        1 1 1
    """
    if len(query_alignment) != len(subject_alignment):
        raise ValueError("Alignment sequences must have the same length")

    score = 0
    identity = 0
    mismatches = 0

    for q_char, s_char in zip(query_alignment, subject_alignment):
        if q_char == '-' or s_char == '-':
            # Gaps don't contribute to score in ungapped extension
            continue

        if q_char.upper() == s_char.upper():
            identity += 1
            score += get_nucleotide_score(q_char, s_char, matrix)
        else:
            mismatches += 1
            score += get_nucleotide_score(q_char, s_char, matrix)

    return score, identity, mismatches


# ============================================================================
# Hit Scoring and Ranking
# ============================================================================

@dataclass
class ScoredHit:
    """A hit with additional scoring information.

    This extends the basic Hit with scoring metrics for ranking.

    Attributes:
        hit: The original Hit object.
        raw_score: Raw alignment score.
        bit_score: Bit score (normalized).
        evalue: E-value (statistical significance).
        identity_percent: Percentage of identical characters.
        similarity_percent: Percentage of similar characters.
        coverage_percent: Percentage of query covered.
        rank: Ranking position among hits.
        normalized_score: Score normalized by alignment length.
    """
    hit: Hit
    raw_score: int = 0
    bit_score: float = 0.0
    evalue: float = 1.0
    identity_percent: float = 0.0
    similarity_percent: float = 0.0
    coverage_percent: float = 0.0
    rank: int = 0
    normalized_score: float = 0.0

    def __repr__(self) -> str:
        return (f"ScoredHit(rank={self.rank}, subject='{self.hit.subject_id}', "
                f"score={self.raw_score}, bitscore={self.bit_score:.1f}, "
                f"evalue={self.evalue:.2e}, identity={self.identity_percent:.1f}%)")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON output."""
        return {
            'rank': self.rank,
            'query_id': self.hit.query_id,
            'subject_id': self.hit.subject_id,
            'raw_score': self.raw_score,
            'bit_score': self.bit_score,
            'evalue': self.evalue,
            'identity_percent': self.identity_percent,
            'similarity_percent': self.similarity_percent,
            'coverage_percent': self.coverage_percent,
            'alignment_length': self.hit.alignment_length,
            'query_start': self.hit.query_start,
            'query_end': self.hit.query_end,
            'subject_start': self.hit.subject_start,
            'subject_end': self.hit.subject_end,
            'query_alignment': self.hit.query_alignment,
            'subject_alignment': self.hit.subject_alignment,
        }


class HitScorer:
    """Score and rank hits from extension results.

    This class handles the scoring, ranking, and filtering of hits,
    including bit score calculation and E-value estimation.

    Attributes:
        config: Alignment configuration.
        matrix: Scoring matrix to use.
        min_score: Minimum score threshold.
        min_identity: Minimum identity threshold.
        top_n: Number of top hits to keep.
    """

    def __init__(
            self,
            config: Optional[AlignmentConfig] = None,
            matrix: Dict[str, Dict[str, int]] = NUCLEOTIDE_SCORES,
            min_score: int = 0,
            min_identity: float = 0.0,
            top_n: int = 10
    ):
        """Initialize the hit scorer.

        Args:
            config: Alignment configuration.
            matrix: Scoring matrix to use.
            min_score: Minimum score to keep.
            min_identity: Minimum identity percentage.
            top_n: Number of top hits to keep.
        """
        self.config = config or AlignmentConfig()
        self.matrix = matrix
        self.min_score = min_score
        self.min_identity = min_identity
        self.top_n = top_n

    def score_hit(self, hit: Hit, query_length: int) -> ScoredHit:
        """Calculate all scoring metrics for a single hit.

        Args:
            hit: Hit object to score.
            query_length: Length of the query sequence.

        Returns:
            ScoredHit with all metrics calculated.
        """
        # Use raw score if available
        raw_score = hit.score

        # Calculate identity percentage
        if hit.alignment_length > 0:
            identity_percent = (hit.identity_percent / 100)  # Already a percentage
            # Recalculate if needed
            if hit.query_alignment and hit.subject_alignment:
                _, identity_count, _ = calculate_alignment_score(
                    hit.query_alignment, hit.subject_alignment, self.matrix
                )
                identity_percent = identity_count / hit.alignment_length * 100
            else:
                identity_percent = hit.identity_percent
        else:
            identity_percent = 0.0

        # Calculate coverage
        if query_length > 0:
            coverage_percent = (hit.alignment_length / query_length) * 100
        else:
            coverage_percent = 0.0

        # Calculate bit score (simplified)
        bit_score = self._calculate_bit_score(raw_score)

        # Calculate normalized score (per base)
        normalized_score = raw_score / hit.alignment_length if hit.alignment_length > 0 else 0

        # Calculate similarity (simplified - same as identity for ungapped)
        similarity_percent = identity_percent

        # Create ScoredHit
        scored = ScoredHit(
            hit=hit,
            raw_score=raw_score,
            bit_score=bit_score,
            evalue=1.0,  # Will be updated separately
            identity_percent=identity_percent,
            similarity_percent=similarity_percent,
            coverage_percent=coverage_percent,
            normalized_score=normalized_score
        )

        return scored

    def score_hits(
            self,
            hits: List[Hit],
            query_length: int,
            db_total_length: Optional[int] = None,
            progress: bool = True
    ) -> List[ScoredHit]:
        """Score and rank a list of hits.

        Args:
            hits: List of Hit objects.
            query_length: Length of the query sequence.
            db_total_length: Total length of the database (for E-value).
            progress: Show progress feedback.

        Returns:
            List of ScoredHit objects sorted by score.
        """
        if not hits:
            return []

        if progress and len(hits) > 100:
            print(f"Scoring {len(hits)} hits...", file=sys.stderr)

        # Score each hit
        scored_hits = []
        for i, hit in enumerate(hits):
            scored = self.score_hit(hit, query_length)
            scored_hits.append(scored)

            if progress and len(hits) > 100 and (i + 1) % 50 == 0:
                print(f"  Scored {i + 1}/{len(hits)} hits", file=sys.stderr)

        # Calculate E-values
        if db_total_length is not None and query_length > 0:
            self._estimate_evalues(scored_hits, query_length, db_total_length)

        # Filter by thresholds
        filtered = self._filter_hits(scored_hits)

        # Sort and rank
        filtered.sort(key=lambda s: (s.raw_score, s.identity_percent), reverse=True)

        # Assign ranks
        for i, scored in enumerate(filtered):
            scored.rank = i + 1

        # Keep top N
        if self.top_n and len(filtered) > self.top_n:
            filtered = filtered[:self.top_n]

        if progress and len(hits) > 100:
            print(f"  Done: kept {len(filtered)} hits", file=sys.stderr)

        return filtered

    def _calculate_bit_score(self, raw_score: int) -> float:
        """Calculate bit score from raw score.

        Bit score is a normalized score that allows comparison between
        different scoring systems.

        Args:
            raw_score: Raw alignment score.

        Returns:
            Bit score (simplified calculation).
        """
        # Simplified bit score calculation
        # In real BLAST: bit_score = (raw_score * λ) / ln(2)
        # where λ is a scaling factor depending on the scoring system
        # We use a simple approximation

        # For our scoring system (match=5, mismatch=-4)
        # The expected score per base for random sequences is ~0.5
        # λ is approximately 0.1 for this system
        lambda_value = 0.1

        if raw_score > 0:
            bit_score = (raw_score * lambda_value) / math.log(2)
        else:
            bit_score = 0.0

        return bit_score

    def _estimate_evalues(
            self,
            scored_hits: List[ScoredHit],
            query_length: int,
            db_total_length: int
    ) -> None:
        """Estimate E-values for a list of hits.

        This is a simplified E-value estimation based on score distribution.

        Args:
            scored_hits: List of ScoredHit objects (modified in place).
            query_length: Length of the query sequence.
            db_total_length: Total length of the database.
        """
        if not scored_hits:
            return

        # Extract raw scores
        raw_scores = [s.raw_score for s in scored_hits]
        max_score = max(raw_scores) if raw_scores else 0
        min_score = min(raw_scores) if raw_scores else 0
        score_range = max_score - min_score if max_score > min_score else 1

        # Estimate expected score for random alignment
        # Using the scoring system properties
        expected_score = self._estimate_expected_score()

        # Estimate E-value for each hit
        for scored in scored_hits:
            # Simple heuristic: E-value decreases exponentially with score
            score_diff = scored.raw_score - expected_score

            if score_diff > 0:
                # More significant = smaller E-value
                # Use exponential decay based on score difference
                evalue = math.exp(-score_diff * 0.5)

                # Adjust for database size and query length
                # More database sequences = more chance matches
                effective_length = db_total_length - (query_length * len(scored_hits))
                if effective_length < 0:
                    effective_length = db_total_length

                evalue = evalue * (effective_length / 1000) * 0.01

                # Cap E-value at 1.0
                scored.evalue = min(evalue, 1.0)
            else:
                scored.evalue = 1.0

    def _estimate_expected_score(self) -> float:
        """Estimate the expected score for random sequences.

        Returns:
            Expected score per alignment.
        """
        # Calculate expected score per base from the scoring matrix
        # For nucleotides, assume uniform distribution
        base_probs = {'A': 0.25, 'T': 0.25, 'C': 0.25, 'G': 0.25}

        expected_per_base = 0
        for base1, prob1 in base_probs.items():
            for base2, prob2 in base_probs.items():
                score = get_nucleotide_score(base1, base2, self.matrix)
                expected_per_base += prob1 * prob2 * score

        # For an alignment length of ~100 bases
        expected_score = expected_per_base * 100

        return expected_score

    def _filter_hits(self, scored_hits: List[ScoredHit]) -> List[ScoredHit]:
        """Filter hits by thresholds.

        Args:
            scored_hits: List of ScoredHit objects.

        Returns:
            Filtered list.
        """
        filtered = []
        for scored in scored_hits:
            # Check score threshold
            if scored.raw_score < self.min_score:
                continue

            # Check identity threshold
            if scored.identity_percent < self.min_identity:
                continue

            # Check if alignment is valid
            if scored.hit.alignment_length < 5:  # Minimum alignment length
                continue

            filtered.append(scored)

        return filtered

    def merge_hits(
            self,
            hits: List[Hit],
            subject_sequences: Optional[Dict[str, str]] = None
    ) -> List[Hit]:
        """Merge overlapping hits on the same subject.

        Args:
            hits: List of Hit objects.
            subject_sequences: Dictionary mapping subject ID to sequence.

        Returns:
            List of merged Hit objects.
        """
        if not hits:
            return []

        # Group hits by subject_id
        by_subject: Dict[str, List[Hit]] = defaultdict(list)
        for hit in hits:
            by_subject[hit.subject_id].append(hit)

        merged = []

        for subject_id, subject_hits in by_subject.items():
            # Sort by query_start
            subject_hits.sort(key=lambda h: h.query_start)

            # Merge overlapping hits
            current = subject_hits[0]
            for hit in subject_hits[1:]:
                # Check if hits overlap or are close on query
                if hit.query_start <= current.query_end + 10:
                    # Merge: extend the current hit
                    current.query_end = max(current.query_end, hit.query_end)
                    current.subject_end = max(current.subject_end, hit.subject_end)
                    current.score = max(current.score, hit.score)
                    current.alignment_length = current.query_end - current.query_start

                    # Combine alignment sequences (simplified)
                    if hit.query_alignment:
                        current.query_alignment += hit.query_alignment
                    if hit.subject_alignment:
                        current.subject_alignment += hit.subject_alignment
                else:
                    merged.append(current)
                    current = hit

            merged.append(current)

        return merged


# ============================================================================
# Output Formatting
# ============================================================================

def format_hit_as_text(scored: ScoredHit, show_alignment: bool = True) -> str:
    """Format a scored hit as human-readable text.

    Args:
        scored: ScoredHit object.
        show_alignment: Whether to show the alignment.

    Returns:
        Formatted string.

    Examples:
        >>> print(format_hit_as_text(scored))
        Rank: 1
        Subject: S1
        Score: 145 (bitscore: 50.0)
        E-value: 1.23e-05
        Identity: 98.5%
        Alignment: 200/200
    """
    lines = []
    lines.append(f"Rank: {scored.rank}")
    lines.append(f"Subject: {scored.hit.subject_id}")
    lines.append(f"Score: {scored.raw_score} (bitscore: {scored.bit_score:.1f})")
    lines.append(f"E-value: {scored.evalue:.2e}")
    lines.append(f"Identity: {scored.identity_percent:.1f}%")
    lines.append(f"Coverage: {scored.coverage_percent:.1f}%")
    lines.append(f"Alignment length: {scored.hit.alignment_length}")
    lines.append(f"Query region: {scored.hit.query_start}-{scored.hit.query_end}")
    lines.append(f"Subject region: {scored.hit.subject_start}-{scored.hit.subject_end}")

    if show_alignment and scored.hit.query_alignment:
        lines.append("")
        lines.append("Query:   " + scored.hit.query_alignment)
        lines.append("         " + "|" * len(scored.hit.query_alignment))
        lines.append("Subject: " + scored.hit.subject_alignment)

    return "\n".join(lines)


def format_hits_as_text(scored_hits: List[ScoredHit], top_n: int = 10) -> str:
    """Format multiple hits as human-readable text.

    Args:
        scored_hits: List of ScoredHit objects.
        top_n: Number of top hits to show.

    Returns:
        Formatted string.
    """
    lines = []
    lines.append("=" * 80)
    lines.append(f"BLAST Search Results ({len(scored_hits)} hits found)")
    lines.append("=" * 80)
    lines.append("")

    for scored in scored_hits[:top_n]:
        lines.append(format_hit_as_text(scored, show_alignment=True))
        lines.append("-" * 80)
        lines.append("")

    if len(scored_hits) > top_n:
        lines.append(f"... and {len(scored_hits) - top_n} more hits")

    return "\n".join(lines)


def format_hits_as_tsv(scored_hits: List[ScoredHit]) -> str:
    """Format hits as TSV with header matching Hit.tsv_header.
    
    The column order is: query_id, subject_id, score, identity_percent,
    alignment_length, query_start, query_end, subject_start, subject_end,
    evalue, bit_score.
    """
    header = "\t".join([
        "query_id", "subject_id", "score", "identity_percent",
        "alignment_length", "query_start", "query_end",
        "subject_start", "subject_end", "evalue", "bit_score"
    ])
    lines = [header]
    for scored in scored_hits:
        line = "\t".join([
            scored.hit.query_id,
            scored.hit.subject_id,
            str(scored.raw_score),
            f"{scored.identity_percent:.2f}",
            str(scored.hit.alignment_length),
            str(scored.hit.query_start),
            str(scored.hit.query_end),
            str(scored.hit.subject_start),
            str(scored.hit.subject_end),
            f"{scored.evalue:.2e}" if scored.evalue is not None else "",
            f"{scored.bit_score:.1f}" if scored.bit_score is not None else "",
        ])
        lines.append(line)
    return "\n".join(lines)


def format_hits_as_json(scored_hits: List[ScoredHit]) -> str:
    """Format hits as JSON.

    Args:
        scored_hits: List of ScoredHit objects.

    Returns:
        JSON formatted string.
    """
    data = {
        "total_hits": len(scored_hits),
        "hits": [scored.to_dict() for scored in scored_hits]
    }
    return json.dumps(data, indent=2)


def format_hits_as_bed(scored_hits: List[ScoredHit]) -> str:
    """Format hits as BED (Browser Extensible Data) format.

    Args:
        scored_hits: List of ScoredHit objects.

    Returns:
        BED formatted string.
    """
    lines = [
        "#track name=blastmini type=bedDetail description=\"BLAST hits\"",
        "#chrom\tsrc\tstart\tend\tname\tscore\tstrand"
    ]

    for scored in scored_hits:
        line = "\t".join([
            scored.hit.subject_id,
            str(scored.hit.subject_start),
            str(scored.hit.subject_end),
            f"{scored.hit.query_id}:{scored.rank}",
            str(scored.raw_score),
            "+"
        ])
        lines.append(line)

    return "\n".join(lines)


# ============================================================================
# Alignment Visualization
# ============================================================================

def visualize_alignment(
        query_alignment: str,
        subject_alignment: str,
        line_width: int = 60,
        show_consensus: bool = True
) -> str:
    """Create a formatted alignment visualization.

    Args:
        query_alignment: Aligned query sequence.
        subject_alignment: Aligned subject sequence.
        line_width: Number of characters per line.
        show_consensus: Show consensus line (| for match, . for mismatch).

    Returns:
        Formatted alignment string.

    Examples:
        >>> print(visualize_alignment("ATCGATCG", "ATCGATCG"))
        Query:   ATCGATCG
                 ||||||||
        Subject: ATCGATCG
    """
    if len(query_alignment) != len(subject_alignment):
        raise ValueError("Alignment sequences must have the same length")

    lines = []
    total_len = len(query_alignment)

    for start in range(0, total_len, line_width):
        end = min(start + line_width, total_len)

        q_seg = query_alignment[start:end]
        s_seg = subject_alignment[start:end]

        # Add sequence lines
        lines.append(f"Query:   {q_seg}")

        # Add consensus line
        if show_consensus:
            consensus = []
            for q, s in zip(q_seg, s_seg):
                if q == '-' or s == '-':
                    consensus.append(' ')
                elif q.upper() == s.upper():
                    consensus.append('|')
                else:
                    consensus.append('.')
            lines.append(f"         {''.join(consensus)}")

        lines.append(f"Subject: {s_seg}")
        lines.append("")

    return "\n".join(lines)


# ============================================================================
# Command Line Interface (for testing)
# ============================================================================

def main():
    """Simple command line interface for testing scoring."""
    import argparse

    parser = argparse.ArgumentParser(description="Hit scoring tool")
    parser.add_argument("-i", "--input", required=True,
                        help="Input hits file (TSV)")
    parser.add_argument("-o", "--output", help="Output file")
    parser.add_argument("--format", choices=["text", "tsv", "json", "bed"], default="text",
                        help="Output format")
    parser.add_argument("--top", type=int, default=10,
                        help="Number of top hits")

    args = parser.parse_args()

    # Load hits
    from .io import load_hits_from_tsv

    hits = load_hits_from_tsv(args.input)
    print(f"Loaded {len(hits)} hits", file=sys.stderr)

    # Create scorer
    scorer = HitScorer(top_n=args.top)

    # Score hits (assuming query length from data)
    query_length = 1000  # Placeholder
    scored_hits = scorer.score_hits(hits, query_length)

    # Format output
    if args.format == "text":
        output = format_hits_as_text(scored_hits, top_n=args.top)
    elif args.format == "tsv":
        output = format_hits_as_tsv(scored_hits)
    elif args.format == "json":
        output = format_hits_as_json(scored_hits)
    elif args.format == "bed":
        output = format_hits_as_bed(scored_hits)
    else:
        output = format_hits_as_text(scored_hits, top_n=args.top)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(output)
        print(f"Output written to {args.output}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()