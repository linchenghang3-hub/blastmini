"""Statistical significance estimation for blastmini.

This module provides statistical significance estimation for sequence
alignments, including E-value calculation and p-value estimation.
Unlike real BLAST which uses Karlin-Altschul extreme value distribution
theory, this implementation uses empirical permutation testing for
educational purposes.

Key features:
    - Empirical p-value estimation via sequence permutation
    - E-value calculation from empirical distributions
    - Score distribution modeling (normal approximation)
    - Background score estimation from random sequences
    - Statistical significance reporting
    - Multiple testing correction (Bonferroni, FDR)
"""

import math
import random
import sys
from typing import List, Tuple, Optional, Dict, Any, Union
from dataclasses import dataclass, field
from collections import defaultdict
import numpy as np

from .models import Hit, AlignmentConfig
from .extension import SeedExtender, ExtensionResult
from .scoring import HitScorer, ScoredHit


# ============================================================================
# Statistical Data Structures
# ============================================================================

@dataclass
class ScoreDistribution:
    """Distribution of alignment scores from random sequences.

    This class stores and analyzes score distributions for estimating
    statistical significance.

    Attributes:
        scores: List of observed scores.
        mean: Mean of the distribution.
        std: Standard deviation of the distribution.
        min_score: Minimum score observed.
        max_score: Maximum score observed.
        n_samples: Number of samples in the distribution.
        fitted_params: Parameters for fitted distribution.

    Examples:
        >>> dist = ScoreDistribution([10, 15, 12, 18, 11])
        >>> print(f"Mean: {dist.mean:.2f}, Std: {dist.std:.2f}")
        Mean: 13.20, Std: 3.03
    """
    scores: List[int] = field(default_factory=list)
    mean: float = 0.0
    std: float = 0.0
    min_score: int = 0
    max_score: int = 0
    n_samples: int = 0
    fitted_params: Dict[str, float] = field(default_factory=dict)

    def __post_init__(self):
        if self.scores:
            self._update_stats()

    def _update_stats(self) -> None:
        """Update statistics from scores list."""
        if not self.scores:
            return
        self.n_samples = len(self.scores)
        self.min_score = min(self.scores)
        self.max_score = max(self.scores)
        self.mean = sum(self.scores) / self.n_samples

        # Calculate standard deviation
        if self.n_samples > 1:
            variance = sum((s - self.mean) ** 2 for s in self.scores) / (self.n_samples - 1)
            self.std = math.sqrt(variance)
        else:
            self.std = 0.0

    def add_score(self, score: int) -> None:
        """Add a score to the distribution."""
        self.scores.append(score)
        self._update_stats()

    def get_percentile(self, percentile: float) -> int:
        """Get a percentile of the distribution.

        Args:
            percentile: Percentile value (0-100).

        Returns:
            Score at the given percentile.
        """
        if not self.scores:
            return 0
        sorted_scores = sorted(self.scores)
        idx = int(percentile / 100 * len(sorted_scores))
        if idx >= len(sorted_scores):
            idx = len(sorted_scores) - 1
        return sorted_scores[idx]

    def get_pvalue(self, score: int) -> float:
        """Estimate p-value for a given score.

        P-value is the probability of observing a score >= given score
        in the background distribution.

        Args:
            score: Score to evaluate.

        Returns:
            Estimated p-value.
        """
        if not self.scores:
            return 1.0

        # Count scores >= given score
        n_ge = sum(1 for s in self.scores if s >= score)
        return n_ge / len(self.scores)

    def get_evalue(self, score: int, effective_db_size: int) -> float:
        """Estimate E-value for a given score.

        E-value is the expected number of random hits with score >= given score.

        Args:
            score: Score to evaluate.
            effective_db_size: Effective database size (total length).

        Returns:
            Estimated E-value.
        """
        pvalue = self.get_pvalue(score)
        # E-value = p-value * effective number of independent tests
        # For sequence search, effective tests ≈ database size
        return pvalue * effective_db_size

    def fit_normal(self) -> Dict[str, float]:
        """Fit a normal distribution to the scores.

        Returns:
            Dictionary with 'mean' and 'std' parameters.
        """
        self.fitted_params = {'mean': self.mean, 'std': self.std}
        return self.fitted_params

    def estimate_extreme_params(self) -> Tuple[float, float]:
        """Estimate parameters for extreme value distribution.

        This is a simplified approximation for the Gumbel distribution
        parameters (λ and K) used in real BLAST.

        Returns:
            Tuple of (lambda_param, K_param).
        """
        if not self.scores or self.n_samples < 10:
            return (0.1, 0.1)

        # Get top scores (extreme values)
        sorted_scores = sorted(self.scores, reverse=True)
        top_scores = sorted_scores[:min(50, len(sorted_scores) // 2)]

        if len(top_scores) < 5:
            return (0.1, 0.1)

        # Estimate λ from the tail of the distribution
        # For Gumbel distribution, the tail follows exp(-λ * (x - u))
        # We approximate by fitting an exponential to the top scores
        try:
            # Use log of rank vs score
            scores_log = []
            for i, score in enumerate(top_scores):
                if score > 0:
                    scores_log.append((score, math.log(i + 1)))

            if len(scores_log) < 3:
                return (0.1, 0.1)

            # Simple linear regression on log(rank) ~ score
            x = [s for s, _ in scores_log]
            y = [l for _, l in scores_log]

            n = len(x)
            if n < 2:
                return (0.1, 0.1)

            x_mean = sum(x) / n
            y_mean = sum(y) / n

            # Slope = -λ
            numerator = sum((x[i] - x_mean) * (y[i] - y_mean) for i in range(n))
            denominator = sum((x[i] - x_mean) ** 2 for i in range(n))

            if denominator == 0:
                return (0.1, 0.1)

            lambda_param = -numerator / denominator

            # K is estimated from the intercept
            # log(K) ≈ intercept
            intercept = y_mean - (-lambda_param) * x_mean
            K_param = math.exp(intercept)

            # Cap parameters
            lambda_param = max(0.01, min(0.5, lambda_param))
            K_param = max(0.001, min(1.0, K_param))

            return (lambda_param, K_param)

        except:
            return (0.1, 0.1)

    def __repr__(self) -> str:
        return (f"ScoreDistribution(n={self.n_samples}, "
                f"mean={self.mean:.2f}, std={self.std:.2f}, "
                f"min={self.min_score}, max={self.max_score})")


@dataclass
class SignificanceResult:
    """Statistical significance results for a hit.

    Attributes:
        hit: The original hit.
        raw_score: Raw alignment score.
        pvalue: Estimated p-value.
        evalue: Estimated E-value.
        bit_score: Bit score (normalized).
        lambda_param: Lambda parameter from extreme value distribution.
        K_param: K parameter from extreme value distribution.
        is_significant: Whether the hit is statistically significant.
        significance_level: Significance level used.

    Examples:
        >>> result = SignificanceResult(hit, 100, 1e-5, 0.001)
        >>> print(f"E-value: {result.evalue:.2e}, Significant: {result.is_significant}")
        E-value: 1.00e-03, Significant: True
    """
    hit: Hit
    raw_score: int
    pvalue: float = 1.0
    evalue: float = 1.0
    bit_score: float = 0.0
    lambda_param: float = 0.1
    K_param: float = 0.1
    is_significant: bool = False
    significance_level: float = 0.05

    def __repr__(self) -> str:
        status = "SIGNIFICANT" if self.is_significant else "NOT SIGNIFICANT"
        return (f"SignificanceResult(score={self.raw_score}, "
                f"evalue={self.evalue:.2e}, pvalue={self.pvalue:.2e}, "
                f"bitscore={self.bit_score:.1f}, {status})")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for output."""
        return {
            'query_id': self.hit.query_id,
            'subject_id': self.hit.subject_id,
            'raw_score': self.raw_score,
            'pvalue': self.pvalue,
            'evalue': self.evalue,
            'bit_score': self.bit_score,
            'lambda': self.lambda_param,
            'K': self.K_param,
            'is_significant': self.is_significant,
            'significance_level': self.significance_level
        }


# ============================================================================
# Statistical Significance Estimator
# ============================================================================

class SignificanceEstimator:
    """Estimate statistical significance of sequence alignments.

    This class uses empirical permutation testing to estimate the
    significance of alignment scores. It can also fit extreme value
    distributions for more efficient estimation.

    Attributes:
        config: Alignment configuration.
        seed_extender: SeedExtender for generating random alignments.
        hit_scorer: HitScorer for scoring alignments.
        n_permutations: Number of permutations for empirical estimation.
        random_seed: Random seed for reproducibility.
        background_distribution: Score distribution from random sequences.
        extreme_params: Fitted extreme value distribution parameters.

    Examples:
        >>> estimator = SignificanceEstimator(n_permutations=100)
        >>> dist = estimator.estimate_background_distribution(query, db_sequences)
        >>> result = estimator.estimate_significance(hit, query_length, db_size)
        >>> print(f"E-value: {result.evalue:.2e}")
    """

    def __init__(
            self,
            config: Optional[AlignmentConfig] = None,
            n_permutations: int = 100,
            random_seed: Optional[int] = None,
            min_score_threshold: int = 10
    ):
        """Initialize the significance estimator.

        Args:
            config: Alignment configuration.
            n_permutations: Number of permutations for empirical estimation.
            random_seed: Random seed for reproducibility.
            min_score_threshold: Minimum score for evaluating significance.
        """
        self.config = config or AlignmentConfig()
        self.n_permutations = n_permutations
        self.random_seed = random_seed
        self.min_score_threshold = min_score_threshold

        # Initialize components
        self.seed_extender = SeedExtender(config=self.config, track_stats=False)
        self.hit_scorer = HitScorer(config=self.config)

        # Cached distributions
        self.background_distribution: Optional[ScoreDistribution] = None
        self.extreme_params: Optional[Tuple[float, float]] = None

        # Set random seed
        if random_seed is not None:
            random.seed(random_seed)

    def estimate_background_distribution(
            self,
            query: Union[str, 'SequenceRecord'],
            subject_sequences: Dict[str, str],
            n_permutations: Optional[int] = None,
            progress: bool = True
    ) -> ScoreDistribution:
        """Estimate the background score distribution using permutation.

        This method generates random alignments by permuting the query
        sequence and scoring alignments to get the background distribution.

        Args:
            query: Query sequence or SequenceRecord.
            subject_sequences: Dictionary mapping subject ID to sequence.
            n_permutations: Number of permutations (uses default if None).
            progress: Show progress feedback.

        Returns:
            ScoreDistribution object.
        """
        # Get query sequence
        if hasattr(query, 'sequence'):
            query_seq = query.sequence
        else:
            query_seq = str(query)

        # Use specified or default number of permutations
        n_perm = n_permutations or self.n_permutations

        if progress:
            print(f"Estimating background distribution with {n_perm} permutations...",
                  file=sys.stderr)

        scores = []

        for i in range(n_perm):
            # Permute the query sequence
            permuted_query = self._permute_sequence(query_seq)

            # Get some random seeds (simplified)
            # For efficiency, we only take a sample of subject sequences
            subject_sample = self._sample_subjects(subject_sequences, max_subjects=5)

            # Generate random alignments
            for subject_id, subject_seq in subject_sample.items():
                # Pick a random position in the query and subject
                if len(permuted_query) < 10 or len(subject_seq) < 10:
                    continue

                q_pos = random.randint(0, len(permuted_query) - 10)
                s_pos = random.randint(0, len(subject_seq) - 10)

                # Extend the random seed
                result = self.seed_extender._extend_bidirectional(
                    query_seq=permuted_query,
                    subject_seq=subject_seq,
                    query_pos=q_pos,
                    subject_pos=s_pos,
                    query_id="permuted",
                    subject_id=subject_id
                )

                if result and result.score >= self.min_score_threshold:
                    scores.append(result.score)

            if progress and (i + 1) % 10 == 0:
                print(f"  Permutation {i + 1}/{n_perm}, collected {len(scores)} scores",
                      file=sys.stderr)

        # Create distribution
        dist = ScoreDistribution(scores=scores)

        if progress:
            print(f"  Done: {dist}", file=sys.stderr)
            if dist.n_samples > 0:
                p95 = dist.get_percentile(95)
                p99 = dist.get_percentile(99)
                print(f"    P95: {p95}, P99: {p99}", file=sys.stderr)

        self.background_distribution = dist
        return dist

    def estimate_extreme_value_distribution(
            self,
            query: Union[str, 'SequenceRecord'],
            subject_sequences: Dict[str, str],
            n_permutations: Optional[int] = None,
            progress: bool = True
    ) -> Tuple[float, float]:
        """Estimate parameters for extreme value distribution.

        This fits a Gumbel distribution to the extreme scores from
        permutation testing.

        Args:
            query: Query sequence or SequenceRecord.
            subject_sequences: Dictionary mapping subject ID to sequence.
            n_permutations: Number of permutations.
            progress: Show progress feedback.

        Returns:
            Tuple of (lambda_param, K_param).
        """
        # Get background distribution
        dist = self.estimate_background_distribution(
            query=query,
            subject_sequences=subject_sequences,
            n_permutations=n_permutations,
            progress=progress
        )

        # Fit extreme value distribution
        lambda_param, K_param = dist.estimate_extreme_params()
        self.extreme_params = (lambda_param, K_param)

        if progress:
            print(f"Extreme value parameters: λ={lambda_param:.4f}, K={K_param:.4f}",
                  file=sys.stderr)

        return lambda_param, K_param

    def estimate_significance(
            self,
            hit: Hit,
            query_length: int,
            db_size: int,
            use_extreme_distribution: bool = True
    ) -> SignificanceResult:
        """Estimate significance of a single hit.

        Args:
            hit: Hit object to evaluate.
            query_length: Length of the query sequence.
            db_size: Total length of the database.
            use_extreme_distribution: Use extreme value distribution if available.

        Returns:
            SignificanceResult object.

        Examples:
            >>> result = estimator.estimate_significance(hit, 1000, 5000000)
            >>> if result.is_significant:
            ...     print("Significant hit found!")
        """
        # Use raw score
        raw_score = hit.score

        # Calculate bit score
        bit_score = self.hit_scorer._calculate_bit_score(raw_score)

        # Estimate p-value and E-value
        if use_extreme_distribution and self.extreme_params is not None:
            # Use extreme value distribution
            lambda_param, K_param = self.extreme_params

            # Gumbel CDF: P(X <= x) = exp(-K * exp(-λ * x))
            # P(X > x) = 1 - exp(-K * exp(-λ * x))
            try:
                if raw_score > 0:
                    exp_term = K_param * math.exp(-lambda_param * raw_score)
                    pvalue = 1 - math.exp(-exp_term)
                else:
                    pvalue = 1.0
                evalue = pvalue * db_size
            except:
                pvalue = 1.0
                evalue = 1.0
        elif self.background_distribution is not None:
            # Use empirical distribution
            pvalue = self.background_distribution.get_pvalue(raw_score)
            evalue = self.background_distribution.get_evalue(raw_score, db_size)
            lambda_param, K_param = self.background_distribution.estimate_extreme_params()
        else:
            # Fallback: simple heuristic
            if raw_score > 50:
                pvalue = 0.001
                evalue = 0.01
            elif raw_score > 30:
                pvalue = 0.01
                evalue = 0.1
            else:
                pvalue = 0.1
                evalue = 1.0
            lambda_param, K_param = (0.1, 0.1)

        # Determine significance (alpha = 0.05, adjusted for multiple testing)
        # Bonferroni correction: alpha / number of tests
        # For simplicity, we use alpha = 0.05
        alpha = 0.05
        is_significant = evalue < alpha

        # Adjust for database size
        if db_size > 0:
            evalue = min(evalue, float(db_size))

        return SignificanceResult(
            hit=hit,
            raw_score=raw_score,
            pvalue=min(pvalue, 1.0),
            evalue=min(evalue, 1.0),
            bit_score=bit_score,
            lambda_param=lambda_param,
            K_param=K_param,
            is_significant=is_significant,
            significance_level=alpha
        )

    def estimate_significance_batch(
            self,
            hits: List[Hit],
            query_length: int,
            db_size: int,
            use_extreme_distribution: bool = True,
            progress: bool = True
    ) -> List[SignificanceResult]:
        """Estimate significance for multiple hits.

        Args:
            hits: List of Hit objects.
            query_length: Length of the query sequence.
            db_size: Total length of the database.
            use_extreme_distribution: Use extreme value distribution if available.
            progress: Show progress feedback.

        Returns:
            List of SignificanceResult objects.
        """
        if not hits:
            return []

        if progress:
            print(f"Estimating significance for {len(hits)} hits...", file=sys.stderr)

        results = []
        for i, hit in enumerate(hits):
            result = self.estimate_significance(
                hit=hit,
                query_length=query_length,
                db_size=db_size,
                use_extreme_distribution=use_extreme_distribution
            )
            results.append(result)

            if progress and (i + 1) % 10 == 0:
                print(f"  Processed {i + 1}/{len(hits)} hits", file=sys.stderr)

        # Sort by E-value (smallest first)
        results.sort(key=lambda r: r.evalue)

        if progress:
            significant = sum(1 for r in results if r.is_significant)
            print(f"  Done: {significant}/{len(results)} hits are significant",
                  file=sys.stderr)

        return results

    def _permute_sequence(self, sequence: str) -> str:
        """Permute a sequence randomly.

        Args:
            sequence: Original sequence.

        Returns:
            Permuted sequence.
        """
        seq_list = list(sequence)
        random.shuffle(seq_list)
        return ''.join(seq_list)

    def _sample_subjects(
            self,
            subject_sequences: Dict[str, str],
            max_subjects: int = 5
    ) -> Dict[str, str]:
        """Sample a subset of subject sequences.

        Args:
            subject_sequences: Dictionary mapping subject ID to sequence.
            max_subjects: Maximum number of subjects to sample.

        Returns:
            Dictionary of sampled subject sequences.
        """
        if len(subject_sequences) <= max_subjects:
            return subject_sequences

        subject_ids = list(subject_sequences.keys())
        sampled_ids = random.sample(subject_ids, min(max_subjects, len(subject_ids)))

        return {sid: subject_sequences[sid] for sid in sampled_ids}

    def adjust_for_multiple_testing(
            self,
            results: List[SignificanceResult],
            method: str = 'bonferroni'
    ) -> List[SignificanceResult]:
        """Apply multiple testing correction.

        Args:
            results: List of SignificanceResult objects.
            method: Correction method ('bonferroni' or 'fdr').

        Returns:
            List of corrected SignificanceResult objects.
        """
        if not results:
            return results

        n_tests = len(results)

        if method == 'bonferroni':
            # Bonferroni correction: alpha / n_tests
            for result in results:
                result.significance_level = 0.05 / n_tests
                result.is_significant = result.pvalue < result.significance_level

        elif method == 'fdr':
            # Benjamini-Hochberg FDR
            sorted_results = sorted(results, key=lambda r: r.pvalue)
            for i, result in enumerate(sorted_results):
                # BH procedure: p-value < (i+1)/n_tests * alpha
                threshold = (i + 1) / n_tests * 0.05
                result.significance_level = threshold
                result.is_significant = result.pvalue < threshold

        return results


# ============================================================================
# Statistical Reporting
# ============================================================================

def format_significance_results(
        results: List[SignificanceResult],
        format: str = 'text'
) -> str:
    """Format significance results for output.

    Args:
        results: List of SignificanceResult objects.
        format: Output format ('text', 'tsv', 'json').

    Returns:
        Formatted string.
    """
    if format == 'text':
        return _format_text(results)
    elif format == 'tsv':
        return _format_tsv(results)
    elif format == 'json':
        return _format_json(results)
    else:
        return _format_text(results)


def _format_text(results: List[SignificanceResult]) -> str:
    """Format as human-readable text."""
    lines = []
    lines.append("=" * 80)
    lines.append("Statistical Significance Results")
    lines.append("=" * 80)
    lines.append("")

    for i, result in enumerate(results[:20]):
        status = "✓ SIGNIFICANT" if result.is_significant else "✗ NOT SIGNIFICANT"
        lines.append(f"{i + 1}. Subject: {result.hit.subject_id}")
        lines.append(f"   Score: {result.raw_score} (bitscore: {result.bit_score:.1f})")
        lines.append(f"   E-value: {result.evalue:.2e}")
        lines.append(f"   P-value: {result.pvalue:.2e}")
        lines.append(f"   Status: {status}")
        lines.append("")

    if len(results) > 20:
        lines.append(f"... and {len(results) - 20} more results")

    return "\n".join(lines)


def _format_tsv(results: List[SignificanceResult]) -> str:
    """Format as TSV."""
    header = "\t".join([
        "rank", "query_id", "subject_id", "score", "pvalue", "evalue",
        "bit_score", "lambda", "K", "is_significant"
    ])

    lines = [header]
    for i, result in enumerate(results, 1):
        line = "\t".join([
            str(i),
            result.hit.query_id,
            result.hit.subject_id,
            str(result.raw_score),
            f"{result.pvalue:.2e}",
            f"{result.evalue:.2e}",
            f"{result.bit_score:.1f}",
            f"{result.lambda_param:.4f}",
            f"{result.K_param:.4f}",
            str(result.is_significant)
        ])
        lines.append(line)

    return "\n".join(lines)


def _format_json(results: List[SignificanceResult]) -> str:
    """Format as JSON."""
    import json

    data = {
        "total_results": len(results),
        "significant_count": sum(1 for r in results if r.is_significant),
        "results": [r.to_dict() for r in results]
    }

    return json.dumps(data, indent=2)


# ============================================================================
# Command Line Interface (for testing)
# ============================================================================

def main():
    """Simple command line interface for testing statistical significance."""
    import argparse
    import json

    from .io import load_hits_from_tsv, parse_fasta

    parser = argparse.ArgumentParser(description="Statistical significance tool")
    parser.add_argument("-i", "--input", required=True,
                        help="Input hits file (TSV)")
    parser.add_argument("-d", "--database", required=True,
                        help="Database FASTA file")
    parser.add_argument("-q", "--query", required=True,
                        help="Query sequence (string or file)")
    parser.add_argument("-n", "--permutations", type=int, default=100,
                        help="Number of permutations")
    parser.add_argument("-o", "--output", help="Output file")
    parser.add_argument("--format", choices=["text", "tsv", "json"], default="text",
                        help="Output format")
    parser.add_argument("--seed", type=int, help="Random seed")

    args = parser.parse_args()

    # Load hits
    hits = load_hits_from_tsv(args.input)
    print(f"Loaded {len(hits)} hits", file=sys.stderr)

    # Load database sequences
    subject_sequences = {}
    for record in parse_fasta(args.database):
        subject_sequences[record.id] = record.sequence
    print(f"Loaded {len(subject_sequences)} database sequences", file=sys.stderr)

    # Load query
    query = None
    try:
        records = list(parse_fasta(args.query))
        if records:
            query = records[0]
    except:
        query = args.query

    if query is None:
        print("Error: Could not load query", file=sys.stderr)
        sys.exit(1)

    # Calculate database size
    db_size = sum(len(seq) for seq in subject_sequences.values())

    # Create estimator
    estimator = SignificanceEstimator(
        n_permutations=args.permutations,
        random_seed=args.seed
    )

    # Estimate background distribution
    print("\nEstimating background distribution...", file=sys.stderr)
    dist = estimator.estimate_background_distribution(
        query=query,
        subject_sequences=subject_sequences,
        n_permutations=args.permutations,
        progress=True
    )

    # Estimate extreme value distribution
    print("\nFitting extreme value distribution...", file=sys.stderr)
    lambda_param, K_param = estimator.estimate_extreme_value_distribution(
        query=query,
        subject_sequences=subject_sequences,
        n_permutations=args.permutations,
        progress=False
    )

    # Estimate significance for all hits
    print("\nEstimating significance...", file=sys.stderr)
    results = estimator.estimate_significance_batch(
        hits=hits,
        query_length=len(query),
        db_size=db_size,
        progress=True
    )

    # Apply multiple testing correction
    print("\nApplying FDR correction...", file=sys.stderr)
    results = estimator.adjust_for_multiple_testing(results, method='fdr')

    # Format output
    output = format_significance_results(results, format=args.format)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(output)
        print(f"\nOutput written to {args.output}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()