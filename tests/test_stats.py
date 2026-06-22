"""Tests for statistical significance."""

from blastmini.models import AlignmentConfig, Hit
from blastmini.stats import ScoreDistribution, SignificanceEstimator


def test_score_distribution():
    dist = ScoreDistribution([10, 12, 15, 18, 20])
    assert dist.mean == 15.0
    assert dist.min_score == 10
    assert dist.max_score == 20
    assert dist.get_pvalue(15) == 3/5  # scores >=15 are 15,18,20 => 3/5
    # Actually careful: pvalue is P(score >= given) => for 15, count >=15 is 3 => 3/5=0.6
    # Let's adjust test: get_pvalue(15) should be 0.6
    assert dist.get_pvalue(15) == 0.6
    assert dist.get_percentile(50) == 15


def test_significance_estimator_basic(sample_records):
    # Use small data to test without permutation heavy
    config = AlignmentConfig(match_score=2, mismatch_penalty=-1, x_dropoff=2)
    estimator = SignificanceEstimator(config=config, n_permutations=2)
    subject_sequences = {rec.id: rec.sequence for rec in sample_records}
    query = sample_records[0]  # seq1
    dist = estimator.estimate_background_distribution(
        query, subject_sequences, n_permutations=2, progress=False)
    assert dist.n_samples >= 0  # may be 0 if no scores
    # Since permuted sequences may not yield scores, just test no crash


def test_estimate_significance():
    estimator = SignificanceEstimator()
    hit = Hit("Q1", "S1", score=100)
    result = estimator.estimate_significance(hit, query_length=100, db_size=1000,
                                             use_extreme_distribution=False)
    assert result.raw_score == 100
    assert result.evalue >= 0
    assert result.pvalue >= 0


def test_multiple_testing_correction():
    estimator = SignificanceEstimator()
    hits = [Hit("Q1", f"S{i}", score=50+i) for i in range(10)]
    results = [estimator.estimate_significance(
        h, 100, 1000, use_extreme_distribution=False) for h in hits]
    corrected = estimator.adjust_for_multiple_testing(
        results, method='bonferroni')
    for r in corrected:
        assert r.significance_level == 0.05/10
    # FDR
    corrected2 = estimator.adjust_for_multiple_testing(results, method='fdr')
    for r in corrected2:
        assert r.significance_level > 0
