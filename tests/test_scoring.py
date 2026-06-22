"""Tests for scoring and formatting."""

import pytest
from blastmini.scoring import HitScorer, ScoredHit, format_hits_as_text, format_hits_as_tsv, format_hits_as_json
from blastmini.models import Hit, AlignmentConfig


def test_hit_scorer_score_hit():
    scorer = HitScorer()
    hit = Hit("Q1", "S1", score=100, identity_percent=95.0, alignment_length=100,
              query_start=0, query_end=100, subject_start=0, subject_end=100)
    scored = scorer.score_hit(hit, query_length=200)
    assert scored.raw_score == 100
    assert scored.identity_percent == 95.0
    assert scored.coverage_percent == 50.0


def test_hit_scorer_score_hits():
    scorer = HitScorer()
    hits = [Hit("Q1", "S1", 100, 95.0, 100), Hit("Q1", "S2", 80, 90.0, 80)]
    scored = scorer.score_hits(hits, query_length=200, db_total_length=500, progress=False)
    assert len(scored) == 2
    assert scored[0].rank == 1
    assert scored[0].raw_score == 100


def test_format_hits_as_text():
    hit = Hit("Q1", "S1", 100, 95.0, 100)
    scored = ScoredHit(hit, raw_score=100, identity_percent=95.0, coverage_percent=50.0)
    text = format_hits_as_text([scored], top_n=1)
    assert "S1" in text
    assert "95.0%" in text


def test_format_hits_as_tsv():
    hit = Hit("Q1", "S1", 100, 95.0, 100)
    scored = ScoredHit(hit, raw_score=100, identity_percent=95.0, coverage_percent=50.0)
    tsv = format_hits_as_tsv([scored])
    assert "Q1" in tsv
    assert "100" in tsv