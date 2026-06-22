"""Tests for data models."""

import pytest

from blastmini.models import AlignmentConfig, Hit, SequenceRecord


def test_sequence_record_creation():
    rec = SequenceRecord("test_id", "ATCG", description="test desc")
    assert rec.id == "test_id"
    assert rec.sequence == "ATCG"
    assert rec.description == "test desc"
    assert len(rec) == 4


def test_sequence_record_normalization():
    rec = SequenceRecord("id", "  atcg  ")
    assert rec.sequence == "ATCG"


def test_hit_creation():
    hit = Hit(query_id="Q1", subject_id="S1", score=100, identity_percent=95.0)
    assert hit.query_id == "Q1"
    assert hit.score == 100
    assert hit.identity_percent == 95.0


def test_hit_tsv():
    hit = Hit("Q1", "S1", 100, 95.0, 200, 0, 200, 0, 200)
    tsv = hit.to_tsv()
    assert "Q1" in tsv
    assert "S1" in tsv
    assert "100" in tsv


def test_alignment_config_validation():
    with pytest.raises(ValueError):
        AlignmentConfig(kmer_size=0)
    with pytest.raises(ValueError):
        AlignmentConfig(match_score=0)
    with pytest.raises(ValueError):
        AlignmentConfig(mismatch_penalty=1)
    with pytest.raises(ValueError):
        AlignmentConfig(x_dropoff=0)
