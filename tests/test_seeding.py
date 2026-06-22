"""Tests for seed finding."""

import pytest
from blastmini.seeding import SeedFinder, Seed, SeedCluster
from blastmini.index import KmerIndex
from blastmini.models import SequenceRecord, AlignmentConfig


def test_seed_creation():
    seed = Seed("Q1", "S1", 10, 100, "ATCG")
    assert seed.query_id == "Q1"
    assert seed.subject_id == "S1"
    assert seed.query_pos == 10
    assert seed.subject_pos == 100


def test_seed_cluster():
    cluster = SeedCluster("Q1", "S1")
    seed1 = Seed("Q1", "S1", 10, 100)
    seed2 = Seed("Q1", "S1", 15, 105)
    cluster.add_seed(seed1)
    cluster.add_seed(seed2)
    assert len(cluster) == 2
    assert cluster.query_range == (10, 15)
    assert cluster.density > 0


def test_seed_finder_find_seeds(sample_records):
    idx = KmerIndex.build(sample_records, k=3)
    finder = SeedFinder(idx)
    query = SequenceRecord("query", "ATCGATCG")
    seeds = finder.find_seeds(query, max_seeds=10, progress=False)
    # Should find seeds for ATC and GAT etc.
    assert len(seeds) > 0
    # Check that all seeds have correct subject IDs
    for seed in seeds:
        assert seed.subject_id in ["seq1", "seq2"]


def test_seed_finder_cluster(sample_records):
    idx = KmerIndex.build(sample_records, k=3)
    finder = SeedFinder(idx)
    query = SequenceRecord("query", "ATCGATCG")
    seeds = finder.find_seeds(query, progress=False)
    clusters = finder.cluster_seeds(seeds)
    # There should be at least one cluster
    assert len(clusters) >= 1
    # Check cluster density
    for cluster in clusters:
        assert cluster.density >= 0


def test_seed_finder_find_best_seeds(sample_records):
    idx = KmerIndex.build(sample_records, k=3)
    finder = SeedFinder(idx)
    query = SequenceRecord("query", "ATCGATCG")
    best = finder.find_best_seeds(query, top_n=5)
    assert len(best) <= 5
    assert len(best) > 0