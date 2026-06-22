"""Tests for k-mer index."""

import pytest
from blastmini.index import KmerIndex, build_index_from_fasta, reverse_complement, canonical_kmer
from blastmini.models import SequenceRecord


def test_reverse_complement():
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("N") == "N"


def test_canonical_kmer():
    assert canonical_kmer("ATCG") == "ATCG"
    assert canonical_kmer("CGAT") == "ATCG"  # reverse complement of CGAT is ATCG
    assert canonical_kmer("AAAA") == "AAAA"


def test_kmer_index_build(sample_records):
    idx = KmerIndex.build(sample_records, k=3)
    assert idx.k == 3
    assert idx.stats.total_sequences == 2
    # Check that some k-mers are present
    assert "ATC" in idx.index
    assert "AGC" in idx.index


def test_kmer_index_lookup(sample_records):
    idx = KmerIndex.build(sample_records, k=3)
    matches = idx.lookup("ATC")
    assert len(matches) >= 1  # seq1 has ATC at 0 and 4
    assert ("seq1", 0) in matches or ("seq1", 4) in matches


def test_kmer_index_build_from_fasta(sample_fasta_file):
    idx = build_index_from_fasta(sample_fasta_file, k=3, progress=False)
    assert idx.stats.total_sequences == 2
    assert idx.stats.unique_kmers > 0


def test_kmer_index_save_load(tmp_path, sample_records):
    idx = KmerIndex.build(sample_records, k=3)
    idx_file = tmp_path / "idx.json"
    idx.save(idx_file)
    loaded = KmerIndex.load(idx_file)
    assert loaded.k == 3
    assert loaded.index == idx.index