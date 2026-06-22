"""Tests for I/O functions."""

import pytest
import gzip
from blastmini.io import parse_fasta, save_index, load_index, save_hits_to_tsv, load_hits_from_tsv
from blastmini.models import SequenceRecord, Hit


def test_parse_fasta(sample_fasta_file):
    records = list(parse_fasta(sample_fasta_file))
    assert len(records) == 2
    assert records[0].id == "seq1"
    assert records[0].sequence == "ATCGATCGATCGATCG"
    assert records[1].id == "seq2"


def test_parse_fasta_gz(tmp_path, sample_fasta_content):
    gz_path = tmp_path / "sample.fa.gz"
    with gzip.open(gz_path, "wt") as f:
        f.write(sample_fasta_content)
    records = list(parse_fasta(gz_path))
    assert len(records) == 2


def test_parse_fasta_file_not_found():
    with pytest.raises(FileNotFoundError):
        list(parse_fasta("nonexistent.fa"))


def test_save_load_index(tmp_path):
    index = {"ATCG": [("seq1", 0), ("seq1", 4)], "GCTA": [("seq2", 0)]}
    index_file = tmp_path / "index.json"
    save_index(index, index_file)
    loaded = load_index(index_file)
    assert loaded == index


def test_save_load_hits(tmp_path):
    hits = [Hit("Q1", "S1", 100, 95.0, 200, 0, 200, 0, 200,evalue=1e-5, bit_score=50.0)]
    tsv_file = tmp_path / "hits.tsv"
    save_hits_to_tsv(hits, tsv_file)
    loaded = load_hits_from_tsv(tsv_file)
    assert len(loaded) == 1
    assert loaded[0].query_id == "Q1"