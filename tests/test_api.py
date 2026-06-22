"""Tests for the high-level API."""

from blastmini.api import (BatchSearchResult, BlastMini, SearchResult,
                           quick_search)


def test_blastmini_from_fasta(sample_fasta_file):
    blast = BlastMini.from_fasta(sample_fasta_file, k=3, verbose=False)
    stats = blast.get_stats()
    assert stats['total_sequences'] == 2
    assert stats['kmer_size'] == 3


def test_blastmini_search(sample_fasta_file, sample_records):
    blast = BlastMini.from_fasta(sample_fasta_file, k=3, verbose=False)
    query = sample_records[0].sequence  # "ATCGATCGATCGATCG"
    result = blast.search(query, top_n=5)
    assert isinstance(result, SearchResult)
    assert result.num_hits >= 1
    # The query is identical to seq1, so should hit seq1
    assert any(hit.hit.subject_id == "seq1" for hit in result.hits)


def test_blastmini_search_with_significance(sample_fasta_file):
    blast = BlastMini.from_fasta(sample_fasta_file, k=3, verbose=False)
    # Need to estimate background first for significance
    query = blast.subject_sequences["seq1"]
    blast.estimate_background(query, n_permutations=2)
    result = blast.search(query, estimate_significance=True,
                          significance_params={'n_permutations': 2})
    # Should have hits with significance info
    if result.hits:
        assert hasattr(result.hits[0], 'evalue')
        assert hasattr(result.hits[0], 'is_significant')


def test_blastmini_search_fasta(sample_fasta_file):
    # Use same file as both database and query file (contains two sequences)
    blast = BlastMini.from_fasta(sample_fasta_file, k=3, verbose=False)
    batch = blast.search_fasta(sample_fasta_file, top_n=2)
    assert isinstance(batch, BatchSearchResult)
    assert batch.total_queries == 2
    assert batch.total_hits >= 2


def test_quick_search(sample_fasta_file):
    result = quick_search("ATCGATCG", sample_fasta_file, k=3, top_n=2)
    assert result.num_hits >= 1
