"""pytest configuration and shared fixtures for blastmini tests."""

import pytest

from blastmini.models import AlignmentConfig, SequenceRecord


@pytest.fixture
def sample_fasta_content():
    """Return a sample FASTA content with two sequences."""
    return """>seq1
ATCGATCGATCGATCG
>seq2
GCTAGCTAGCTAGCTA
"""


@pytest.fixture
def sample_fasta_file(tmp_path, sample_fasta_content):
    """Create a temporary FASTA file with sample content."""
    file_path = tmp_path / "sample.fa"
    file_path.write_text(sample_fasta_content)
    return file_path


@pytest.fixture
def sample_records():
    """Return a list of SequenceRecord objects for testing."""
    return [
        SequenceRecord("seq1", "ATCGATCGATCGATCG"),
        SequenceRecord("seq2", "GCTAGCTAGCTAGCTA"),
    ]


@pytest.fixture
def default_config():
    """Return a default AlignmentConfig."""
    return AlignmentConfig()
