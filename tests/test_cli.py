"""Tests for the command-line interface (cli.py)."""

import json
import subprocess
import sys


def run_cli(*args: str) -> subprocess.CompletedProcess:
    """Run the CLI with given arguments using subprocess."""
    cmd = [sys.executable, "-m", "blastmini.cli"] + list(args)
    return subprocess.run(cmd, capture_output=True, text=True)


def test_cli_info():
    """Test the 'info' command."""
    result = run_cli("info")
    assert result.returncode == 0
    assert "blastmini v" in result.stdout
    assert "k-mer indexing" in result.stdout
    assert "X-dropoff extension" in result.stdout


def test_cli_version():
    """Test the '--version' flag."""
    result = run_cli("--version")
    assert result.returncode == 0
    assert "blastmini" in result.stdout


def test_cli_build(sample_fasta_file, tmp_path):
    """Test the 'build' command."""
    index_file = tmp_path / "index.json"
    result = run_cli(
        "build",
        "-i", str(sample_fasta_file),
        "-o", str(index_file),
        "-k", "3",
        "-q"  # quiet mode
    )
    assert result.returncode == 0
    assert index_file.exists()
    # Check that index file contains valid JSON
    with open(index_file) as f:
        data = json.load(f)
    assert len(data) > 0


def test_cli_build_with_filters(sample_fasta_file, tmp_path):
    """Test the 'build' command with occurrence filters."""
    index_file = tmp_path / "index.json"
    result = run_cli(
        "build",
        "-i", str(sample_fasta_file),
        "-o", str(index_file),
        "-k", "3",
        "--min-occurrences", "1",
        "--max-occurrences", "5",
        "-q"
    )
    assert result.returncode == 0
    assert index_file.exists()


def test_cli_stats(sample_fasta_file, tmp_path):
    """Test the 'stats' command."""
    # First build an index
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    # Then get stats
    result = run_cli("stats", "-i", str(index_file))
    assert result.returncode == 0
    assert "k-mer size: 3" in result.stdout
    assert "Sequences: 2" in result.stdout
    assert "Unique k-mers:" in result.stdout


def test_cli_stats_verbose(sample_fasta_file, tmp_path):
    """Test the 'stats' command with verbose flag."""
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    result = run_cli("stats", "-i", str(index_file), "-v")
    assert result.returncode == 0
    assert "Average positions per k-mer" in result.stdout
    assert "Most frequent k-mers:" in result.stdout


def test_cli_search(sample_fasta_file, tmp_path):
    """Test the 'search' command."""
    # Build index
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    # Search
    output_file = tmp_path / "results.tsv"
    result = run_cli(
        "search",
        # Use same file as query (contains seq1 and seq2)
        "-q", str(sample_fasta_file),
        "-d", str(index_file),
        "--database-fasta", str(sample_fasta_file),
        "-o", str(output_file),
        "--top", "5",
        "--min-score", "0",
        "--format", "tsv",
        "--quiet"
    )
    assert result.returncode == 0
    assert output_file.exists()
    # Check that output has content
    with open(output_file) as f:
        content = f.read()
    assert len(content) > 0
    assert "query_id" in content  # Header present


def test_cli_search_text_output(sample_fasta_file, tmp_path):
    """Test the 'search' command with text output (stdout)."""
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    result = run_cli(
        "search",
        "-q", str(sample_fasta_file),
        "-d", str(index_file),
        "--database-fasta", str(sample_fasta_file),
        "--top", "5",
        "--format", "text",
        "--quiet"
    )
    assert result.returncode == 0
    assert "BLAST Search Results" in result.stdout


def test_cli_search_json_output(sample_fasta_file, tmp_path):
    """Test the 'search' command with JSON output."""
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    output_file = tmp_path / "results.json"
    result = run_cli(
        "search",
        "-q", str(sample_fasta_file),
        "-d", str(index_file),
        "--database-fasta", str(sample_fasta_file),
        "-o", str(output_file),
        "--top", "5",
        "--format", "json",
        "--quiet"
    )
    assert result.returncode == 0
    assert output_file.exists()
    with open(output_file) as f:
        data = json.load(f)
    assert "total_hits" in data
    assert "hits" in data


def test_cli_search_bed_output(sample_fasta_file, tmp_path):
    """Test the 'search' command with BED output."""
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    output_file = tmp_path / "results.bed"
    result = run_cli(
        "search",
        "-q", str(sample_fasta_file),
        "-d", str(index_file),
        "--database-fasta", str(sample_fasta_file),
        "-o", str(output_file),
        "--top", "5",
        "--format", "bed",
        "--quiet"
    )
    assert result.returncode == 0
    assert output_file.exists()
    with open(output_file) as f:
        content = f.read()
    assert "#track" in content
    assert "blastmini" in content


def test_cli_view(sample_fasta_file, tmp_path):
    """Test the 'view' command."""
    # Build index and search first
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    results_file = tmp_path / "results.tsv"
    run_cli(
        "search",
        "-q", str(sample_fasta_file),
        "-d", str(index_file),
        "--database-fasta", str(sample_fasta_file),
        "-o", str(results_file),
        "--top", "5",
        "--format", "tsv",
        "--quiet"
    )

    # View results
    result = run_cli("view", "-r", str(results_file), "--top", "3")
    assert result.returncode == 0
    assert "BLAST Search Results" in result.stdout


def test_cli_view_with_alignment(sample_fasta_file, tmp_path):
    """Test the 'view' command with alignment visualization."""
    # Build index and search first
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    results_file = tmp_path / "results.tsv"
    run_cli(
        "search",
        "-q", str(sample_fasta_file),
        "-d", str(index_file),
        "--database-fasta", str(sample_fasta_file),
        "-o", str(results_file),
        "--top", "5",
        "--format", "tsv",
        "--quiet"
    )

    result = run_cli("view", "-r", str(results_file), "--show-alignment")
    assert result.returncode == 0
    assert "Alignment Visualization" in result.stdout
    assert "Query:" in result.stdout or "No alignment data" in result.stdout


def test_cli_view_json_format(sample_fasta_file, tmp_path):
    """Test the 'view' command with JSON format."""
    # Build index and search first
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    results_file = tmp_path / "results.tsv"
    run_cli(
        "search",
        "-q", str(sample_fasta_file),
        "-d", str(index_file),
        "--database-fasta", str(sample_fasta_file),
        "-o", str(results_file),
        "--top", "5",
        "--format", "tsv",
        "--quiet"
    )

    result = run_cli("view", "-r", str(results_file), "--format", "json")
    assert result.returncode == 0
    data = json.loads(result.stdout)
    assert "total_hits" in data


def test_cli_significance(sample_fasta_file, tmp_path):
    """Test the 'significance' command (with small permutation count)."""
    # Build index and search first
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    results_file = tmp_path / "results.tsv"
    run_cli(
        "search",
        "-q", str(sample_fasta_file),
        "-d", str(index_file),
        "--database-fasta", str(sample_fasta_file),
        "-o", str(results_file),
        "--top", "5",
        "--format", "tsv",
        "--quiet"
    )

    # Estimate significance with small number of permutations
    output_file = tmp_path / "significance.txt"
    result = run_cli(
        "significance",
        "-r", str(results_file),
        "-d", str(sample_fasta_file),
        "-q", str(sample_fasta_file),
        "-n", "2",  # Very small for testing
        "-o", str(output_file),
        "--quiet"  # quiet mode
    )
    assert result.returncode == 0
    assert output_file.exists()
    with open(output_file) as f:
        content = f.read()
    assert "Statistical Significance Results" in content or "SIGNIFICANT" in content


def test_cli_significance_with_correction(sample_fasta_file, tmp_path):
    """Test the 'significance' command with FDR correction."""
    # Build index and search first
    index_file = tmp_path / "index.json"
    run_cli("build", "-i", str(sample_fasta_file),
            "-o", str(index_file), "-k", "3", "-q")

    results_file = tmp_path / "results.tsv"
    run_cli(
        "search",
        "-q", str(sample_fasta_file),
        "-d", str(index_file),
        "--database-fasta", str(sample_fasta_file),
        "-o", str(results_file),
        "--top", "5",
        "--format", "tsv",
        "--quiet"
    )

    result = run_cli(
        "significance",
        "-r", str(results_file),
        "-d", str(sample_fasta_file),
        "-q", str(sample_fasta_file),
        "-n", "2",
        "--correction", "fdr",
        "--quiet"
    )
    assert result.returncode == 0
    # Check that output contains expected format (text by default)
    assert "Statistical Significance Results" in result.stdout or "SIGNIFICANT" in result.stdout


def test_cli_error_file_not_found():
    """Test CLI error handling for missing files."""
    result = run_cli("build", "-i", "nonexistent.fa", "-o", "out.json")
    assert result.returncode == 1
    assert "ERROR:" in result.stderr


def test_cli_error_missing_required_argument():
    """Test CLI error handling for missing required arguments."""
    result = run_cli("build")  # Missing -i and -o
    assert result.returncode == 2
    assert "error:" in result.stderr


def test_cli_help():
    """Test the help message."""
    result = run_cli("--help")
    assert result.returncode == 0
    assert "blastmini" in result.stdout
    assert "build" in result.stdout
    assert "search" in result.stdout
    assert "stats" in result.stdout
    assert "view" in result.stdout
    assert "significance" in result.stdout


def test_cli_command_help():
    """Test help for individual commands."""
    result = run_cli("build", "--help")
    assert result.returncode == 0
    assert "-i" in result.stdout
    assert "-o" in result.stdout
    assert "--kmer" in result.stdout

    result = run_cli("search", "--help")
    assert result.returncode == 0
    assert "-q" in result.stdout
    assert "-d" in result.stdout
