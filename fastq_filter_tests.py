import sys
import pytest
from Bio import SeqIO
from loguru import logger
from fastq_filter import filter_fastq


def setup_logging(log_file: str = "tests.log") -> None:
    """
    Configure loguru to write to both a file and the console.

    Args:
        log_file (str): Path to the log file. Default is 'tests.log'.
    """
    logger.remove()

    logger.add(
        sys.stderr,
        level="DEBUG",
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> [<level>{level}</level>] {message}",
        colorize=True,
    )

    logger.add(
        log_file,
        level="DEBUG",
        format="{time:YYYY-MM-DD HH:mm:ss} [{level}] {message}",
        rotation="10 MB",
        retention="7 days",
        encoding="utf-8",
    )


setup_logging()


def write_fastq(path, records: list[tuple[str, str, str]]) -> None:
    """
    Write a test FASTQ file.

    Args:
        path: pathlib.Path to write to.
        records: list of (name, sequence, quality_string) tuples.
                 Quality string must be the same length as the sequence.
    """
    with open(path, "w") as tf:
        for name, seq, qual in records:
            tf.write(f"@{name}\n{seq}\n+\n{qual}\n")


def read_fastq_ids(path) -> list[str]:
    """Return the list of sequence IDs present in a FASTQ file."""
    return [rec.id for rec in SeqIO.parse(path, "fastq")]


SEQ_SHORT = ("short", "AATTTAAACC", chr(10 + 33) * 10)
SEQ_MEDIUM = ("medium", "AT" * 10 + "GC" * 10 + "AT" * 5, chr(30 + 33) * 50)
SEQ_LONG = ("long", "GC" * 35 + "AT" * 15, chr(40 + 33) * 100)


@pytest.fixture()
def sample_fastq(tmp_path):
    """FASTQ file with three representative records."""
    path = tmp_path / "test.fastq"
    write_fastq(path, [SEQ_SHORT, SEQ_MEDIUM, SEQ_LONG])
    logger.info("Fixture test_fastq created at {}", path)
    return path


@pytest.fixture()
def output_path(tmp_path):
    """Empty output path inside a temp directory."""
    return tmp_path / "output_test.fastq"


class TestLengthFilter:
    def test_upper_bound_excludes_long_sequence(self, sample_fastq, output_path):
        """Sequences longer than the upper bound must be excluded."""
        logger.info("Start test_upper_bound_excludes_long_sequence")
        filter_fastq(str(sample_fastq), str(output_path), length_bounds=(0, 60))
        ids = read_fastq_ids(output_path)
        assert "long" not in ids
        assert "short" in ids
        assert "medium" in ids

    def test_lower_bound_excludes_short_sequence(self, sample_fastq, output_path):
        """Sequences shorter than the lower bound must be excluded."""
        logger.info("Start test_lower_bound_excludes_short_sequence")
        filter_fastq(str(sample_fastq), str(output_path), length_bounds=(20, 2**32))
        ids = read_fastq_ids(output_path)
        assert "short" not in ids
        assert "medium" in ids
        assert "long" in ids

    def test_single_int_treated_as_upper_bound(self, sample_fastq, output_path):
        """Passing a single int as length_bounds should act as (0, int)."""
        logger.info("Start test_single_int_treated_as_upper_bound")
        filter_fastq(str(sample_fastq), str(output_path), length_bounds=60)
        ids = read_fastq_ids(output_path)
        assert "long" not in ids


class TestQualityFilter:
    def test_quality_threshold_excludes_low_quality(self, sample_fastq, output_path):
        """Sequences with average quality below the threshold must be excluded."""
        logger.info("Start test_quality_threshold_excludes_low_quality")
        filter_fastq(str(sample_fastq), str(output_path), quality_threshold=25)
        ids = read_fastq_ids(output_path)
        assert "short" not in ids
        assert "medium" in ids
        assert "long" in ids


class TestGCFilter:
    def test_gc_bounds_exclude_high_gc(self, sample_fastq, output_path):
        """Sequences with GC content above the upper bound must be excluded."""
        logger.info("Start test_gc_bounds_exclude_high_gc")
        filter_fastq(str(sample_fastq), str(output_path), gc_bounds=55)
        ids = read_fastq_ids(output_path)
        assert "long" not in ids
        assert "short" in ids
        assert "medium" in ids


class TestFileIO:
    def test_output_file_is_created(self, sample_fastq, output_path):
        """Output file must be created even when all sequences pass the filter."""
        logger.info("Start test_output_file_is_created")
        filter_fastq(str(sample_fastq), str(output_path))
        assert output_path.exists(), "Output file was not created"

    def test_output_file_is_empty_when_nothing_passes(self, sample_fastq, output_path):
        """Output file must be empty (zero records) when no sequence passes."""
        logger.info("Start test_output_file_is_empty_when_nothing_passes")
        filter_fastq(str(sample_fastq), str(output_path), quality_threshold=99)
        assert read_fastq_ids(output_path) == []


class TestErrorHandling:
    def test_missing_input_file_raises_file_not_found(self, output_path):
        """filter_fastq must raise FileNotFoundError for a non-existent input."""
        logger.error(
            "Intentionally triggering FileNotFoundError (non-existent input file test)"
        )
        with pytest.raises(FileNotFoundError):
            filter_fastq("non_existent_file.fastq", str(output_path))

    def test_wrong_parameter_raises_value_error(self, sample_fastq, output_path):
        """filter_fastq must raise ValueError for a wrong value of a parameter."""
        logger.error("Intentionally triggering ValueError (invalid parameter test)")
        with pytest.raises(ValueError):
            filter_fastq(str(sample_fastq), str(output_path), length_bounds=(-30, -60))


def main():
    """Main function"""
    try:
        pytest.main([__file__, "-v"])
    except FileNotFoundError as e:
        logger.error("Input file not found: {}", e)
        sys.exit(1)
    except Exception as e:
        logger.exception("Unexpected error: {}", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
