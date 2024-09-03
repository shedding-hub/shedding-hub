from pathlib import Path
import subprocess


def test_summarize(tmp_path: Path) -> None:
    # Summarize Wölfel et al. (2020) which is a well established and tested dataset.
    # This test may need to be updated if we make changes to the dataset.
    tmp_file = tmp_path / "output.md"
    subprocess.check_output(
        [
            "python",
            ".github/workflows/summarize.py",
            "-o",
            str(tmp_file),
            "data/woelfel2020.yaml",
        ]
    )
    result = tmp_file.read_text()
    expected = """
             n_samples  n_unique_participants  n_negative  n_positive  n_quantified
sputum             147                      9          24           0           123
stool               82                      9          13           0            69
throat_swab        153                      9          57           0            96
"""
    for line in expected.splitlines():
        assert line.strip() in result
