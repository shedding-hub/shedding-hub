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
            "data/woelfel2020virological/woelfel2020virological.yaml",
        ]
    )
    result = tmp_file.read_text()
    expected = """
             n_samples  n_unique_participants  n_negative  n_positive  n_quantified        min        median           max
sputum           147.0                    9.0        24.0         0.0         123.0  93.830687  29089.805626  4.473311e+08
stool             82.0                    9.0        13.0         0.0          69.0   5.736962  12422.291543  3.225616e+07
throat_swab      153.0                    9.0        57.0         0.0          96.0   3.977220    461.750274  5.835205e+08
"""  # noqa: E501
    for line in expected.splitlines():
        assert line.strip() in result
