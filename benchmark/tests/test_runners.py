from pathlib import Path

import pytest

from benchmark.runners.base import ToolRunResult
from benchmark.runners.chainsaw import ChainsawRunner
from benchmark.runners.merizo import MerizoRunner
from benchmark.runners.sword2_orig import Sword2OriginalRunner
from benchmark.runners.sword2_rust import Sword2RustRunner
from benchmark.run_benchmark import _check_run_or_summary


def _result(command: list[str], cwd: Path, stdout: str = "") -> ToolRunResult:
    return ToolRunResult(
        command=command,
        cwd=cwd,
        returncode=0,
        runtime_s=0.1,
        peak_rss_mb=1.0,
        stdout=stdout,
        stderr="",
    )


def test_chainsaw_runner_passes_absolute_input_and_output_paths(monkeypatch, tmp_path):
    calls = []

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        Path(command[command.index("--output") + 1]).write_text("chain_id\tchopping\nA\t1-10\n")
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.chainsaw.run_timed", fake_run_timed)

    runner = ChainsawRunner(chainsaw_dir=tmp_path / "chainsaw", python=tmp_path / "py")
    runner.run(Path("benchmark/cache/chains/19hcA.pdb"), Path("benchmark/results/chainsaw/19hcA.tsv"))

    command, _ = calls[0]
    assert Path(command[3]).is_absolute()
    assert Path(command[5]).is_absolute()


def test_chainsaw_runner_preserves_venv_python_symlink(monkeypatch, tmp_path):
    calls = []
    real_python = tmp_path / "uv-python"
    real_python.touch()
    venv_python = tmp_path / "chswEnv" / "bin" / "python"
    venv_python.parent.mkdir(parents=True)
    venv_python.symlink_to(real_python)

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        Path(command[command.index("--output") + 1]).write_text("chain_id\tchopping\nA\t1-10\n")
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.chainsaw.run_timed", fake_run_timed)

    runner = ChainsawRunner(chainsaw_dir=tmp_path / "chainsaw", python=venv_python)
    runner.run(Path("benchmark/cache/chains/19hcA.pdb"), tmp_path / "out.tsv")

    command, _ = calls[0]
    assert command[0] == str(venv_python.absolute())


def test_merizo_runner_passes_absolute_input_path_and_writes_stdout(monkeypatch, tmp_path):
    calls = []

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        return _result(
            command,
            cwd,
            stdout="input\tnres\tnres_dom\tnres_ndr\tndom\tpIoU\truntime\tresult\n"
            "19hcA.pdb\t10\t10\t0\t2\t0.9\t0.1\t1-5,6-10\n",
        )

    monkeypatch.setattr("benchmark.runners.merizo.run_timed", fake_run_timed)

    output = tmp_path / "out.tsv"
    runner = MerizoRunner(merizo_dir=tmp_path / "merizo", python=tmp_path / "py")
    runner.run(Path("benchmark/cache/chains/19hcA.pdb"), output)

    command, _ = calls[0]
    assert Path(command[3]).is_absolute()
    assert output.read_text().startswith("input\tnres")


def test_merizo_runner_preserves_venv_python_symlink(monkeypatch, tmp_path):
    calls = []
    real_python = tmp_path / "uv-python"
    real_python.touch()
    venv_python = tmp_path / "merizo" / "bin" / "python"
    venv_python.parent.mkdir(parents=True)
    venv_python.symlink_to(real_python)

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        return _result(
            command,
            cwd,
            stdout="input\tnres\tnres_dom\tnres_ndr\tndom\tpIoU\truntime\tresult\n"
            "19hcA.pdb\t10\t10\t0\t2\t0.9\t0.1\t1-5,6-10\n",
        )

    monkeypatch.setattr("benchmark.runners.merizo.run_timed", fake_run_timed)

    runner = MerizoRunner(merizo_dir=tmp_path / "merizo-src", python=venv_python)
    runner.run(Path("benchmark/cache/chains/19hcA.pdb"), tmp_path / "out.tsv")

    command, _ = calls[0]
    assert command[0] == str(venv_python.absolute())


def test_merizo_cuda_runner_checks_cuda_before_launch(monkeypatch, tmp_path):
    checks = []
    calls = []

    def fake_ensure_cuda_available(python, cwd):
        checks.append((python, cwd))

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        return _result(
            command,
            cwd,
            stdout="input\tnres\tnres_dom\tnres_ndr\tndom\tpIoU\truntime\tresult\n"
            "19hcA.pdb\t10\t10\t0\t2\t0.9\t0.1\t1-5,6-10\n",
        )

    monkeypatch.setattr("benchmark.runners.merizo.ensure_cuda_available", fake_ensure_cuda_available)
    monkeypatch.setattr("benchmark.runners.merizo.run_timed", fake_run_timed)

    python = tmp_path / "py"
    merizo_dir = tmp_path / "merizo"
    runner = MerizoRunner(merizo_dir=merizo_dir, python=python, device="cuda")
    runner.run(Path("benchmark/cache/chains/19hcA.pdb"), tmp_path / "out.tsv")

    assert checks == [(python, merizo_dir)]
    assert calls


def test_merizo_batch_runner_passes_multiple_inputs_and_splits_stdout(monkeypatch, tmp_path):
    calls = []
    input_a = tmp_path / "inputs" / "19hcA.pdb"
    input_b = tmp_path / "inputs" / "1a59A.pdb"
    input_a.parent.mkdir()
    input_a.write_text("ATOM\n")
    input_b.write_text("ATOM\n")

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        return _result(
            command,
            cwd,
            stdout="input\tnres\tnres_dom\tnres_ndr\tndom\tpIoU\truntime\tresult\n"
            "19hcA.pdb\t10\t10\t0\t2\t0.9\t0.11\t1-5,6-10\n"
            "1a59A.pdb\t12\t12\t0\t1\t0.8\t0.22\t1-12\n",
        )

    monkeypatch.setattr("benchmark.runners.merizo.run_timed", fake_run_timed)

    output_a = tmp_path / "raw" / "19hcA.tsv"
    output_b = tmp_path / "raw" / "1a59A.tsv"
    runner = MerizoRunner(merizo_dir=tmp_path / "merizo", python=tmp_path / "py")
    result = runner.run_batch(
        [
            (input_a, output_a, "A"),
            (input_b, output_b, "A"),
        ]
    )

    command, _ = calls[0]
    input_index = command.index("-i")
    assert Path(command[input_index + 1]).is_absolute()
    assert Path(command[input_index + 2]).is_absolute()
    assert [Path(command[input_index + 1]).name, Path(command[input_index + 2]).name] == [
        "19hcA.pdb",
        "1a59A.pdb",
    ]
    assert output_a.read_text().splitlines() == [
        "input\tnres\tnres_dom\tnres_ndr\tndom\tpIoU\truntime\tresult",
        "19hcA.pdb\t10\t10\t0\t2\t0.9\t0.11\t1-5,6-10",
    ]
    assert output_b.read_text().splitlines() == [
        "input\tnres\tnres_dom\tnres_ndr\tndom\tpIoU\truntime\tresult",
        "1a59A.pdb\t12\t12\t0\t1\t0.8\t0.22\t1-12",
    ]
    assert result is not None


def test_merizo_runner_fails_when_process_succeeds_without_tsv(monkeypatch, tmp_path):
    def fake_run_timed(command, cwd):
        return _result(command, cwd, stdout="")

    monkeypatch.setattr("benchmark.runners.merizo.run_timed", fake_run_timed)

    runner = MerizoRunner(merizo_dir=tmp_path / "merizo", python=tmp_path / "py")
    with pytest.raises(RuntimeError, match="produced no TSV"):
        runner.run(Path("benchmark/cache/chains/missing.pdb"), tmp_path / "missing.tsv")


def test_merizo_runner_does_not_reuse_stale_tsv(monkeypatch, tmp_path):
    output = tmp_path / "stale.tsv"
    output.write_text("old\tdata\n")

    def fake_run_timed(command, cwd):
        assert not output.exists()
        return _result(command, cwd, stdout="")

    monkeypatch.setattr("benchmark.runners.merizo.run_timed", fake_run_timed)

    runner = MerizoRunner(merizo_dir=tmp_path / "merizo", python=tmp_path / "py")
    with pytest.raises(RuntimeError, match="produced no TSV"):
        runner.run(Path("benchmark/cache/chains/19hcA.pdb"), output)
    assert not output.exists()


def test_chainsaw_runner_does_not_reuse_stale_tsv(monkeypatch, tmp_path):
    output = tmp_path / "stale.tsv"
    output.write_text("old\tdata\n")

    def fake_run_timed(command, cwd):
        assert not output.exists()
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.chainsaw.run_timed", fake_run_timed)

    runner = ChainsawRunner(chainsaw_dir=tmp_path / "chainsaw", python=tmp_path / "py")
    with pytest.raises(RuntimeError, match="produced no TSV"):
        runner.run(Path("benchmark/cache/chains/19hcA.pdb"), output)
    assert not output.exists()


def test_chainsaw_batch_runner_uses_structure_directory_and_splits_output(monkeypatch, tmp_path):
    calls = []
    input_a = tmp_path / "inputs" / "19hcA.pdb"
    input_b = tmp_path / "inputs" / "1a59A.pdb"
    input_a.parent.mkdir()
    input_a.write_text("ATOM\n")
    input_b.write_text("ATOM\n")

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        output_path = Path(command[command.index("--output") + 1])
        output_path.write_text(
            "chain_id\tsequence_md5\tnres\tndom\tchopping\tconfidence\ttime_sec\n"
            "19hcA\tmd5a\t10\t2\t1-5|6-10\t0.9\t0.31\n"
            "1a59A\tmd5b\t12\t1\t1-12\t0.8\t0.42\n"
        )
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.chainsaw.run_timed", fake_run_timed)

    output_a = tmp_path / "raw" / "19hcA.tsv"
    output_b = tmp_path / "raw" / "1a59A.tsv"
    runner = ChainsawRunner(chainsaw_dir=tmp_path / "chainsaw", python=tmp_path / "py")
    result = runner.run_batch(
        [
            (input_a, output_a),
            (input_b, output_b),
        ],
        stage_dir=tmp_path / "stage",
        combined_output=tmp_path / "combined.tsv",
    )

    command, _ = calls[0]
    assert "--structure_directory" in command
    stage_dir = Path(command[command.index("--structure_directory") + 1])
    assert stage_dir.is_absolute()
    assert (stage_dir / "19hcA.pdb").exists()
    assert (stage_dir / "1a59A.pdb").exists()
    assert output_a.read_text().splitlines() == [
        "chain_id\tsequence_md5\tnres\tndom\tchopping\tconfidence\ttime_sec",
        "19hcA\tmd5a\t10\t2\t1-5|6-10\t0.9\t0.31",
    ]
    assert output_b.read_text().splitlines() == [
        "chain_id\tsequence_md5\tnres\tndom\tchopping\tconfidence\ttime_sec",
        "1a59A\tmd5b\t12\t1\t1-12\t0.8\t0.42",
    ]
    assert result is not None


def test_chainsaw_runner_can_require_cuda_before_launch(monkeypatch, tmp_path):
    checks = []
    calls = []

    def fake_ensure_cuda_available(python, cwd):
        checks.append((python, cwd))

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        Path(command[command.index("--output") + 1]).write_text("chain_id\tchopping\nA\t1-10\n")
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.chainsaw.ensure_cuda_available", fake_ensure_cuda_available)
    monkeypatch.setattr("benchmark.runners.chainsaw.run_timed", fake_run_timed)

    python = tmp_path / "py"
    chainsaw_dir = tmp_path / "chainsaw"
    runner = ChainsawRunner(chainsaw_dir=chainsaw_dir, python=python, require_cuda=True)
    runner.run(Path("benchmark/cache/chains/19hcA.pdb"), tmp_path / "out.tsv")

    assert checks == [(python, chainsaw_dir)]
    assert calls


def test_sword2_rust_runner_passes_absolute_paths(monkeypatch, tmp_path):
    calls = []

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.sword2_rust.run_timed", fake_run_timed)

    runner = Sword2RustRunner(repo_dir=tmp_path / "repo", binary=tmp_path / "sword2")
    runner.run(Path("benchmark/cache/chains/19hcA.pdb"), Path("benchmark/results/sword2-rust/19hcA"))

    command, _ = calls[0]
    assert Path(command[2]).is_absolute()
    assert Path(command[4]).is_absolute()


def test_sword2_rust_runner_cleans_stale_output_directory(monkeypatch, tmp_path):
    def fake_run_timed(command, cwd):
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.sword2_rust.run_timed", fake_run_timed)

    output_dir = tmp_path / "raw" / "sword2-rust" / "19hcA"
    stale_file = output_dir / "old" / "summary.json"
    stale_file.parent.mkdir(parents=True)
    stale_file.write_text("{}")

    runner = Sword2RustRunner(repo_dir=tmp_path / "repo", binary=tmp_path / "sword2")
    runner.run(Path("benchmark/cache/chains/19hcA.pdb"), output_dir)

    assert not stale_file.exists()
    assert output_dir.exists()


def test_sword2_rust_runner_passes_experiment_env(monkeypatch, tmp_path):
    calls = []

    def fake_run_timed(command, cwd, env=None):
        calls.append((command, cwd, env))
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.sword2_rust.run_timed", fake_run_timed)

    runner = Sword2RustRunner(
        repo_dir=tmp_path / "repo",
        binary=tmp_path / "sword2",
        experiments="semantic-fixes",
    )
    runner.run(Path("benchmark/cache/chains/19hcA.pdb"), tmp_path / "out")

    assert calls[0][2] == {"SWORD2_EXPERIMENTS": "semantic-fixes"}


def test_sword2_original_runner_uses_configured_conda_prefix_and_absolute_output(monkeypatch, tmp_path):
    calls = []

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.sword2_orig.run_timed", fake_run_timed)

    conda_env = Path("/opt/apps/pkgs/anaconda3/2023.09-0/kvhhznl/envs/sword2")
    runner = Sword2OriginalRunner(repo_dir=tmp_path / "orig", conda_env=conda_env)
    runner.run("19hc", Path("benchmark/results/sword2-original/19hcA"))

    command, _ = calls[0]
    assert command[:4] == ["conda", "run", "-p", str(conda_env)]
    assert Path(command[command.index("-o") + 1]).is_absolute()


def test_sword2_original_runner_can_use_local_structure_input(monkeypatch, tmp_path):
    calls = []

    def fake_run_timed(command, cwd):
        calls.append((command, cwd))
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.sword2_orig.run_timed", fake_run_timed)

    runner = Sword2OriginalRunner(repo_dir=tmp_path / "orig")
    runner.run(
        "19hc",
        Path("benchmark/results/sword2-original/19hcA"),
        structure_file=Path("benchmark/cache/chains/19hcA.pdb"),
        chain_id="A",
    )

    command, _ = calls[0]
    assert "-p" not in command[4:]
    assert "-i" in command
    assert Path(command[command.index("-i") + 1]).is_absolute()
    assert command[command.index("-c") + 1] == "A"


def test_sword2_original_runner_cleans_stale_output_directory(monkeypatch, tmp_path):
    def fake_run_timed(command, cwd):
        return _result(command, cwd)

    monkeypatch.setattr("benchmark.runners.sword2_orig.run_timed", fake_run_timed)

    output_dir = tmp_path / "raw" / "sword2-original" / "19hcA"
    stale_file = output_dir / "old" / "partial.txt"
    stale_file.parent.mkdir(parents=True)
    stale_file.write_text("stale")

    runner = Sword2OriginalRunner(repo_dir=tmp_path / "orig")
    runner.run("19hc", output_dir, structure_file=Path("benchmark/cache/chains/19hcA.pdb"))

    assert not stale_file.exists()
    assert output_dir.exists()


def test_nonzero_original_run_is_usable_when_summary_exists(tmp_path):
    raw_dir = tmp_path / "raw"
    summary = raw_dir / "nested" / "SWORD2_summary.json"
    summary.parent.mkdir(parents=True)
    summary.write_text("{}")
    result = ToolRunResult(
        command=["cmd"],
        cwd=tmp_path,
        returncode=1,
        runtime_s=0.1,
        peak_rss_mb=None,
        stdout="",
        stderr="post-summary failure",
    )

    assert _check_run_or_summary(result, "sword2-original", raw_dir, "SWORD2_summary.json") == summary
