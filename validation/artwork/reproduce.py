"""Recreate current figure assets, then compare their bytes with the manuscript."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]


def digest(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    if out == REPO or REPO in out.parents or out.exists():
        parser.error("Choose a new output directory outside the repository.")
    # The renderer imports matplotlib; keep all generated files outside the release.
    out.mkdir(parents=True)
    os.environ["MPLCONFIGDIR"] = str(out / "matplotlib-cache")
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["PYTHONOPTIMIZE"] = "0"
    sys.dont_write_bytecode = True
    provenance = json.loads((ROOT / "PROVENANCE.json").read_text(encoding="utf-8"))
    for archive in provenance.values():
        for name, source in archive["members"].items():
            if digest(ROOT / name) != source["sha256"]:
                raise RuntimeError("Changed source file " + name)

    import generate_current_figures as earlier
    earlier.configure()
    outputs = []
    for key, draw in (("packet", earlier.packet_plot), ("alpha", earlier.alpha_plot),
                      ("continuum", earlier.continuum_plot)):
        data = json.loads((ROOT / "data" / (key + ".json")).read_text(encoding="utf-8"))
        name = earlier.NAMES[key]
        draw(data, out / name)
        outputs.append(name)

    ref = ROOT.parent / "reference/tools"
    for script, name in (
        ("make_reference_figure.py", "paper1_sectionVI_fig4_reference_partition_v2.pdf"),
        ("make_packet_diagnostics_figure.py", "paper1_sectionVII_fig6b_selector_controls_v2.pdf"),
    ):
        subprocess.run([sys.executable, "-B", ref / script, "--output", out / name], check=True)
        outputs.append(name)

    name = "paper1_sectionIV_fig3_bulk_modes_symmetry_v2.pdf"
    subprocess.run([
        sys.executable, "-B", ROOT / "section4/tools/revise_figure3.py",
        ROOT / "section4/baseline/paper1_figure3_v0.pdf", out / name,
    ], check=True)
    outputs.append(name)
    name = "paper1_low_frequency_ratio.pdf"
    subprocess.run([
        sys.executable, "-B", ROOT / "low_frequency_ratio/render.py",
        "--check-sources", "--repository", REPO, "--output", out / name,
    ], check=True)
    outputs.append(name)
    checks = {}
    for name in outputs:
        expected = digest(REPO / "manuscript/source/figures" / name)
        actual = digest(out / name)
        checks[name] = {"expected": expected, "actual": actual, "match": actual == expected}
    (out / "ARTWORK_HASHES.json").write_text(json.dumps(checks, indent=2) + "\n", encoding="utf-8")
    if not all(row["match"] for row in checks.values()):
        raise SystemExit("Rendered files differ from the current publication assets. See ARTWORK_HASHES.json.")
    artwork = REPO / "manuscript/source/figures"
    historical = ROOT.parent / "reference/context/figures"
    for name in (
        "paper1_frequency_partition.pdf",
        "paper1_sectionVII_fig6a_qT_window_v2.pdf",
        "paper1_sectionVII_fig6b_selector_controls_v2.pdf",
    ):
        if digest(artwork / name) != digest(historical / name):
            raise RuntimeError("Vector comparison input differs from the current artwork: " + name)
    subprocess.run([sys.executable, "-B", ref / "verify_frequency_figure.py"], check=True)
    subprocess.run([sys.executable, "-B", ref / "verify_packet_figure.py"], check=True)
    print(f"PASS. {len(checks)} reproduced assets match byte-for-byte; frequency and packet vectors also checked.")


if __name__ == "__main__":
    main()
