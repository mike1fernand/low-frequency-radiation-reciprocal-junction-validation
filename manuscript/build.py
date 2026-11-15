"""Compile both Paper 1 documents with pdfLaTeX. Requires Python 3.8 or later."""
import argparse
from contextlib import contextmanager
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
DOCUMENTS = ("paper1_main", "paper1_supplement")
REVTEX_FILES = ("revtex4-2.cls", "aps4-2.rtx", "aps10pt4-2.rtx", "revsymb4-2.sty")
UNRESOLVED = re.compile(
    r"There were undefined (?:references|citations)|"
    r"(?:Reference|Citation) [^\n]*?undefined|"
    r"Label\(s\) may have changed|Rerun to get (?:outlines|cross-references) right|"
    r"There were multiply-defined labels",
    re.IGNORECASE,
)


@contextmanager
def workspace(keep):
    if keep is None:
        with tempfile.TemporaryDirectory(prefix="paper1-latex-") as directory:
            yield Path(directory)
    else:
        path = keep.resolve()
        path.mkdir(parents=True, exist_ok=False)
        yield path


def stage_inputs(work):
    sources = work / "inputs"
    shutil.copytree(ROOT / "source", sources)
    with ZipFile(ROOT / "third_party" / "revtex-4.2f-tds.zip") as archive:
        for name in REVTEX_FILES:
            (sources / name).write_bytes(archive.read("tex/latex/revtex/" + name))
    return sources


def project_inputs(recorder, sources):
    found = set()
    for record in recorder.read_text(encoding="utf-8", errors="replace").splitlines():
        if not record.startswith("INPUT "):
            continue
        path = (sources / record[6:]).resolve()
        try:
            relative = path.relative_to(sources).as_posix()
        except ValueError:
            continue
        if relative in REVTEX_FILES:
            found.add("third_party/revtex-4.2f-tds.zip!tex/latex/revtex/" + relative)
        else:
            found.add("source/" + relative)
    return sorted(found)


def compile_one(engine, stem, sources, work):
    job = work / stem
    job.mkdir()
    cmd = [
        engine, "-no-shell-escape", "-interaction=nonstopmode", "-halt-on-error",
        "-file-line-error", "-recorder", "-output-directory=" + str(job), stem + ".tex",
    ]
    for number in range(1, 4):
        print("Compiling {} (pass {}/3)".format(stem, number), flush=True)
        run = subprocess.run(cmd, cwd=sources, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if run.returncode:
            details = run.stdout.decode("utf-8", errors="replace").splitlines()[-30:]
            raise RuntimeError("{} failed on pass {}.\n{}".format(stem, number, "\n".join(details)))
    log = (job / (stem + ".log")).read_text(encoding="utf-8", errors="replace")
    if UNRESOLVED.search(log):
        raise RuntimeError("{} still has unresolved or duplicate references. Rebuild with --work for its full log.".format(stem))
    inputs = project_inputs(job / (stem + ".fls"), sources)
    (job / (stem + "-inputs.txt")).write_text("\n".join(inputs) + "\n", encoding="utf-8")
    return job


def publish_file(source, target):
    # Each output takes its final name only once its complete copy is on the destination volume.
    with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as file:
        staged = Path(file.name)
    try:
        shutil.copyfile(source, staged)
        staged.replace(target)
    finally:
        if staged.exists():
            staged.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "build", help="PDFs and final logs (default: build/)")
    parser.add_argument("--work", type=Path, help="retain full TeX recorder paths and auxiliaries in a new directory")
    args = parser.parse_args()
    engine = shutil.which("pdflatex")
    if not engine:
        parser.exit(1, "Install pdfLaTeX and the standard packages listed in README.md, then run this command again.\n")
    out = args.out.resolve()
    for path in (out, args.work.resolve() if args.work else None):
        if path and (path == ROOT / "source" or ROOT / "source" in path.parents):
            parser.error("Keep output and work directories outside the editable source/ folder.")
    try:
        with workspace(args.work) as work:
            sources = stage_inputs(work)
            jobs = {stem: compile_one(engine, stem, sources, work) for stem in DOCUMENTS}
            logs = out / "logs"
            logs.mkdir(parents=True, exist_ok=True)
            for stem, job in jobs.items():
                publish_file(job / (stem + ".pdf"), out / (stem + ".pdf"))
                for ending in (".log", ".fls", "-inputs.txt"):
                    publish_file(job / (stem + ending), logs / (stem + ending))
                print(out / (stem + ".pdf"))
    except (OSError, RuntimeError) as exc:
        parser.exit(1, str(exc) + "\n")


if __name__ == "__main__":
    main()
