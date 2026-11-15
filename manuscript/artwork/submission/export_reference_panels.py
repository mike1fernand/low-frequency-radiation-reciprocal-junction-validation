"""Recreate Figure 4 and its SM4 companion from the archived vector drawing."""
import argparse
import hashlib
from pathlib import Path

import pymupdf


SOURCE_NAME = "paper1_sectionVI_fig4_reference_partition_v2.pdf"
SOURCE_SHA256 = "b0d17bcc337d7b2271998a29bc93094889114a932c318eca0f0531aa15f948fb"
PANELS = (
    ("a", 0, 250, "Outgoing power", "figure4_reference_power.pdf"),
    ("b", 255, 504, "qT conversion by method", "sm4_method_comparison.pdf"),
)


def export_panel(source, output, letter, left, right, title, filename):
    with pymupdf.open(source) as original, pymupdf.open() as result:
        page = original[0]
        viewport = pymupdf.Rect(left, 0, right, page.rect.height)
        # The original embedded font stores these characters as big-endian UTF-16.
        # In a PDF literal string, the parentheses around the letter are escaped.
        prefix = b"\x00\\(\x00" + letter.encode("ascii") + b"\x00\\)\x00 \x00 "
        heading = title.encode("utf-16be")
        content = page.read_contents()
        if content.count(prefix + heading) != 1:
            raise ValueError(f"Cannot uniquely locate the title of panel {letter}.")
        stream, = page.get_contents()
        original.update_stream(stream, content.replace(prefix + heading, heading))
        panel = result.new_page(width=viewport.width, height=viewport.height)
        panel.show_pdf_page(panel.rect, original, 0, clip=viewport)
        result.save(output / filename, garbage=4, deflate=True, no_new_id=True)


def main():
    figures = Path(__file__).resolve().parents[2] / "source" / "figures"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=figures / SOURCE_NAME)
    parser.add_argument("--output-dir", type=Path, default=figures)
    args = parser.parse_args()
    if hashlib.sha256(args.source.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("This recipe is for the original two-panel reference drawing.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for settings in PANELS:
        export_panel(args.source, args.output_dir, *settings)


if __name__ == "__main__":
    main()
