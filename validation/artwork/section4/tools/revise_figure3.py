"""Revise only the panel-(a) annotations in the original vector Figure 3.

Usage: python tools/revise_figure3.py baseline/paper1_figure3_v0.pdf OUTPUT.pdf
The project's Figure 4 supplies the matching embedded DejaVu Sans subset for
the new legend wording (the Figure 3 font subset lacks the regular letter q).
All dispersion curves, shading, axes, marked levels, and panel (b) are copied
as original PDF drawing instructions.
"""
import argparse
import hashlib
from pathlib import Path
import re

from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject, NameObject

ROOT = Path(__file__).resolve().parents[1]
SOURCE_SHA256 = '9882955f2a34701f5154841c2325fff32d3344f619e7a1438d2a8672c56575b8'
FONT_FIGURE = ROOT / 'context/figures/paper1_sectionVI_fig4_reference_partition_v2.pdf'
FONT_FIGURE_SHA256 = '7c5313c544de0373b3e234355c143d404ec3c569efcd946d7e477b47798f67c1'


def checked_pdf(path, sha256):
    if hashlib.sha256(path.read_bytes()).hexdigest() != sha256:
        raise ValueError(f'Input PDF differs from the reviewed source: {path}')
    return PdfReader(path)


def text_lines(*runs):
    """Place same-family text runs, including math subscripts, in page points."""
    commands = ['q', '0 g', 'BT']
    for text, font, size, x, y in runs:
        encoded = text.encode('utf-16be').hex()
        commands.extend([f'/{font} {size} Tf', f'1 0 0 1 {x} {y} Tm', f'<{encoded}> Tj'])
    commands.extend(['ET', 'Q'])
    return '\n'.join(commands).encode('ascii')


def replace_label(stream, old_origin, new_text):
    # Each annotation is an independent text group in the Matplotlib source.
    prefix = ('q\n1 0 -0 1 ' + old_origin + ' cm\nBT\n').encode('ascii')
    result, count = re.subn(re.escape(prefix) + rb'.*?\nET\nQ', lambda _: new_text,
                            stream, flags=re.DOTALL)
    if count != 1:
        raise ValueError(f'Expected one text group at {old_origin}; found {count}')
    return result


def revise(source, target):
    original = checked_pdf(source, SOURCE_SHA256)
    font_figure = checked_pdf(FONT_FIGURE, FONT_FIGURE_SHA256)
    original_stream = original.pages[0]['/Contents'].get_data()
    stream = original_stream
    # All coordinates refer to the original PDF page. Keep the label between
    # the drive and upper-minimum lines, horizontally centred on the same spot.
    stream = replace_label(stream, '89.1520535714 79.3738127964', text_lines(
        ('upper minimum at fixed ', 'F1', 9, 83.54, 79.37381),
        ('k', 'F3', 9, 195.35, 79.37381), ('x', 'F3', 6.3, 200.56, 77.89725)))
    stream = replace_label(stream, '124.3395535714 69.2863127964', text_lines(
        ('1.60', 'F1', 9, 133.82, 69.28631)))
    stream = replace_label(stream, '39.768125 246.55636', text_lines(
        ('fixed-', 'F1', 9, 39.768125, 246.71), ('k', 'F3', 9, 64.95, 246.71),
        ('x', 'F3', 6.3, 70.16, 245.23), (' slice', 'F1', 9, 74.05, 246.71)))
    stream = replace_label(stream, '39.768125 236.15636', text_lines(
        ('finite contact permits other ', 'F1', 9, 39.768125, 236.31),
        ('k', 'F3', 9, 165.59, 236.31), ('x', 'F3', 6.3, 170.80, 234.83)))

    # Replace the grey legend swatch with a stroke matching the blue data
    # boundaries. Adjust both handles to leave room for the wider label.
    legend_start = b'/A2 gs 0.82 g 0 j 0.82 G 0.82 g'
    legend_end = b'q\n1 0 -0 1 285.468125 '
    start = stream.index(legend_start)
    end = stream.index(legend_end, start)
    if stream.count(legend_start) != 1 or stream.count(legend_end) != 1:
        raise ValueError('Figure legend boundary is ambiguous')
    legend = (b'q\n/A3 gs 2 J 1.8 w 0.1215686275 0.4666666667 0.7058823529 RG\n'
              b'135 218.46 m 149 218.46 l S\n1 0.4980392157 0.0549019608 RG\n'
              b'135 207.24 m 149 207.24 l S\nQ\n')
    legend += text_lines(('two lower frequencies', 'FLegend', 8.5, 153.5, 215.31),
                         ('upper branch', 'F1', 8.5, 153.5, 204.09)) + b'\n'
    stream = stream[:start] + legend + stream[end:]

    # Everything preceding the first annotation, and the complete right-hand
    # panel after the legend, must be byte-identical in the uncompressed stream.
    plot_end = original_stream.index(b'q\n1 0 -0 1 89.1520535714 ')
    if stream[:plot_end] != original_stream[:plot_end] or (
        stream[stream.index(legend_end):] != original_stream[original_stream.index(legend_end):]
    ):
        raise ValueError('A plotted path, axis, or right-panel drawing was altered')

    writer = PdfWriter(clone_from=original)
    page = writer.pages[0]
    extra_font = font_figure.pages[0]['/Resources']['/Font'].raw_get('/F1')
    page['/Resources']['/Font'][NameObject('/FLegend')] = extra_font.clone(writer)
    new_content = DecodedStreamObject()
    new_content.set_data(stream)
    page.replace_contents(new_content)
    page.compress_content_streams()
    writer.write(target)
    if PdfReader(target).pages[0]['/Contents'].get_data() != stream:
        raise ValueError('Written page differs from the verified vector content')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('original', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    revise(args.original, args.output)
