"""Extract multiple DXA meshes and build an offline interactive comparison page."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from html import escape
import math
from pathlib import Path, PurePosixPath
import sys
from urllib.parse import quote, urlsplit

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))


def _relative_link(output: Path, value: object) -> str:
    """Validate a local artifact path and encode it for an HTML attribute."""
    raw = str(value).replace("\\", "/")
    relative = PurePosixPath(raw)
    if urlsplit(raw).scheme or raw.startswith("//") or relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"comparison artifacts must have paths relative to the output directory: {raw}")
    target = (output / Path(relative)).resolve()
    if not target.is_relative_to(output.resolve()) or not target.is_file():
        raise ValueError(f"comparison artifact does not exist inside the output directory: {raw}")
    return escape(quote(relative.as_posix(), safe="/"), quote=True)


def build_index(output_dir: str | Path, records: Sequence[Mapping[str, object]]) -> Path:
    """Build ``index.html`` from existing per-structure views without running DXA.

    Each record supplies ``name``, ``native_faces``, ``display_faces``,
    ``cell_lengths`` (three values), and root-relative ``html`` and ``png``
    paths. Optional ``surface_smoothing`` defaults to zero; optional
    ``length_unit`` defaults to ``input_length_unit``. All referenced files
    must exist below ``output_dir``. The page has no network dependencies.
    """
    if not records:
        raise ValueError("at least one result record is required")
    output = Path(output_dir)
    cards = []
    for record in records:
        name = escape(str(record["name"]), quote=True)
        native = int(record["native_faces"])
        displayed = int(record["display_faces"])
        lengths = tuple(float(value) for value in record["cell_lengths"])
        smoothing = int(record.get("surface_smoothing", 0))
        if native < 0 or displayed < 0 or smoothing < 0:
            raise ValueError("triangle counts and smoothing iterations must be nonnegative")
        if len(lengths) != 3 or any(not math.isfinite(value) or value <= 0 for value in lengths):
            raise ValueError("cell_lengths must contain three finite positive lengths")
        unit = escape(str(record.get("length_unit", "input_length_unit")), quote=True)
        dimensions = " × ".join(f"{value:.6g}" for value in lengths)
        html_link = _relative_link(output, record["html"])
        png_link = _relative_link(output, record["png"])
        cards.append(f"""<article class="card">
<header><h2>{name}</h2><nav><a href="{html_link}" target="_blank" rel="noopener">单独打开交互视图</a><a href="{png_link}" target="_blank" rel="noopener">PNG</a></nav></header>
<dl><div><dt>原始周期三角形</dt><dd>{native:,}</dd></div><div><dt>显示三角形</dt><dd>{displayed:,}</dd></div><div><dt>盒矢量长度</dt><dd>{dimensions} {unit}</dd></div><div><dt>平滑次数</dt><dd>{smoothing}</dd></div></dl>
<iframe src="{html_link}" title="{name} DXA defect mesh" loading="lazy" allowfullscreen></iframe>
</article>""")
    page = """<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>DXA defect mesh 对比</title>
<style>
*{box-sizing:border-box}body{margin:0;background:#eef2f6;color:#172333;font:15px/1.5 system-ui,sans-serif}
main{max-width:1600px;margin:auto;padding:24px}h1{margin:0 0 6px;font-size:25px}main>p{margin:0 0 20px;color:#526174}
.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}.card{background:white;border:1px solid #dce3ec;border-radius:12px;overflow:hidden}
header{padding:14px 18px 8px;display:flex;align-items:center;justify-content:space-between;gap:12px;flex-wrap:wrap}h2{font-size:18px;margin:0}
nav{display:flex;gap:14px;font-size:13px}a{color:#175f9f;text-decoration:none}a:hover{text-decoration:underline}
dl{margin:0;padding:0 18px 12px;display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:6px 14px;font-size:13px}
dl>div{min-width:0}dt{color:#69788a}dd{margin:0;overflow-wrap:anywhere;font-variant-numeric:tabular-nums}
iframe{display:block;width:100%;height:420px;border:0;background:#fff}
@media(max-width:800px){main{padding:14px}.grid{grid-template-columns:1fr}iframe{height:360px}}
</style></head><body><main><h1>DXA defect mesh 对比</h1>
<p>在每个视图中拖动旋转、滚轮缩放；可单独打开交互视图或查看 PNG。显示网格保留周期边界截断，不添加封口。</p>
<section class="grid">""" + "\n".join(cards) + "\n</section></main></body></html>\n"
    output.mkdir(parents=True, exist_ok=True)
    destination = output / "index.html"
    destination.write_text(page, encoding="utf-8")
    return destination


def extract_dxa_defect_mesh_batch(
    filenames: Sequence[str | Path],
    *,
    output_dir: str | Path,
    surface_smoothing: int = 0,
    crystal_structure: str = "fcc",
    length_unit: str = "input_length_unit",
    show_atoms: bool = False,
) -> tuple[list[dict[str, object]], Path]:
    """Extract each input into its stem directory and return records/index path.

    All inputs and case-insensitive stem collisions are checked before any
    extraction. An extraction or export failure stops the batch immediately.
    Each case receives a PNG and a self-contained interactive HTML view.
    """
    files = [Path(filename).resolve() for filename in filenames]
    if not files:
        raise ValueError("at least one input structure is required")
    missing = [str(filename) for filename in files if not filename.is_file()]
    if missing:
        raise FileNotFoundError("input structures do not exist: " + ", ".join(missing))
    stems = [filename.stem.casefold() for filename in files]
    if len(stems) != len(set(stems)):
        raise ValueError("input files must have distinct stems to avoid output-directory collisions")

    from example.crystalops.analysis.extract_dxa_defect_mesh import extract_dxa_defect_mesh_file

    output = Path(output_dir).resolve()
    records = []
    for filename in files:
        case_output = output / filename.stem
        result, paths = extract_dxa_defect_mesh_file(
            filename, output_dir=case_output, crystal_structure=crystal_structure,
            surface_smoothing=surface_smoothing, length_unit=length_unit,
            plot=True, show_atoms=show_atoms, interactive=True,
        )
        record = {
            "name": filename.name,
            "native_faces": len(result.faces),
            "display_faces": len(result.clipped_faces),
            "cell_lengths": [math.sqrt(sum(float(value) ** 2 for value in vector)) for vector in result.cell],
            "surface_smoothing": surface_smoothing,
            "length_unit": length_unit,
            "html": Path(paths["html"]).resolve().relative_to(output).as_posix(),
            "png": Path(paths["png"]).resolve().relative_to(output).as_posix(),
        }
        records.append(record)
        print(f"{filename.name}: {record['native_faces']} periodic / {record['display_faces']} display triangles")
    index = build_index(output, records)
    return records, index


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="complete atomistic structure files")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--surface-smoothing", type=int, default=0)
    parser.add_argument(
        "--crystal-structure", default="fcc",
        choices=("fcc", "hcp", "bcc", "cubic_diamond", "hexagonal_diamond"),
        help="reference atomic lattice for DXA",
    )
    parser.add_argument("--length-unit", default="input_length_unit", help="coordinate-unit label only, no conversion")
    parser.add_argument("--show-atoms", action="store_true", help="overlay the original atom positions in every view")
    options = parser.parse_args()
    _, index = extract_dxa_defect_mesh_batch(
        options.inputs, output_dir=options.output_dir,
        surface_smoothing=options.surface_smoothing, crystal_structure=options.crystal_structure,
        length_unit=options.length_unit, show_atoms=options.show_atoms,
    )
    print(f"Comparison page: {index.resolve()}")


if __name__ == "__main__":
    main()
