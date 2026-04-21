#!/usr/bin/env python3

import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from xml.etree import ElementTree as ET

import yaml
from PIL import Image, UnidentifiedImageError

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
NS = {"svg": SVG_NS}

SCRIPT_DIR = Path(__file__).resolve().parent

DEFAULT_CONFIG_PATH = Path(".nf-core.yml")
DEFAULT_OUTPUT_DIR = Path("docs/images")
DEFAULT_TEMPLATE_PATH = SCRIPT_DIR / "sanger-tol.template.svg"
DEFAULT_HEIGHT = "4cm"
ACCENT_COLOR = "#597fba"
MODE_TEXT_COLORS = {
    "dark": "#b2c9d3",
    "light": "#232642",
}

ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)


def non_empty_string(value: str) -> str:
    if value is None:
        raise argparse.ArgumentTypeError("value must be a non-empty string")
    val = str(value).strip()
    if not val:
        raise argparse.ArgumentTypeError("value must be a non-empty string")
    return val


def length_arg_type(value: str) -> str:
    try:
        parse_number_with_unit(value)
    except ValueError as exc:  # pragma: no cover - validation
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render Sanger ToL logos from an SVG template by replacing the "
            "placeholder text, widening the SVG when the new name needs more room, "
            "and exporting SVG and PNG outputs for dark and light modes."
        )
    )
    parser.add_argument(
        "name",
        nargs="?",
        type=non_empty_string,
        help="Replacement text for the placeholder label. Defaults to template.name from .nf-core.yml.",
    )
    parser.add_argument(
        "--org",
        type=non_empty_string,
        default=None,
        help="Organisation used for output filenames. Defaults to template.org from .nf-core.yml.",
    )
    parser.add_argument(
        "--pipeline-dir",
        type=Path,
        default=Path("."),
        help=(
            "Pipeline directory used to resolve defaults for --config and --output-dir. "
            "Defaults to the current directory."
        ),
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="Path to the .nf-core.yml file used to discover template.org and template.name.",
    )
    parser.add_argument(
        "--template",
        type=Path,
        default=DEFAULT_TEMPLATE_PATH,
        help="Path to the input SVG template used as the geometry source.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Directory for the generated SVG and PNG files. If not provided, "
            "this is resolved relative to --pipeline-dir (defaults to docs/images)."
        ),
    )
    parser.add_argument(
        "--mode",
        choices=["dark", "light", "all"],
        default="all",
        help="Which color mode to generate. 'all' writes both dark and light variants.",
    )
    parser.add_argument(
        "--height",
        type=length_arg_type,
        default=DEFAULT_HEIGHT,
        help=(
            "Output logo height to write into the SVG (e.g. '4cm', '120px'). "
            "Width will be computed to preserve aspect ratio."
        ),
    )
    parser.add_argument(
        "--inkscape",
        default="inkscape",
        help=("Path or command name for the Inkscape executable. Defaults to 'inkscape', which is looked up on $PATH."),
    )

    args = parser.parse_args()

    if not args.inkscape or not str(args.inkscape).strip():
        parser.error("--inkscape must be a non-empty string")
    if shutil.which(str(args.inkscape)) is None:
        parser.error(f"Inkscape executable not found in $PATH: {args.inkscape}")
    if not args.template.exists():
        parser.error(f"Template not found: {args.template}")

    # Store the resolved paths back on the args namespace for downstream use
    args.config = args.pipeline_dir / args.config
    args.output_dir = args.pipeline_dir / args.output_dir
    if not args.config.exists() and (args.name is None or args.org is None):
        parser.error(f"Metadata file not found and overrides are incomplete: {args.config}")

    return args


def parse_number_with_unit(value: str) -> tuple[float, str]:
    match = re.fullmatch(r"\s*([0-9]*\.?[0-9]+)([A-Za-z%]*)\s*", value)
    if match is None:
        raise ValueError(f"Unsupported SVG length: {value!r}")
    return float(match.group(1)), match.group(2)


def parse_scale(transform: str | None) -> float:
    if not transform:
        return 1.0
    match = re.fullmatch(r"\s*scale\(([-+]?[0-9]*\.?[0-9]+)\)\s*", transform)
    if match is None:
        raise ValueError(f"Unsupported group transform: {transform!r}")
    return float(match.group(1))


def read_project_metadata(config_path: Path) -> tuple[str, str]:
    with config_path.open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}

    template = payload.get("template")
    if not isinstance(template, dict):
        raise ValueError(f"Missing template section in {config_path}")

    org = template.get("org")
    name = template.get("name")
    if not isinstance(org, str) or not org.strip():
        raise ValueError(f"Missing template.org in {config_path}")
    if not isinstance(name, str) or not name.strip():
        raise ValueError(f"Missing template.name in {config_path}")

    return org.strip(), name.strip()


def resolve_project_metadata(
    config_path: Path,
    name_override: str | None,
    org_override: str | None,
) -> tuple[str, str]:
    config_org: str | None = None
    config_name: str | None = None
    if config_path.exists():
        config_org, config_name = read_project_metadata(config_path)

    org = org_override or config_org
    name = name_override or config_name
    if org is None or name is None:
        raise ValueError("Both org and name must be provided, either via .nf-core.yml or CLI overrides")
    return org, name


def resolve_modes(mode: str) -> list[str]:
    if mode == "all":
        return ["dark", "light"]
    return [mode]


def measure_text_width(
    text: str,
    font_size: float,
    font_family: str,
    font_weight: str,
    inkscape_cmd: str,
) -> float:
    probe_padding = int(max(4, round(font_size * 2.0)))
    probe_baseline = int(max(4, round(font_size * 2.0)))
    probe_width = int(max(probe_padding * 2 + 16, round(font_size * max(len(text), 1) * 2.2)))
    probe_height = int(max(probe_padding * 2 + 16, round(font_size * 4.0)))

    probe_svg = (
        f'<svg xmlns="{SVG_NS}" width="{probe_width}" height="{probe_height}" '
        f'viewBox="0 0 {probe_width} {probe_height}">'
        f'<text x="{probe_padding}" y="{probe_baseline}" '
        f'font-family="{font_family}" font-size="{font_size}" font-weight="{font_weight}" '
        f'fill="#ffffff">{escape_xml(text)}</text></svg>'
    )

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            svgtmp_path = Path(tmpdir) / "probe.svg"
            pngtmp_path = Path(tmpdir) / "probe.png"
            svgtmp_path.write_bytes(probe_svg.encode("utf-8"))
            export_svg_to_png(inkscape_cmd, svgtmp_path, pngtmp_path)

            image = Image.open(pngtmp_path).convert("RGBA")
            alpha = image.split()[3]
            bbox = alpha.getbbox()
            if bbox is None:
                raise ValueError("Rendered text did not produce a bounding box")
            left, _, right, _ = bbox
            return float(right - left)
    except (
        OSError,
        UnidentifiedImageError,
        ValueError,
        FileNotFoundError,
        subprocess.CalledProcessError,
    ):
        return font_size * 0.6 * max(len(text), 1)


def escape_xml(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def set_svg_dimensions(root: ET.Element, view_width: float, view_height: float, img_height: str | None) -> None:
    if img_height:
        height_value, unit = parse_number_with_unit(img_height)
        width_value = height_value * (view_width / view_height)
        root.attrib["width"] = format_length(width_value, unit)
        root.attrib["height"] = format_length(height_value, unit)
        return

    root.attrib["width"] = format_float(view_width)
    root.attrib["height"] = format_float(view_height)


def apply_mode_style(tree: ET.ElementTree, mode: str, name: str) -> None:
    root = tree.getroot()
    root.attrib["aria-label"] = f"Logo for {name} ({mode})"

    group = root.find("svg:g", NS)
    if group is None:
        raise ValueError("Template SVG does not contain the expected <g> element")

    text_nodes = group.findall("svg:text", NS)
    if len(text_nodes) < 2:
        raise ValueError("Template SVG does not contain the expected text nodes")

    brand_node = text_nodes[0]
    name_node = text_nodes[1]
    tspans = brand_node.findall(".//svg:tspan", NS)
    if len(tspans) < 2:
        raise ValueError("Template SVG does not contain the expected brand tspans")

    tspans[0].attrib["fill"] = ACCENT_COLOR
    tspans[1].attrib["fill"] = MODE_TEXT_COLORS[mode]
    name_node.attrib["fill"] = MODE_TEXT_COLORS[mode]


def update_logo_text(
    tree: ET.ElementTree[ET.Element],
    name: str,
    inkscape_cmd: str,
    mode: str,
    img_height: str | None = None,
) -> None:
    apply_mode_style(tree, mode, name)

    root = tree.getroot()
    group = root.find("svg:g", NS)
    if group is None:
        raise ValueError("Template SVG does not contain the expected <g> element")

    text_nodes = group.findall("svg:text", NS)
    if len(text_nodes) < 2:
        raise ValueError("Template SVG does not contain the expected text nodes")
    name_node = text_nodes[1]

    image_node = group.find("svg:image", NS)
    if image_node is None:
        raise ValueError("Template SVG does not contain the expected image node")

    original_name = "".join(name_node.itertext())
    font_size = float(name_node.attrib["font-size"])
    font_family = name_node.attrib["font-family"]
    font_weight = name_node.attrib["font-weight"]

    name_node.text = name

    original_width = measure_text_width(original_name, font_size, font_family, font_weight, inkscape_cmd)
    new_width = measure_text_width(name, font_size, font_family, font_weight, inkscape_cmd)
    text_x = float(name_node.attrib["x"])
    image_x = float(image_node.attrib["x"])
    image_width = float(image_node.attrib["width"])

    view_box = root.attrib.get("viewBox")
    if not view_box:
        raise ValueError("Template SVG is missing a viewBox")
    min_x, min_y, view_width, view_height = parse_view_box(view_box)

    scale = parse_scale(group.attrib.get("transform"))
    inner_canvas_width = view_width / scale
    original_content_right = max(image_x + image_width, text_x + original_width)
    right_padding = inner_canvas_width - original_content_right
    if right_padding < 0:
        raise ValueError("Template content exceeds the original canvas width")

    new_content_right = max(image_x + image_width, text_x + new_width)
    required_inner_width = new_content_right + right_padding
    if required_inner_width <= inner_canvas_width:
        set_svg_dimensions(root, view_width, view_height, img_height)
        return

    new_view_width = required_inner_width * scale
    root.attrib["viewBox"] = format_view_box(min_x, min_y, new_view_width, view_height)
    set_svg_dimensions(root, new_view_width, view_height, img_height)


def parse_view_box(value: str) -> tuple[float, float, float, float]:
    parts = value.replace(",", " ").split()
    if len(parts) != 4:
        raise ValueError(f"Invalid viewBox: {value!r}")
    return tuple(float(part) for part in parts)  # type: ignore[return-value]


def format_view_box(min_x: float, min_y: float, width: float, height: float) -> str:
    return f"{format_float(min_x)} {format_float(min_y)} {format_float(width)} {format_float(height)}"


def format_float(value: float) -> str:
    text = f"{value:.4f}"
    text = text.rstrip("0").rstrip(".")
    return text or "0"


def format_length(value: float, unit: str) -> str:
    return f"{format_float(value)}{unit}"


def export_svg_to_png(inkscape_cmd: str, input_path: Path, output_path: Path) -> None:
    """Run Inkscape to export an SVG file to a PNG file."""
    subprocess.check_call(
        [
            inkscape_cmd,
            str(input_path),
            "--export-type=png",
            "--export-filename",
            str(output_path),
        ]
    )


def build_logo(
    name: str,
    template_path: Path,
    svg_output_path: Path,
    png_output_path: Path,
    img_height: str,
    inkscape_cmd: str,
    mode: str,
) -> None:
    tree = ET.parse(template_path)
    update_logo_text(tree, name, inkscape_cmd, mode, img_height)
    ET.indent(tree, space=" ")

    svg_output_path.parent.mkdir(parents=True, exist_ok=True)
    png_output_path.parent.mkdir(parents=True, exist_ok=True)

    tree.write(svg_output_path, encoding="utf-8", xml_declaration=False)
    export_svg_to_png(inkscape_cmd, svg_output_path, png_output_path)


def build_logos(
    org: str,
    name: str,
    template_path: Path,
    output_dir: Path,
    img_height: str,
    inkscape_cmd: str,
    modes: list[str],
) -> list[Path]:
    output_paths: list[Path] = []
    base_name = f"{org}-{name}_logo"
    for mode in modes:
        stem = f"{base_name}_{mode}"
        svg_output_path = output_dir / f"{stem}.svg"
        png_output_path = output_dir / f"{stem}.png"
        build_logo(
            name,
            template_path,
            svg_output_path,
            png_output_path,
            img_height,
            inkscape_cmd,
            mode,
        )
        output_paths.extend([svg_output_path, png_output_path])
    return output_paths


def main() -> int:
    args = parse_args()
    org, name = resolve_project_metadata(args.config, args.name, args.org)
    output_paths = build_logos(
        org=org,
        name=name,
        template_path=args.template,
        output_dir=args.output_dir,
        img_height=args.height,
        inkscape_cmd=args.inkscape,
        modes=resolve_modes(args.mode),
    )

    for path in output_paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
