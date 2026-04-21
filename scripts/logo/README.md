# nf-logo

This repo contains a small utility to render Sanger ToL logos from an SVG
template by replacing a placeholder label and exporting SVG and PNG assets
for dark and light backgrounds.

## Simple usage

Run the script from a pipeline directory that contains `.nf-core.yml` and it
will read `template.org` and `template.name`, then write all four default
outputs into `docs/images/`:

```bash
python /path/to/nf-logo/render_logo.py
```

That writes:

- `docs/images/<org>-<name>_logo_dark.svg`
- `docs/images/<org>-<name>_logo_dark.png`
- `docs/images/<org>-<name>_logo_light.svg`
- `docs/images/<org>-<name>_logo_light.png`

You can override the discovered metadata, choose a single mode, change the
output directory, or point at a different `.nf-core.yml`:

```bash
python render_logo.py custom-name \
  --org sanger-tol \
  --config ../sequencecomposition/.nf-core.yml \
  --mode light \
  --output-dir docs/images \
  --height 4cm \
  --inkscape /usr/bin/inkscape
```

## Install

```bash
pip install -e .
```

## Notes

- The script is `render_logo.py`.
- The script reads `template.org` and `template.name` from `.nf-core.yml` by default.
- The bundled geometry source is `sanger-tol-dark.template.svg`; the script applies dark and light text colors itself.
- The script will increase the SVG canvas width only if the replacement
  text requires more horizontal space to avoid colliding with the
  embedded image.
- PNG export is performed by calling Inkscape (`--inkscape`); Inkscape
  must be available on your `PATH` (or supply an explicit path).
- Text-width measurement probes a rendered glyph bounding box using
  Inkscape when possible; if rendering or Inkscape is unavailable the
  script falls back to a heuristic.
- The CLI exposes `--height` (e.g. `4cm`) for the output SVG.
- On success the script prints the paths of the generated files.
