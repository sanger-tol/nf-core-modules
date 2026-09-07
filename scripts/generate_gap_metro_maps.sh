#!/bin/bash
MMD="$1"
NAME="${MMD//.mmd}"
LOGO="${NAME//metro_map/logo}"
render () {
  nf-metro render "${MMD}" -o "${NAME}_$1.svg" --theme "$2" --logo "${LOGO}_$1.png" --no-chrome-css
  cairosvg "${NAME}_$1.svg" -o "${NAME}_$1.png"
  nf-metro render "${MMD}" -o "${NAME}_$1.svg" --theme "$2" --logo "${LOGO}_$1.png"
}
render dark nfcore
render light light
