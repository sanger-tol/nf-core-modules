#!/usr/bin/env python3

"""
Read the pipeline manifest and generate the CITATION.cff file
"""

import datetime
import logging
import operator
import os
from pathlib import Path

import rich_click as click

from nf_core.pipelines.lint_utils import dump_yaml_with_prettier
from nf_core.pipelines.rocrate import ROCrate
from nf_core.utils import Pipeline

log = logging.getLogger(__name__)


# Only update the dictionary if there's a value
def set_if_set(d, k, v):
    if v is not None:
        sv = v.strip()
        if sv:
            d[k] = sv


## Hijack nf-core's ROCrate code to parse manifest.contributors


# Object with a dummy "append_to" method
class DummyCrateEntity:
    def append_to(self, mode, author_entity):
        pass


# Object with a dummy "resolve_id" method, and an "add" method that records all values
class DummyCrate:
    def __init__(self):
        self.elements = []

    def add(self, obj):
        self.elements.append(obj)

    def resolve_id(self, _id):
        return _id


def get_contributors(pipeline_dir):
    # The nf-core ROCrate object
    rocrate_obj = ROCrate(pipeline_dir)
    # "add_main_authors" needs "self.crate" with "add" and "resolve_id"
    rocrate_obj.crate = DummyCrate()
    # "add_main_authors" needs "wf_file" with "append_to"
    dummy_wf_file = DummyCrateEntity()
    rocrate_obj.add_main_authors(dummy_wf_file)
    return rocrate_obj.crate.elements


message = (
    "If you use this software, please cite it using the metadata from this file and all references from CITATIONS.md."
)


def get_pipeline(path):
    pipeline_obj = Pipeline(path)
    pipeline_obj._load()
    return pipeline_obj


# The release name isn't in the manifest. Look for it in the CHANGELOG
def find_release_name(pipeline_dir, version):
    changelog = pipeline_dir / "CHANGELOG.md"
    with changelog.open() as f:
        for line in f:
            if (
                line.startswith("#")
                and version in line
                and ("2024" in line or "2025" in line or "2026" in line or "2027" in line)
            ):
                # We have a variety of line structures (and hyphen styles !)
                ## [[0.7.1](https://github.com/sanger-tol/blobtoolkit/releases/tag/0.7.1)] – Psyduck (patch 1) – [2025-03-29]
                ## [[0.7.0](https://github.com/sanger-tol/blobtoolkit/releases/tag/0.7.0)] – Psyduck – [2025-03-19]
                ## [1.2.2] - Ancient Destiny (H2)- [2025-01-30]
                ## [1.2.0] - Ancient Destiny - [2024-11-15]
                ## [[1.3.1](https://github.com/sanger-tol/curationpretext/releases/tag/1.3.1)] - UNSC Pillar-of-Autumn (H1) - [2025-04-02]
                ## [[1.3.0](https://github.com/sanger-tol/curationpretext/releases/tag/1.3.0)] - UNSC Pillar-of-Autumn - [2025-02-27]
                ## v0.7.0 - Raymond Carhart [08/03/2025]
                ## [1.1.0dev] - unnamed - [2025-03-31]
                ## [[1.0.1](https://github.com/sanger-tol/metagenomeassembly/releases/tag/1.0.1)] - Scarborough Fair (patch 1) - [2025-03-31]
                line = line.replace("–", "-")
                name = line.split("- ")[1].strip()
                return name


# Intro: https://citation-file-format.github.io/
# Schema: https://github.com/citation-file-format/citation-file-format/blob/main/schema-guide.md
# Validate with `cffconvert`
def build_cff(pipeline_dir):

    pipeline_obj = get_pipeline(pipeline_dir)
    manifest = pipeline_obj.nf_config.get("manifest")
    if not manifest:
        log.error("No manifest in nextflow.config")
        return

    pipeline_name = manifest["name"]
    pipeline_version = manifest["version"]
    release_name = find_release_name(pipeline_obj.wf_path, pipeline_version)

    content = {
        "cff-version": "1.2.0",
        "message": message,
        "type": "software",  # it's either that or "dataset"
        "url": f"https://pipelines.tol.sanger.ac.uk/{pipeline_name.split('/')[1]}",
        "license": "MIT",
        "title": f"{pipeline_name} v{pipeline_version}{' - ' + release_name if release_name else ''}",
        "version": pipeline_version,
        "date-released": datetime.date.today().isoformat(),
    }
    set_if_set(content, "repository-code", manifest.get("homePage"))
    set_if_set(content, "doi", manifest.get("doi"))

    contributors = get_contributors(pipeline_dir)
    log.debug("Parsed contributors: %s", contributors)
    if not contributors:
        log.error("Empty list of contributors in manifest of nextflow.config")
        return
    log.info(f"Found {len(contributors)} contributors")

    authors = []
    for contributor in contributors:
        # Note: contributor is a ROCrate "Person" object, which has a dict interface
        log.debug(f"Adding author: {contributor}")
        author = {}
        set_if_set(author, "affiliation", contributor.get("affiliation"))
        set_if_set(author, "orcid", contributor.get("@id"))
        set_if_set(author, "email", contributor.get("email"))
        set_if_set(author, "website", contributor.get("url"))
        if "," in contributor["name"]:
            # "Family name, Given name"
            (family, given) = contributor["name"].split(",", 1)
        elif " " in contributor["name"]:
            # "Given Family" (no whitespace allowed in the given name)
            (given, family) = contributor["name"].split(maxsplit=1)
        else:
            # "Given"
            given = contributor["name"]
            family = None
        set_if_set(author, "given-names", given)
        set_if_set(author, "family-names", family)
        authors.append(author)
    # Authors ordered by alphabetical order
    authors.sort(key=operator.itemgetter("family-names"))
    content["authors"] = authors
    return content


@click.command()
@click.argument(
    "pipeline_dir",
    type=click.Path(exists=True),
    default=Path.cwd(),
    required=True,
    metavar="<pipeline directory>",
)
def cff(pipeline_dir):
    pipeline_dir = Path(pipeline_dir)
    content = build_cff(pipeline_dir)
    # dump_yaml_with_prettier expects to be run from the directory of the repository
    os.chdir(pipeline_dir)
    # Use the function from nf-core to make sure the CFF file is pretty
    dump_yaml_with_prettier("CITATION.cff", content)


if __name__ == "__main__":
    cff()
