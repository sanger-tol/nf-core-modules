#!/usr/bin/env python3

"""
Read the pipeline manifest and generate the ro-crate-metadata.json file
"""

import logging
from pathlib import Path

import requests
import rich_click as click
import rocrate.rocrate

from nf_core.pipelines.rocrate import ROCrate

log = logging.getLogger(__name__)


class SangerToLROCrate(ROCrate):
    """
    Class to generate an RO Crate for a pipeline
    Overrides and complements the code written by nf-core
    """

    def make_workflow_rocrate(self) -> None:
        super().make_workflow_rocrate()

        # Link to the pipelines website instead of the nf-core website
        self.crate.mainEntity["url"] = [
            self.crate.mainEntity["url"][0],
            self.crate.mainEntity["url"][1].replace(
                "https://nf-co.re/sanger-tol", "https://pipelines.tol.sanger.ac.uk"
            ),
        ]

        # Change the Organization object and relink sdPublisher
        self.crate.delete("https://nf-co.re/")
        self.crate.add_jsonld(
            {
                "@id": "https://pipelines.tol.sanger.ac.uk/",
                "@type": "Organization",
                "name": "Sanger Tree of Life programme",
                "url": "https://pipelines.tol.sanger.ac.uk/",
            }
        )
        self.crate.mainEntity["sdPublisher"] = {"@id": "https://pipelines.tol.sanger.ac.uk/"}

        # Same code as upstream, but fetching keywords from the pipelines website
        remote_workflows = requests.get("https://pipelines.tol.sanger.ac.uk/pipelines.json").json()["remote_workflows"]
        # Also put "nextflow" first in the list
        topics = ["nextflow", "nf-core"]
        for remote_wf in remote_workflows:
            assert self.pipeline_obj.pipeline_name is not None  # mypy
            if remote_wf["name"] == self.pipeline_obj.pipeline_name:
                topics = topics + remote_wf["topics"]
                break

        log.debug(f"Adding topics: {topics}")
        self.crate.mainEntity["keywords"] = topics


@click.command()
@click.argument(
    "pipeline_dir",
    type=click.Path(exists=True),
    default=Path.cwd(),
    required=True,
    metavar="<pipeline directory>",
)
def rocrate(pipeline_dir):
    pipeline_dir = Path(pipeline_dir)
    rocrate_obj = SangerToLROCrate(pipeline_dir)
    rocrate_obj.create_rocrate(json_path=pipeline_dir)


if __name__ == "__main__":
    rocrate()
