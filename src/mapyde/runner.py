"""
File containing functionality for running the various steps in the workflow.
"""
from __future__ import annotations

import sys
import typing as T
from pathlib import Path

from mapyde.backends import madgraph
from mapyde.container import Container


def run_madgraph(config: dict[str, T.Any]) -> tuple[bytes, bytes]:
    """
    Run madgraph.
    """
    # ./test/wrapper_mgpy.py config_file
    madgraph.generate_mg5config(config)

    image = f"ghcr.io/scipp-atlas/mario-mapyde/{config['madgraph']['version']}"
    command = b"mg5_aMC /data/run.mg5 && rsync -a PROC_madgraph /data/madgraph\n"

    with Container(
        image=image,
        name=f"{config['base']['output']}__mgpy",
        mounts=[
            (str(Path(config["base"]["path"]).joinpath("cards").resolve()), "/cards"),
            (
                str(
                    Path(config["base"]["path"])
                    .joinpath(config["base"]["output"])
                    .resolve()
                ),
                "/data",
            ),
        ],
        stdout=sys.stdout,
    ) as container:
        stdout, stderr = container.process.communicate(command)

    return stdout, stderr


def run_delphes(config: dict[str, T.Any]) -> tuple[bytes, bytes]:
    """
    Run delphes.
    """
    # ./test/wrapper_delphes.py config_file
    image = f"ghcr.io/scipp-atlas/mario-mapyde/{config['delphes']['version']}"
    command = bytes(
        f"""pwd && ls -lavh && ls -lavh /data && cp $(find /data/ -name "*hepmc.gz") hepmc.gz && \
gunzip hepmc.gz && \
/bin/ls -ltrh --color && \
/usr/local/share/delphes/delphes/DelphesHepMC2 /cards/delphes/{config['delphes']['card']} delphes.root hepmc && \
rsync -rav --exclude hepmc . /data/delphes""",
        "utf-8",
    )

    with Container(
        image=image,
        name=f"{config['base']['output']}__delphes",
        mounts=[
            (str(Path(config["base"]["path"]).joinpath("cards").resolve()), "/cards"),
            (
                str(
                    Path(config["base"]["path"])
                    .joinpath(config["base"]["output"])
                    .resolve()
                ),
                "/data",
            ),
        ],
        stdout=sys.stdout,
    ) as container:
        stdout, stderr = container.process.communicate(command)

    return stdout, stderr


def run_ana(config: dict[str, T.Any]) -> tuple[bytes, bytes]:
    """
    Run analysis.
    """
    # ./test/wrapper_delphes.py config_file
    # recompute XS, override XS if needed,
    # add XS to config, re-dump config file
    # probably want to move this to wrapper_ana.py script

    xsec = 10  # TODO

    image = f"ghcr.io/scipp-atlas/mario-mapyde/{config['delphes']['version']}"
    command = bytes(
        f"""/scripts/{config['analysis']['script']} --input /data/delphes/delphes.root --output {config['analysis']['output']} --lumi {config['analysis']['lumi']} --XS {xsec} && rsync -rav . /data/analysis""",
        "utf-8",
    )

    with Container(
        image=image,
        name=f"{config['base']['output']}__hists",
        mounts=[
            (str(Path(config["base"]["path"]).joinpath("cards").resolve()), "/cards"),
            (
                str(Path(config["base"]["path"]).joinpath("scripts").resolve()),
                "/scripts",
            ),
            (
                str(
                    Path(config["base"]["path"])
                    .joinpath(config["base"]["output"])
                    .resolve()
                ),
                "/data",
            ),
        ],
        stdout=sys.stdout,
    ) as container:
        stdout, stderr = container.process.communicate(command)

    return stdout, stderr


def run_pyhf(config: dict[str, T.Any]) -> None:
    """
    Run statistical inference via pyhf.
    """
    # ./test/wrapper_pyhf.py config_file
    assert config
