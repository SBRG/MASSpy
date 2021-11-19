# MASSpy - Mass Action Stoichiometric Simulation in Python

[![MASSpy release (latest SemVer)][1]](https://github.com/SBRG/MASSpy/releases)
[![PyPI - Python Version][2]](https://pypi.org/project/masspy/)
[![GitHub Workflow Status - Tests][3]](https://github.com/SBRG/MASSpy/actions/workflows/main.yml)
[![Read the Docs][4]](https://masspy.readthedocs.io/)
[![Codecov][5]](https://app.codecov.io/gh/SBRG/MASSpy)
[![GitHub License][6]](https://github.com/SBRG/MASSpy/blob/main/LICENSE)
[![Docker Image Size (latest semver)][7]](https://hub.docker.com/r/sbrg/masspy)
[![Code style: black][8]](https://github.com/psf/black)
[![pre-commit][9]](https://github.com/pre-commit/pre-commit)
[![MASSpy Publication][10]](https://doi.org/10.1371/journal.pcbi.1008208)

<p align="center"><img src="https://raw.githubusercontent.com/SBRG/MASSpy/main/docs/images/masspy-logo.svg" alt="MASSpy-Symbol" width="250" height="250"/></p>

## What is MASSpy?

The **M**ass **A**ction **S**toichiometric **S**imulation **py**thon
(**MASSpy**) package contains modules for the construction, simulation, and
analysis of kinetic models of biochemical reaction systems.

**MASSpy** is built to integrate seemlessly with [**COBRApy**][11], a widely used
modeling software package for constraint-based reconstruction and analysis of
biochemical reaction systems. Therefore **MASSpy**  can be used seperately from
or in conjuction with [**COBRApy**][11], thereby providing a wide range of modeling
workflows and techniques. Additional information about [**COBRApy**][11] can be found in its
[documentation](https://cobrapy.readthedocs.io/en/latest/index.html>).

## Installation

Check out the following information from the [Quick Start](https://masspy.readthedocs.io/en/latest/installation/quickstart.html) guide to get started using MASSpy!

### With Python 3.6+

The recommended method is to install **MASSpy** is to use ``pip`` to
install the software from the [Python Package Index](https://pypi.org/project/masspy/)
It is recommended to do this inside a [virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/)::

    pip install "masspy"

To install **MASSpy** with visualization dependencies:

    pip install "masspy[visualization]"

To install **MASSpy** with all optional dependencies:

    pip install "masspy[all]"

### With Docker
To quickly get started with the latest version of MASSpy using Docker, run the following commands in a shell:

    docker pull sbrg/masspy
    docker run --rm \
        --mount type=volume,src=licenses,dst=/home/masspy_user/opt/licenses \
        --mount type=volume,src=mass_project,dst=/home/masspy_user/mass_project \
        --publish 8888:8888 \
        -it sbrg/masspy

From within the container, either run ``python`` or ``jupyter notebook --ip=0.0.0.0 --port=8888`` depending on
the desired Python workspace. Don't forget to change the port number if it is already being used!

### Additional installation details

For additional details about how to set up an environment for MASSpy, including how to set up commercial optimizers and  optional dependencies, check out the detailed [Installation guide](https://masspy.readthedocs.io/en/latest/installation/quickstart.html) in the documentation!

## Documentation
The documentation for **MASSpy** is found at https://masspy.rtfd.io with installation instructions and several tutorials for getting started. All documentation is generated using [Sphinx](https://www.sphinx-doc.org/) and hosted by [ReadTheDocs](https://readthedocs.org/).

### Cite

To cite the **MASSpy** software publication:

> Haiman ZB, Zielinski DC, Koike Y, Yurkovich JT, Palsson BO (2021)
> MASSpy: Building, simulating, and visualizing dynamic biological models in Python using mass action kinetics.
> PLOS Computational Biology 17(1): e1008208. https://doi.org/10.1371/journal.pcbi.1008208

Additionally, please consider citing **COBRApy**, **libRoadRunner**, and other software dependencies of MASSpy! Citations and links to several dependencies as well as other useful literature sources are found in the [Works Cited](https://masspy.readthedocs.io/en/latest/references.html) and [Code Repositories](https://masspy.readthedocs.io/en/latest/additional/code_repositories.html) sections of the documentation.

## Support

Not sure how to [file an issue](.github/SUPPORT.md), want to [contribute](.github/CONTRIBUTING.md) to MASSpy, or just looking for some [general guidance](.github/FAQ.md)? Check out the [Support page](.github/SUPPORT.md)!

## License

The **MASSpy** source is released under the [MIT license](https://github.com/SBRG/MASSpy/blob/main/LICENSE). However, please note that binary packages (e.g., GLPK, CPLEX, etc.) and other dependencies (e.g. [openCOBRA packages](https://opencobra.github.io/), [libRoadRunner](http://libroadrunner.org/), etc.) will be bound by their respective license agreements as well.

[1]: https://img.shields.io/github/v/release/sbrg/masspy?label=MASSpy&sort=semver&style=plastic
[2]: https://img.shields.io/pypi/pyversions/masspy?logo=python&style=plastic
[3]: https://img.shields.io/github/workflow/status/sbrg/masspy/CI-CD?label=Tests&logo=GitHub%20Actions&style=plastic
[4]: https://img.shields.io/readthedocs/masspy?label=docs&logo=Read%20the%20Docs&style=plastic
[5]: https://img.shields.io/codecov/c/github/sbrg/masspy?logo=codecov&style=plastic
[6]: https://img.shields.io/github/license/sbrg/masspy?style=plastic
[7]: https://img.shields.io/docker/image-size/sbrg/masspy?label=Docker%20Img&logo=Docker&sort=semver&style=plastic
[8]: https://img.shields.io/badge/code%20style-black-000000.svg?style=plastic
[9]: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white&style=plastic
[10]: https://img.shields.io/badge/DOI-10.1371%2Fjournal.pcbi.1008208-blue?style=plastic
[11]: https://github.com/opencobra/cobrapy
