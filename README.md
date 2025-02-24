# MUSD

Python package for shot design on MAST-U

## Purpose

Tool has been created to help design coil current paths for scenario development on MAST-U. This has been adapted from a code made by James Harrison.

## Usage

Tool opens a gui with sliders to adjust coil currents. Development is on going to be able to lock ratios with a virtual circuit and then drive that instead.

To run this code you will need access to the freegsnke repo. For the moment on friea this can be found at ~cvincent/FreeGS/freegsnke. You will also need to run this module using python 3.9

To install first create a python virtual environment:
`python3 -m venv venv`

Then activate the virtual environment:
`source venv/bin/activate`

First install freegsnke, this can likely be done by going to ~cvincent/FreeGS/freegsnke and running
`pip install ".[freegsfast]"`

Upgrade your pip and install the package:
```
python3 -m pip install --upgrade pip
python3 -m pip install -e .

It has been pointed out to me that install freegsnke from my directory doesn't work. Alternatively you can source the virtual environment from my area and that should hopefully work.
`source ~cvincent/MUSD/musd/venv/bin/activate`

## Metadata

Important information for running this code on the scheduler is found in
`scheduler_dependencies.yml`. Python dependencies when using this code as a library are
specified in `pyproject.toml`. Pinned dependencies are specified in
`requirements.txt` for reproducible execution.
