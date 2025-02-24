#!/bin/bash

if [ ! -d "venv" ]; then
  python -m venv venv
fi

source venv/bin/activate

pip install freegsnke/.
pip install .
pip install freegs4e
