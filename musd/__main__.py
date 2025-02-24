"""This is the main module which runs the musd analysis code."""
import fire
import logging
import pathlib
from musd import type_checks
from musd import calculation
from musd import dependency
from musd import widget
import sys
import pickle
import pyqtgraph as pg

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QLabel, QSizePolicy, QSlider, QSpacerItem, \
     QVBoxLayout, QWidget, QLineEdit, QFrame, QToolBar, QDoubleSpinBox, QPushButton, QProgressBar, \
     QCheckBox, QComboBox, QTableWidget, QHeaderView, QTableWidgetItem, QFileDialog
from PyQt5.QtGui import QColor
from PyQt5 import QtCore

logger = logging.getLogger("musd")
PathStr = type_checks.PathStr

@type_checks.enforce_types
def run_gui(
        shot: int = 51118,
        time: float = 0.3,
        config: PathStr = None,
        vc_dir: PathStr = None,
        debug: bool = False,
        verbose: bool = False,
):
    """Run the MUSD inter-shot analysis

    Parameters
    ----------
    shot : int
        Shot/pulse number
    time : float
        Time of equilibrium
    config : PathStr, optional
        Directory containing configuration files, by default "{module_path}/musd_config"
    vc_dir : PathStr, optional
        Directory containing the virtual circuits you want to test, by default your current directory
    debug : bool, optional
        Show debugging messages (--debug) or not (--nodebug, default)
    verbose : bool, optional
        Show all messages (--verbose) or not (--noverbose, default)
    """
    # setup
    if debug:
        logger.setLevel("DEBUG")
    elif verbose:
        logger.setLevel("INFO")
    else:
        logger.setLevel("WARNING")
    # output_dir = pathlib.Path(output)
    # if not output_dir.is_dir():
    #     logger.info(f"Creating output directory {output_dir}")
    #     output_dir.mkdir(parents=True)
    # main components of analysis

    app = QApplication(sys.argv)
    w = widget.Widget(shot=shot, time=time, config=config)
    w.show()
    sys.exit(app.exec_())


def main():
    """Interface to run the analysis from the command line"""
    fire.Fire(run_gui)


if __name__ == "__main__":
    main()
