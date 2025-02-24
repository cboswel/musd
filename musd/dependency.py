"""Routines for retrieving the signal dependencies of the analysis code."""
from typing import List
import logging


logger = logging.getLogger("musd")


def check_base_signal(signal_code: str) -> bool:
    """Check the base signal is part of the declared dependencies

    Parameters
    ----------
    signal_code: str
        The 3-letter code/tag for the signal being retrieved.
    
    Returns
    -------
    bool:
        True if the signal code is in the existing dependencies, False if not.
    """
    # TODO check the signal code against the yaml file metadata. See
    # https://jwodder.github.io/kbits/posts/pypkg-data/#accessing-package-data-at-runtime
    return True


def get_dependencies() -> List[str]:
    """Generate the fully qualified signal dependencies needed for our analysis

    Returns
    -------
    List[str]
        The list of fully qualified signal dependencies

    Raises
    ------
    RuntimeError
        If the any of the base dependencies used are not in the declared
        dependencies of this package
    """
    # xbm = bolometry
    base_signal = "XBM"
    leaf_signals = ["PHASE", "AMP", "POWER"]
    if not check_base_signal(base_signal):
        message = f"Signal dependency {base_signal} not in declared signal dependencies."
        logger.error(message)
        raise RuntimeError(message)
    dependencies = [f"/{base_signal}/CORE/F05/{leaf}" for leaf in leaf_signals]
    return dependencies
