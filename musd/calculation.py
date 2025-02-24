"""Analyse the data to achieve the purpose of musd"""
from typing import Dict
import logging
import numpy as np
from numpy.typing import ArrayLike
from .type_checks import PathStr


logger = logging.getLogger("musd")


def do_science(
    dependency_data: Dict[str, ArrayLike], config_dir: PathStr, calib_dir: PathStr
) -> Dict[str, ArrayLike]:
    """Manipulate the input data to extract some desired quantity or derived
    data
    
    NOTE Nothing is done in the simple calculation below with regards to
    configuration and calibration (i.e. the config_dir and calib_dir arguments).
    You will need to implement this for your particular use case. You can either
    put some additional functions within this module, or create your own
    calibration and configuration modules.

    Parameters
    ----------
    dependency_data : Dict[str, ArrayLike]
        The input data from the dependent signals
    config_dir : PathStr
        Path to the configuration directory
    calib_dir : PathStr
        Path to the calibration directory

    Returns
    -------
    Dict[str, ArrayLike]
        The output data for our calculation. The keys are the signal name, and
        the values are the data as Numpy arrays.
    """
    logger.info("Starting calculation...")
    logger.debug(
        f"Calculation arguments: \n\tdependency_data = {dependency_data}\n\t"
        f"config_dir = {config_dir}\n\t"
        f"calib_dir = {calib_dir}"
    )
    output_data = {}
    output_data["signals"] = []
    output_data["time"] = dependency_data["time"]
    data_keys = [x for x in dependency_data.keys() if x != "time"]
    for key in data_keys:
        output_key = key.replace("XBM", "MUSD")
        output_data["signals"].append(output_key)
        # You might want to take absolute normalisation your data. This is
        # nonsense in this example, but you get the point.
        data = np.abs(dependency_data[key])
        output_data[output_key] = data / data.max()
    logger.info("Completed calculation...")
    logger.debug(f"Calculation output: {output_data}")

    return output_data
