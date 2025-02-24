"""Tests of the read_data module

NOTE: These are very basic tests that you should improve for your particular
case. And indeed as you add and modify the functions in the module, the tests
below will also need to be updated.
"""
from musd import read_data
import pyuda
import numpy as np


def test_get_data_single_signal():
    """Test that the `get_data` function returns the object we expect for a
    single signal
    """
    signal = "/XBM/CORE/F05/PHASE"
    phase = read_data.get_data(signal, 47160)
    assert isinstance(phase, dict)
    assert isinstance(phase[signal], np.ndarray)
    assert isinstance(phase["time"], np.ndarray)
    assert len(phase["time"]) == 50000


def test_get_data_list_of_signals():
    """Test that the `get_data` routine works with a list of signals
    """
    signals = [
        "/XBM/CORE/F05/PHASE",
        "/XBM/CORE/F05/POWER",
        "/XBM/CORE/F05/AMP",
    ]
    core_f05 = read_data.get_data(signals, 47160)
    assert isinstance(core_f05, dict)
    assert np.all(np.array([isinstance(x, np.ndarray) for x in core_f05.values()]))
