"""Tests of the calculation module

NOTE: These are very basic tests that you should improve for your particular
case. And indeed as you add and modify the functions in the module, the tests
below will also need to be updated.
"""
from musd import calculation
import numpy as np
import numpy.testing as npt
import pytest
import json
import pathlib


@pytest.fixture
def xbm_example_data():
    """Some cached example data as would be returned by read_data.get_data()

    This data corresponds to the first 1000 elements of
    /XBM/CORE/F05/{AMP,PHASE,POWER} signals as returned by read_data.get_data()
    for pulse 47160/0.
    It can be regenerated as follows in Python:

    .. code:: python

       from musd.read_data import get_data
       import json
       signals = [ "/XBM/CORE/F05/PHASE", "/XBM/CORE/F05/POWER", "/XBM/CORE/F05/AMP"]
       core_f05 = get_data(signals, 47160)
       core_f05 = {key: value[:1000].tolist() for key, value in core_f05.items()}
       with open("core_f05.json", "w") as fp:
           json.dump(core_f05, fp)
    
    """
    data_file = pathlib.Path(__file__).parent / "core_f05.json"
    with open(data_file, "r") as fp:
        core_f05 = json.load(fp)
    return {key: np.array(value) for key, value in core_f05.items()}



def test_do_science(xbm_example_data):
    """Test that the do_science function returns the correct calculated values

    NOTE The check below if very simplistic: an absolute normalised quantity
    should have a maximum of 1. You will need to select properties to match your
    use case and test against those. Alternative, you could go more of a
    regression test route and directly save to disk the values that you would
    expect from the calculation and compare against those. Generally, a mix of
    both types of test is best.
    """
    results = calculation.do_science(dependency_data=xbm_example_data, config_dir=None, calib_dir=None)
    results_max = [value.max() for key, value in results.items() if key not in ("time", "signals")]
    npt.assert_allclose(results_max, 1.0)
