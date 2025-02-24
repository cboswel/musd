"""
Test that the analysis code is able to produce netCDF outputs for the test shots
provided in the scheduler information file
"""
from musd import run_analysis
import yaml
import pathlib
import subprocess
import pytest
import time
import warnings


# Ignore code coverage due to integration tests (this would artificially inflate
# the coverage of the individual units of code)
pytestmark = pytest.mark.no_cover


def get_test_shots():
    metadata_file = pathlib.Path(__file__).parents[2] / "scheduler_dependencies.yml"
    with open(metadata_file, "r") as fp:
        metadata = yaml.safe_load(fp)
    custom_error = RuntimeError(
        "You must provide shot numbers with which to test your analysis code."
    )
    try:
        test_shots = metadata["test_shots"]
    except KeyError as exception:
        raise custom_error from exception
    else:
        if test_shots is None or test_shots is []:
            raise custom_error
    return test_shots

# Loop through the test shots, generate netCDF files for each, and test them
@pytest.mark.parametrize("test_shot", get_test_shots())
def test_run_analysis(test_shot):
    output_dir = pathlib.Path("./analysis_output")
    # start timer
    tic = time.perf_counter()
    #######################################################################
    #######################################################################
    # NOTE code authors will need to adjust this function call based on any
    # modifications they have made
    #######################################################################
    #######################################################################
    run_analysis(
        shotno=test_shot,
        passno=0,
        config="../musd_config",
        calib="../musd_calib",
        output=output_dir,
        verbose=True,
    )
    # stop timer
    toc = time.perf_counter()
    # issue any warning about elapsed time
    time_elapsed = toc - tic
    if time_elapsed > 60.0:
        raise RuntimeError(f"Analysis script has taken longer than 1 minute for shot {test_shot}")
    elif time_elapsed > 15.0:
        warnings.warn(f"Analysis script has taken longer than 15 seconds for shot {test_shot}")
    # check the netCDF file produced by analysis
    subprocess.check_call(["test", "-e", output_dir / f"musd{test_shot:06d}.nc"])
