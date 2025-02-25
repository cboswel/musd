"""Analyse the data to achieve the purpose of musd"""
from typing import Dict
import logging
import numpy as np
from numpy.typing import ArrayLike
from .type_checks import PathStr
import os
import pickle
import pdb

logger = logging.getLogger("musd")

#import freegsfast
from freegsnke import equilibrium_update, GSstaticsolver
from freegsnke import jtor_update # for initialising the profile object
from copy import deepcopy
import pyqtgraph as pg
from scipy.ndimage import zoom

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QLabel, QSizePolicy, QSlider, QSpacerItem, \
     QVBoxLayout, QWidget, QLineEdit, QFrame, QToolBar, QDoubleSpinBox, QPushButton, QProgressBar, \
     QCheckBox, QComboBox, QTableWidget, QHeaderView, QTableWidgetItem, QFileDialog
from PyQt5.QtGui import QColor
from PyQt5 import QtCore

def set_env(config=None, symmetric_machine=False):
    # freegsnke builds the tokamak using freegs and environment variables
    # that point to pickle files. You can't import the build_machine
    # function until you've set those environment variables

    if config is None:
        config_path = os.path.abspath(__file__)
        config_str = config_path.split('/')

        config = "/"
        for s, string in enumerate(config_str[0:-2]):
            config += string
            config += "/"

        config = config + "musd_config/"

    config = "./musd_config/"
    # Now that the files are made we can set the environment variables
    # Personally I'm not sure I like this way of setting up the machine, but
    # I don't really have the time to faff with it right now
    pckl_files = [os.path.normpath(f"{config}/MAST-U_active_coils.pickle"),
                  os.path.normpath(f"{config}/MAST-U_active_coils_nonsym.pickle"),
                  os.path.normpath(f"{config}/MAST-U_passive_coils.pickle"),
                  os.path.normpath(f"{config}/MAST-U_wall.pickle"),
                  os.path.normpath(f"{config}/MAST-U_limiter.pickle"),
                  os.path.normpath(f"{config}/MAST-U_magnetic_probes.pickle")]
    print(pckl_files)
    pckl_exist = True
    for i in pckl_files:
        if pckl_exist:
            pckl_exist = os.path.isfile(i)
        else:
            print("Missing some pickle files, generating before setting env")
            break

    if pckl_exist is False:
        # First have a function to write the pickle files, can check this by selecting
        # plotting true
        from .create_MASTU_pickles import build_machine_pickles
        build_machine_pickles(path=config, split_passives=True)

    if symmetric_machine:
        os.environ["ACTIVE_COILS_PATH"] = f"{config}/MAST-U_active_coils.pickle"
    else:
        os.environ["ACTIVE_COILS_PATH"] = f"{config}/MAST-U_active_coils_nonsym.pickle"

    os.environ["PASSIVE_COILS_PATH"] = f"{config}/MAST-U_passive_coils.pickle"
    os.environ["WALL_PATH"] = f"{config}/MAST-U_wall.pickle"
    os.environ["LIMITER_PATH"] = f"{config}/MAST-U_limiter.pickle"
    os.environ["PROBE_PATH"] = f"{config}/MAST-U_magnetic_probes.pickle"

# Function to caluculate a modified equilibrium using additional currents from the sliders
def recalculate_equilibrium(self):

    coil_names = ['Solenoid', ['p4_upper', 'p4_lower'], ['p5_upper', 'p5_lower'],
                  ['px_lower', 'px_upper'], ['d1_lower', 'd1_upper'], ['d2_lower', 'd2_upper'],
                  ['d3_lower', 'd3_upper'], ['d5_lower', 'd5_upper'], ['d6_lower', 'd6_upper'],
                  ['d7_lower', 'd7_upper'], ['dp_lower', 'dp_upper'], None]

    if self.equil_orig is not None:
        if self.equil_lock is not None:
            self.equil_mod = deepcopy(self.equil_lock)
            tmp_equil = self.equil_lock
            print("Modifying new equil")
        else:
            self.equil_mod = deepcopy(self.equil_orig)
            tmp_equil = self.equil_orig
            print("Modifying original equil")


        for i in np.arange(len(coil_names)):

            if type(coil_names[i]) is str:
                while i < len(self.coil_list):
                    if self.coil_list[i] == 'p1':
                        c = i
                        break
                    i += 1

                coil_drive = self.vc_coeffs[c]
                self.equil_mod.tokamak[coil_names[i]].current = tmp_equil.tokamak[coil_names[i]].current + coil_drive*1000

            if type(coil_names[i]) is list:
                for j in np.arange(len(coil_names[i])):
                    coil_label = coil_names[i][j].split('_')

                    while i < len(self.coil_list):
                        if self.coil_list[i] == coil_label[0]:
                            c = i
                            break
                        i += 1

                    coil_drive = self.vc_coeffs[c]
                    self.equil_mod.tokamak[coil_names[i][j]].current = tmp_equil.tokamak[coil_names[i][j]].current + coil_drive*1000

    # carry out the forward solve
    self.solver.solve(eq=self.equil_mod,
                      profiles=self.profiles,
                      constrain=None,
                      target_relative_tolerance=1e-6)

    self.plot_equilibrium()
    self.update_coil_currents_label()

# Loading the efit equilibrium data and running the freegs dynamic solve
def load_equilibrium(self):

    # Read the shot number and time
    self.shot = self.shot_edit.text()

    if (self.shot == '0') or (self.shot_time == '0'):
        print('No shot or time entered!')
    else:

        print('-- starting freeGSNKE reconstruction --')
        print('Shot:'+self.shot)
        print('Time:'+str(self.shot_time))

        # Build the tokamak object for freegs
        self.tokamak = self.build_machine.tokamak()

        # load the parameters
        Ip, fvac, alpha, beta, alpha_logic, beta_logic, currents_sym, currents_nonsym, _ = self.load_static_solver_inputs()


        # load the efit reconstruction timestamps and convergence status
        # Note: efit_times that have efit_status=-1 did not converge so we will exclude those time slices
        efit_times = self.efit_times
        efit_status = self.efit_status

        i = np.argmin(np.abs(float(self.shot_time) - efit_times))

        # If there are mark points then want to make sure that they are carried forward, we're assuming
        # that we're maintaining the other plasma parameters because the shape changes are small
        coil_names = ['Solenoid', ['p4_upper', 'p4_lower'], ['p5_upper', 'p5_lower'],
                      ['px_lower', 'px_upper'], ['d1_lower', 'd1_upper'], ['d2_lower', 'd2_upper'],
                      ['d3_lower', 'd3_upper'], ['d5_lower', 'd5_upper'], ['d6_lower', 'd6_upper'],
                      ['d7_lower', 'd7_upper'], ['dp_lower', 'dp_upper'], None]

        if len(self.mark_dict.keys()) > 0:
            # We have some time keys so let's update the currents
            prev_mark = None
            for k, key in enumerate(self.mark_dict.keys()):
                if float(key) > self.shot_time:
                    break
                else:
                    prev_mark = key

            if prev_mark is not None:
                for c in np.arange(len(coil_names)):
                    if type(coil_names[c]) is str:
                        if self.symmetric_machine:
                            currents_sym[coil_names[c]][i::] += self.mark_dict[prev_mark]['coil_diffs'][coil_names[c]]
                        else:
                            currents_nonsym[coil_names[c]][i::] += self.mark_dict[prev_mark]['coil_diffs'][coil_names[c]]
                    if type(coil_names[c]) is list:
                        for j in np.arange(len(coil_names[c])):
                            if self.symmetric_machine:
                                currents_sym[coil_names[c][j]][i::] +=self.mark_dict[prev_mark]['coil_diffs'][coil_names[c][j]]
                            else:
                                currents_nonsym[coil_names[c][j]][i::] +=self.mark_dict[prev_mark]['coil_diffs'][coil_names[c][j]]

        # figure out which time slices did not converge
        time_indices = np.where(efit_status == 1)[0]
        time_slices_excluded = np.where(efit_status == -1)[0]
        print(f"{len(time_slices_excluded)} time slices (out of total {len(efit_times)}) excluded from simulation.")
        print(f"The excluded time slices are: {efit_times[time_slices_excluded]} seconds.")

        if self.shot_time in efit_times[time_slices_excluded]:
            print(f"Selected time is excluded, not bothering to simulate")
            self.equil_valid = False
            self.plot_equilibrium()
            return None
        else:
            self.equil_valid = True

        times = efit_times[time_indices]  # these are the slices we simulate
        print(f"Total time slices to be simulated: {len(times)}.")

        #time_indices = time_indices[0:-1]
        time_indices = [np.argmin(np.abs(float(self.shot_time) - efit_times))]
        times = efit_times[time_indices]
        n = len(times)

        # equilibrium object (note that both nx and ny have to be of the form 2**n + 1 with n being an integer)
        eq = equilibrium_update.Equilibrium(
        tokamak=self.tokamak,        # sets up the object with the MAST-U tokamak
        Rmin=0.06, Rmax=2.0,         # computational grid radial limits (same as EFIT++)
        Zmin=-2.2, Zmax=2.2,         # computational grid vertical limits (same as EFIT)
        nx=65,                       # number of grid points in the radial direction
        ny=65,                       # number of grid points in the vertical direction
        psi=None                     # initial guess for the plasma flux (can provide one if available)
        )

        # static solver object (used for solving later on)
        self.solver = GSstaticsolver.NKGSsolver(eq)

        # strings for dictionary keys
        time_str = str(efit_times[i]).replace(".", "")
        name = f"shot_{self.shot}_time_{time_str}"

        print(f"-----Solving equilibrium for shot {self.shot} at time {efit_times[i]}-----")

        # initialise profile object
        self.profiles = jtor_update.Lao85(
            eq=eq,                                          # equilibrium object
            limiter=self.tokamak.limiter,                   # plasma limiter
            Ip=Ip[i],                                       # total plasma current
            fvac=fvac[i],                                   # f vacuum parameter (R*Bt)
            alpha=alpha[i,:],                               # p' profile coefficients
            beta=beta[i,:],                                 # ff' profile coefficients
            alpha_logic=bool(alpha_logic[i]),               # logic parameters from above
            beta_logic=bool(beta_logic[i]),
        )

        # set coil currents in eq object
        # checks if machine is up/down symmetric or not
        if self.symmetric_machine:
            for key in currents_sym.keys():
                eq.tokamak[key].current = currents_sym[key][i]
        else:
            for key in currents_nonsym.keys():
                eq.tokamak[key].current = currents_nonsym[key][i]

        # carry out the forward solve
        self.solver.solve(eq=eq,
                          profiles=self.profiles,
                          constrain=None,
                          target_relative_tolerance=1e-6)

        if len(self.mark_dict.keys()) > 0:
            self.equil_lock = deepcopy(eq)
        else:
            self.equil_orig = deepcopy(eq)

        print("-----SIMULATIONS COMPLETE!-----")

        self.plot_equilibrium()

def plot_equilibrium(self):

    if self.p1 is None:
    # Set up the plot
        self.p1 = self.win.addPlot(title="<b>Equilibrium Reconstruction<\b>", col=0, row=0, rowspan=6)
        vb = self.p1.getViewBox()
        vb.setAspectLocked(lock=True)
        self.p1.clear()
        self.efit_error_label = pg.TextItem("", anchor=(-1.0,-1.0))
        self.efit_error_label.setParentItem(self.p1.graphicsItem())

    if self.equil_lim is None:
        with open(os.path.normpath("./musd_config/equil_r_lim.pickle"), "rb") as pickle_file:
            self.equil_r_lim = list(pickle.load(pickle_file))
        with open(os.path.normpath("./musd_config/equil_z_lim.pickle"), "rb") as pickle_file:
            self.equil_z_lim = list(pickle.load(pickle_file))
        self.equil_lim = self.p1.plot(self.equil_r_lim, self.equil_z_lim, pen='grey')

    if self.equil_valid is False:
        self.efit_error_label.setText("EFIT time invalid", color='red')
        return None
    else:
        self.efit_error_label.setText("EFIT time valid", color='blue')

    if self.equil_orig is not None:

        self.p1.setXRange(np.min(self.equil_orig.R), np.max(self.equil_orig.R))
        self.p1.setYRange(np.min(self.equil_orig.Z), np.max(self.equil_orig.Z))

        #print('Equil_r')
        #print(self.equil_orig.R)
        #print('Equil_z')
        #print(self.equil_orig.Z)

        equil_aspect_ratio = (np.max(self.equil_orig.Z) - np.min(self.equil_orig.Z)) / (np.max(self.equil_orig.R) - np.min(self.equil_orig.R))
        dz = np.abs(self.equil_orig.Z[0,1]-self.equil_orig.Z[0,0])

        #print('Aspect ratio '+str(equil_aspect_ratio))
        #print('dz '+str(dz))

        psi = self.equil_orig.psi()
        psi_n = (psi-self.equil_orig.psi_axis) / (self.equil_orig.psi_bndry - self.equil_orig.psi_axis)

        if self.iso1 is None:
            self.iso1 = []
            for i in np.arange(len(self.equil_psin)):

                self.iso1 = self.iso1 + [pg.IsocurveItem(data=zoom(psi_n,[1,1.0*equil_aspect_ratio]), level=self.equil_psin[i])]

                self.iso1[i].resetTransform()

                self.iso1[i].setScale((np.max(self.equil_orig.R) - np.min(self.equil_orig.R)) / len(self.equil_orig.R))
                self.iso1[i].setTransformOriginPoint(np.min(self.equil_orig.R), np.min(self.equil_orig.Z)-dz)
                self.iso1[i].setParentItem(self.p1)

                self.iso1[i].pen.setColor(pg.mkColor('k'))
                if self.equil_psin[i] > 1.0:
                    self.iso1[i].pen.setStyle(Qt.DashLine)
                elif self.equil_psin[i] < 1.0:
                    self.iso1[i].pen.setStyle(Qt.DashDotLine)

                self.p1.addItem(self.iso1[i])

        else:
            for i in np.arange(len(self.equil_psin)):
                self.iso1[i].setData(zoom(psi_n,[1,1.0*equil_aspect_ratio]))

    if self.equil_mod is not None:

        psi = self.equil_mod.psi()
        psi_n = (psi-self.equil_mod.psi_axis) / (self.equil_mod.psi_bndry - self.equil_mod.psi_axis)

        if self.iso2 is None:
            self.iso2 = []
            for i in np.arange(len(self.equil_psin)):
                self.iso2 = self.iso2 + [pg.IsocurveItem(data=zoom(psi_n,[1,1.0*equil_aspect_ratio]), level=self.equil_psin[i])]

                self.iso2[i].resetTransform()

                self.iso2[i].setScale((np.max(self.equil_mod.R) - np.min(self.equil_mod.R)) / len(self.equil_mod.R))
                self.iso2[i].setTransformOriginPoint(np.min(self.equil_mod.R), np.min(self.equil_mod.Z)-dz)
                self.iso2[i].setParentItem(self.p1)

                self.iso2[i].pen.setColor(pg.mkColor('r'))
                if self.equil_psin[i] > 1.0:
                    self.iso2[i].pen.setStyle(Qt.DashLine)
                elif self.equil_psin[i] < 1.0:
                    self.iso2[i].pen.setStyle(Qt.DashDotLine)

                self.p1.addItem(self.iso2[i])

        else:
            for i in np.arange(len(self.equil_psin)):
                self.iso2[i].setData(zoom(psi_n,[1,1.0*equil_aspect_ratio]))

    if self.equil_lock is not None:

        psi = self.equil_lock.psi()
        psi_n = (psi-self.equil_lock.psi_axis) / (self.equil_lock.psi_bndry - self.equil_lock.psi_axis)

        if self.iso3 is None:
            self.iso3 = []
            for i in np.arange(len(self.equil_psin)):
                self.iso3 = self.iso3 + [pg.IsocurveItem(data=zoom(psi_n,[1,1.0*equil_aspect_ratio]), level=self.equil_psin[i])]

                self.iso3[i].resetTransform()

                self.iso3[i].setScale((np.max(self.equil_lock.R) - np.min(self.equil_lock.R)) / len(self.equil_lock.R))
                self.iso3[i].setTransformOriginPoint(np.min(self.equil_lock.R), np.min(self.equil_lock.Z)-dz)
                self.iso3[i].setParentItem(self.p1)

                self.iso3[i].pen.setColor(pg.mkColor('b'))
                if self.equil_psin[i] > 1.0:
                    self.iso3[i].pen.setStyle(Qt.DashLine)
                elif self.equil_psin[i] < 1.0:
                    self.iso3[i].pen.setStyle(Qt.DashDotLine)

                self.p1.addItem(self.iso3[i])

        else:
            for i in np.arange(len(self.equil_psin)):
                self.iso3[i].setData(zoom(psi_n,[1,1.0*equil_aspect_ratio]))



def load_efit_times_and_status(self):
    """
    Extract the shot status (converged or not) and times from the EFIT data.

    """

    # load data
    #status = client.get('/epm/equilibriumstatusinteger', self.shot)
    status = 1

    return status.time.data, status.data

def load_static_solver_inputs(self, zero_passives=False):
    """
    Extract the EFIT++ data required for the static solve.

    """

    # load data
    #Ip = client.get('/epm/input/constraints/plasmacurrent/computed', self.shot).data                           # plasma current
    #fvac = client.get('/epm/input/bvacradiusproduct', self.shot).data                                          # fvac
    #alpha = client.get('/epm/output/numericaldetails/degreesoffreedom/pprimecoeffs', self.shot).data           # pprime coefficients
    #beta = client.get('/epm/output/numericaldetails/degreesoffreedom/ffprimecoeffs', self.shot).data           # ffprime coefficients
    #alpha_logic = client.get('/epm/input/numericalcontrols/pp/edge', self.shot).data                           # pprime logical
    #beta_logic = client.get('/epm/input/numericalcontrols/ffp/edge', self.shot).data                           # ffprime logical

    # active coil\passive structure currents need to be done carefully
    #current_labels = client.get('/epm/input/constraints/pfcircuits/shortname', self.shot).data          # active/passive coil current names
    #currents_values = client.get('/epm/input/constraints/pfcircuits/computed', self.shot).data          # active/passive coil current values

    with open(os.path.normpath("./musd_config/Ip.pickle"), "rb") as pickle_file:
        Ip = np.vstack(list(pickle.load(pickle_file).values()))
    with open(os.path.normpath("./musd_config/fvac.pickle"), "rb") as pickle_file:
        fvac = np.vstack(list(pickle.load(pickle_file).values()))
    with open(os.path.normpath("./musd_config/alpha.pickle"), "rb") as pickle_file:
        alpha = np.vstack(list(pickle.load(pickle_file).values()))
    with open(os.path.normpath("./musd_config/beta.pickle"), "rb") as pickle_file:
        beta = np.vstack(list(pickle.load(pickle_file).values()))
    with open(os.path.normpath("./musd_config/alpha_logic.pickle"), "rb") as pickle_file:
        alpha_logic = list(pickle.load(pickle_file).values())
    with open(os.path.normpath("./musd_config/beta_logic.pickle"), "rb") as pickle_file:
        beta_logic = list(pickle.load(pickle_file).values())
    with open(os.path.normpath("./musd_config/current_lables.pickle"), "rb") as pickle_file:
        current_labels = list(pickle.load(pickle_file).values())
    with open(os.path.normpath("./musd_config/currents_values.pickle"), "rb") as pickle_file:
        currents_values = np.vstack(list(pickle.load(pickle_file).values()))


    # Active coils
    currents = {}
    currents_nonsym = {}
    currents_discrepancy = {}

    with open(os.environ["ACTIVE_COILS_PATH"], 'rb') as file:
        active_coils = pickle.load(file)

    efit_names = current_labels[0:24]   # active coil names in efit

    # loop through active coil names in freegsnke to set them with UDA data
    for active_coil_name in active_coils.keys():

        # special cases for specific coil names
        if active_coil_name == "Solenoid":
            # find indices for Solenoid coils in efit_names
            indices = [i for i, efit_name in enumerate(efit_names) if "p1" in efit_name]
            currents["Solenoid"] = currents_values[:, indices[0]] if indices else None
            currents_nonsym["Solenoid"] = currents_values[:, indices[0]] if indices else None
            currents_discrepancy["Solenoid"] = 0.0
        else:
            # find indices for other coils in efit_names
            indices = [i for i, efit_name in enumerate(efit_names) if active_coil_name in efit_name]
            if indices:
                polarity = np.sign(currents_values[:, indices[0]])
                average_current = np.sum(np.abs(currents_values[:, indices]), axis=1) / len(indices)
                currents[active_coil_name] = polarity*average_current
                currents_discrepancy[active_coil_name] = polarity*np.abs(np.abs(currents_values[:, indices[0]]) - average_current)
                for j in indices:
                    currents_nonsym[efit_names[j]] = currents_values[:, j]

            else:
                currents[active_coil_name] = None

    # passive structures
    with open(os.environ["PASSIVE_COILS_PATH"], 'rb') as file:
        passive_coils = pickle.load(file)

    # Passive structures
    for i in range(0, len(passive_coils)):
        coil = passive_coils[i]

        if "efitGroup" in coil:
            group_name = coil['efitGroup']
            ind = current_labels.index(group_name)

            if zero_passives:
                currents[coil["name"]] = 0.0
                currents_nonsym[coil["name"]] = 0.0
            else:
                currents[coil["name"]] = currents_values[:,ind]*coil["current_multiplier"]
                currents_nonsym[coil["name"]] = currents_values[:,ind]*coil["current_multiplier"]
            # currents[coil["name"]] = efit_currents['currents_input'][t,ind]
        else:
            group_name = coil['element']
            ind = current_labels.index(group_name)
            if zero_passives:
                currents[coil["name"]] = 0.0
                currents_nonsym[coil["name"]] = 0.0
            else:
                currents[coil["name"]] = currents_values[:,ind]*coil["current_multiplier"]
                currents_nonsym[coil["name"]] = currents_values[:,ind]*coil["current_multiplier"]
            # currents[coil["name"]] = efit_currents['currents_input'][t,ind]


    return Ip, fvac, alpha, beta, alpha_logic, beta_logic, currents, currents_nonsym, currents_discrepancy

