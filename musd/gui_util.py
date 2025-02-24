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

def update_coil_currents_label(self):
   if self.equil_mod is None:
       current = self.equil_orig.tokamak.getCurrents()
   else:
       current = self.equil_mod.tokamak.getCurrents()
   coil_current = 0
   for c, coil in enumerate(self.coil_list):
       if coil == 'p1':
           coil_current = current['Solenoid']/1000
           # Set the coil label to the correct value
           self.coil_slider_dict[coil]['label'].setText(f'{coil}: {coil_current:4.2f}kA')
           if (coil_current > self.coil_maxs[c]) or (coil_current < self.coil_mins[c]):
              self.coil_slider_dict[coil]['label'].setStyleSheet("background-color: red")
           else:
              self.coil_slider_dict[coil]['label'].setStyleSheet("background-color: none")

       elif coil == 'pc':
           pass
       else:
           components = 0
           for I, Icoil in enumerate(current.keys()):
               if coil.lower() in Icoil.lower():
                   if 'case' in Icoil.lower():
                       pass
                   else:
                       coil_current += current[Icoil]
                       components += 1
           coil_current = (coil_current / components)/1000
           self.coil_slider_dict[coil]['label'].setText(f'{coil}: {coil_current:4.2f}kA')
           if (coil_current > self.coil_maxs[c]) or (coil_current < self.coil_mins[c]):
               self.coil_slider_dict[coil]['label'].setStyleSheet("background-color: red")
           else:
               self.coil_slider_dict[coil]['label'].setStyleSheet("background-color: none")

       coil_current = 0

def update_coil_currents_slider(self):
    for c, coil in enumerate(self.coil_list):
        if coil == 'pc' or coil == 'p1':
           continue
        current_val = float(self.coil_slider_dict[coil]['edit'].text()) # In kA
        if abs(current_val) > self.coil_dIdt[c]:
            if current_val == abs(current_val):
                current_val = self.coil_dIdt[c]
            else:
                current_val = -1*self.coil_dIdt[c]
        current_val = current_val/self.coil_dIdt[c]
        slider_val = (current_val*5000) + 5000
        self.coil_slider_dict[coil]['slider'].setValue(int(slider_val))

def update_coil_currents_edit(self):
    for c, coil in enumerate(self.coil_list):
        if coil == 'pc' or coil == 'p1':
            continue
        current_val = float(self.coil_slider_dict[coil]['slider'].value()) # In kA
        text_val = ((current_val/5000)-1)*self.coil_dIdt[c]
        self.coil_slider_dict[coil]['edit'].setText(f'{text_val:4.2f}')
        self.vc_coeffs[c] = text_val

def update_shot_time(self):
    if (self.efit_times is not None) and (self.efit_status) is not None:
        self.shot_time = self.efit_times[self.shot_time_slider.value()]
        self.shot_time_label.setText("Time: {:5.3f} s".format(self.shot_time))
    else:  # No equilibrium available
        self.shot_time_label.setText("Time:")

def update_shot_time_slider(self):
   shot_time_text = self.shot_time_edit.text()
   indx = np.where(self.efit_times >= float(shot_time_text))[0][0]
   self.shot_time_slider.setValue(indx)

def reset_equil(self):
    print( '----- Resetting modified equilibrium -----')
    # Clear the modified equilibriums
    self.equil_mod = None
    self.equil_lock = None

    # Remove the modified equilibrium from the plot
    if self.iso2 is not None:
        for i in np.arange(len(self.iso2)):
            self.p1.removeItem(self.iso2[i])
    self.iso2 = None

    # Reset the coil current labels
    self.update_coil_currents_label()

    # Delete the mark points for this time
    if f'{self.shot_time:4.3f}' in self.mark_dict.keys():
           del self.mark_dict[f'{self.shot_time:4.3f}']

    if len(self.mark_dict.keys()) == 0:
        self.equil_lock = None
        self.load_equilibrium()


def click_load_vc(self):

    # Invoke the file open dialog box to get the file name
    try:
       vc_filename = self.vc_file_edit.text()
       vc_coeffs = np.loadtxt(vc_filename)
    except:
       vc_filename = self.showOpenFileDialog()
       if vc_filename is None:
          return None

       # Open the VC file
       vc_coeffs = np.loadtxt(vc_filename)

    # Scale VC's to kAs
    vc_coeffs = vc_coeffs/1000

    for c, coil in enumerate(self.coil_list):
        self.vc_coeffs[c] = vc_coeffs[c]
        self.coil_slider_dict[coil]['edit'].setText(f'{vc_coeffs[c]:4.3f}')

    self.update_coil_currents_slider()

def set_vc_file(self):
    # Set the virtual circuit file name from the vc edit box
    vc_file_path = self.vc_file_edit.text()

    try:
        f = open(vc_file_path)
    except FileNotFoundError:
        # Need to add a proper exception here
        print('Please enter a proper file name')
        return None

    vc_file_str = vc_file_path.split('/')
    self.vc_dir = "/"
    for s, string in enumerate(vc_file_str[0:-1]):
        self.vc_dir += string
        self.vc_dir += "/"

def showOpenFileDialog(self):
    if self.vc_dir is None:
        self.vc_dir = os.path.abspath('./')

    fileName, filter = QFileDialog.getOpenFileName(self, 'Open file',
        self.vc_dir, 'Text files (*.txt)')
    if fileName:
        self.vc_file_edit.setText(f'{fileName:s}')
        self.vc_dir = "/"
        for s, string in enumerate(fileName[0:-1]):
           self.vc_dir += string
           self.vc_dir += "/"
        return fileName

def lock_coil_ratio(self):
    # Function to lock the current modification to make
    # changes on top before marking the equil state
    if self.equil_mod is not None:

        # Get the VC coefficients from the file to calculate the drive value
        # Invoke the file open dialog box to get the file name
        try:
           vc_filename = self.vc_file_edit.text()
           vc_coeffs = np.loadtxt(vc_filename)
        except:
           print('Saving manual coil ratios as virtual circuit')
           vc_filename = 'Manual'
           vc_coeffs = []
           for c, coil in enumerate(self.coil_list):
              vc_coeffs.append(float(self.coil_slider_dict[coil]['edit'].text()))

        # Next need to make a dictionary for this time if there isn't one already
        if f'{self.shot_time:4.3f}' not in self.mark_dict.keys():
           self.mark_dict[f'{self.shot_time:4.3f}'] = {'vc_conf': []}
           vc_configs = []
        else:
           # Checking if there are already conifgurations mark at this time
           vc_configs = self.mark_dict[f'{self.shot_time:4.3f}']['vc_conf']

        # Need to get the drive so get that from the VC_drive slider
        vc_drive = self.vc_drive_slider_dict['edit'].text()

        mark_dict = {'vc_file':vc_filename,
                     'vc_coeffs':vc_coeffs,
                     'vc_drive': vc_drive,
                     }
        vc_configs.append(mark_dict)
        self.mark_dict[f'{self.shot_time:4.3f}']['vc_conf'] = vc_configs

        # The current state of VC drive should now be saved in the mark_dict
        # Now need to save the current mod point so that it can be iterated on

        self.equil_lock = deepcopy(self.equil_mod)
        print('------Modification locked, continue editing-------')

    else:
        print('No point locking because equilibrium hasn\'t changed')

def drive_vc_currents_edit(self):
    # Function to use the drive edit button to drive the locked coil ratios
    driver_val = float(self.vc_drive_slider_dict['edit'].text()) # Unitless

    if abs(driver_val) > 2:
       if abs(driver_val) == driver_val:
          slider_val = 10000
       else:
          slider_val = 0
    else:
       slider_val = ((driver_val/2) * 5000) + 5000
    
    # update the vc_driver slider
    self.vc_drive_slider_dict['slider'].setValue(slider_val)

    for c, coil in enumerate(self.coil_list):
        current_val = float(self.coil_slider_dict[coil]['edit'].text())
        Icoil = current_val * driver_val
        self.vc_coeffs[c] = Icoil


def drive_vc_currents_slider(self):
    # Function to use the drive slider to drive the locked coil ratios
    slider_val = float(self.vc_drive_slider_dict['slider'].value())

    driver_val = (slider_val/5000) - 1
    print(slider_val, driver_val)

    self.vc_drive_slider_dict['edit'].setText(f'{driver_val:4.2f}')
    for c, coil in enumerate(self.coil_list):
        current_val = float(self.coil_slider_dict[coil]['edit'].text())
        Icoil = current_val * driver_val
        self.vc_coeffs[c] = Icoil


def mark_equil_state(self):
    # Function to mark the total modified state at this time so it can be
    # held for modification of the shot later in time
    # Specifically at this point we need to get the absolute current difference
    # between the marked state and the original state of the coil currents

    if self.equil_mod is None:
        print("No modifed equilibrium to mark, doing nothing")
        return None
    else:
        print(f'Marking modified equilibrium at time {self.shot_time:4.3f}')

    self.mark_dict[f'{self.shot_time:4.3f}']['coil_diffs'] = {}

    # In the mark dictionary setting the absolute changes of every coil
    coil_names = ['Solenoid', ['p4_upper', 'p4_lower'], ['p5_upper', 'p5_lower'],
                  ['px_lower', 'px_upper'], ['d1_lower', 'd1_upper'], ['d2_lower', 'd2_upper'],
                  ['d3_lower', 'd3_upper'], ['d5_lower', 'd5_upper'], ['d6_lower', 'd6_upper'],
                  ['d7_lower', 'd7_upper'], ['dp_lower', 'dp_upper'], None]
    for i in np.arange(len(coil_names)):
        if type(coil_names[i]) is str:
            old_current = self.equil_orig.tokamak[coil_names[i]].current
            new_current = self.equil_mod.tokamak[coil_names[i]].current
            diff_current = new_current - old_current
            self.mark_dict[f'{self.shot_time:4.3f}']['coil_diffs'][coil_names[i]] = diff_current

        if type(coil_names[i]) is list:
            for j in np.arange(len(coil_names[i])):
                old_current = self.equil_orig.tokamak[coil_names[i][j]].current
                new_current = self.equil_mod.tokamak[coil_names[i][j]].current
                diff_current = new_current - old_current
                self.mark_dict[f'{self.shot_time:4.3f}']['coil_diffs'][coil_names[i][j]] = diff_current

    # Now that I've marked my currents I can set the label in the marks box
    label_text = f'{self.shot_time:4.3f}:'
    for c in range(len(self.mark_dict[f'{self.shot_time:4.3f}']['vc_conf'])):
        conf = self.mark_dict[f'{self.shot_time:4.3f}']['vc_conf'][c]
        file_name = conf['vc_file']
        drive = float(conf['vc_drive'])
        label_text += f'\nVC file: {file_name:s}'
        label_text += f'\n    Drive: {drive:4.2f}'

    if self.mark_labels[self.shot_time_slider.value()] is None:
       self.mark_labels[self.shot_time_slider.value()] = QLabel(label_text)
       self.marks_box_layout.addWidget(self.mark_labels[self.shot_time_slider.value()])
    else:
       self.mark_labels[self.shot_time_slider.value()].setText(label_text)
