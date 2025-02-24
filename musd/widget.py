"""Analyse the data to achieve the purpose of musd"""
from typing import Dict
import logging
import numpy as np
from numpy.typing import ArrayLike
from .type_checks import PathStr
from . import freegsnke_util
import pickle
import pdb

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QLabel, QSizePolicy, QSlider, QSpacerItem, \
     QVBoxLayout, QWidget, QLineEdit, QFrame, QToolBar, QDoubleSpinBox, QPushButton, QProgressBar, \
     QCheckBox, QComboBox, QTableWidget, QHeaderView, QTableWidgetItem, QFileDialog, QScrollArea
from PyQt5.QtGui import QColor
from PyQt5 import QtCore

from scipy.ndimage import zoom

import importlib

import pyqtgraph as pg
from scipy.ndimage import zoom

#import freegsfast
from freegsnke import equilibrium_update, GSstaticsolver
from freegsnke import jtor_update # for initialising the profile object
from copy import deepcopy

logger = logging.getLogger("musd")
from copy import deepcopy


from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QFrame, QLabel,
                             QLineEdit, QSlider, QPushButton, QSplitter, QScrollArea,
                             QSizePolicy)
from PyQt5.QtCore import Qt
import pyqtgraph as pg
import numpy as np
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QFrame, QLabel,
                             QLineEdit, QSlider, QPushButton, QSplitter, QScrollArea,
                             QSizePolicy, QSpacerItem)  # Added QSpacerItem
from PyQt5.QtCore import Qt
import pyqtgraph as pg
import numpy as np

class Widget(QWidget):
    def __init__(self, shot=None, time=None, parent=None, config=None, vc_dir=None):
        super(Widget, self).__init__(parent=parent)

        self.num_vcs = 4
        self.num_coils = 12
        self.equil_psin = [0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2]
        self.config = config
        self.vc_dir = vc_dir

        # Some variables that will be needed throughout
        self.symmetric_machine = False
        self.shot = shot
        self.shot_time = time
        self.vc_coeffs = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.coil_list = ['p1', 'p4', 'p5', 'px', 'd1', 'd2', 'd3', 'd5', 'd6', 'd7', 'dp', 'pc']
        self.coil_mins = [-48.0, -11.0, -11.0, -5.4, -8.0, -8.0, -6.5, -4.0, -3.3, -4.6, -6.4, 0]
        self.coil_maxs = [48.0, 0.0, 0.0, 5.4, 8.0, 8.0, 6.5, 4.0, 3.3, 4.6, 6.4, 0]
        self.coil_dIdt = np.array([1820.0,
                                   800.0, 900.0, 800.0,
                                   800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0,
                                   0]) / 364
        self.equil_valid = True

        # Sets the top level layout, will first organise things vertically
        self.verticalLayout = QVBoxLayout(self)
        # NEW: Set spacing and margins to 0
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)

        # Adding a bar across the top to define the shot number and shot time
        self.shot_toolbar = QFrame()
        # NEW: Set fixed height for the toolbar
        self.shot_toolbar.setFixedHeight(50)
        self.toolbar_layout = QHBoxLayout(self.shot_toolbar)
        # NEW: Set content margins for the toolbar
        self.toolbar_layout.setContentsMargins(5, 5, 5, 5)

        # Adds a box to select the shot number
        self.shot_label = QLabel('Shot: ')
        self.shot_edit = QLineEdit(str(self.shot))
        self.shot_edit.returnPressed.connect(self.load_equilibrium)

        # Adds a slider to select the time
        self.shot_time_label = QLabel('Time (s): ')
        self.shot_time_slider = QSlider(Qt.Horizontal)
        self.shot_time_edit = QLineEdit(str(self.shot_time))

        with open("/home/charlie/CLEANmusd/musd_config/efit_times.pickle", "rb") as pickle_file:
            self.efit_times = pickle.load(pickle_file)
        with open("/home/charlie/CLEANmusd/musd_config/efit_status.pickle", "rb") as pickle_file:
            self.efit_status = pickle.load(pickle_file).data
        self.shot_time_slider.setMinimum(0)
        self.shot_time_slider.setMaximum(len(self.efit_times) - 1)
        #indx = np.where(self.efit_times >= self.shot_time)[0][0]
        indx = 0
        self.shot_time_slider.setValue(indx)
        self.shot_time_slider.valueChanged.connect(self.update_shot_time)
        self.shot_time_slider.sliderReleased.connect(self.load_equilibrium)
        self.shot_time_slider.sliderReleased.connect(self.update_coil_currents_label)
        self.shot_time_edit.returnPressed.connect(self.update_shot_time_slider)
        self.shot_time_edit.returnPressed.connect(self.update_shot_time)
        self.shot_time_edit.returnPressed.connect(self.load_equilibrium)
        self.shot_time_edit.returnPressed.connect(self.update_coil_currents_label)

        # Adding a shot recalculate button
        self.calc_equilibrium = QPushButton('Recalculate')
        self.calc_equilibrium.clicked.connect(self.recalculate_equilibrium)

        # Setting the upper toolbar layout
        self.toolbar_layout.addWidget(self.shot_label)
        self.toolbar_layout.addWidget(self.shot_edit)
        self.toolbar_layout.addWidget(self.shot_time_label)
        self.toolbar_layout.addWidget(self.shot_time_slider)
        self.toolbar_layout.addWidget(self.shot_time_edit)
        self.toolbar_layout.addWidget(self.calc_equilibrium)

        # Adding the toolbar to the top of the widget
        self.verticalLayout.addWidget(self.shot_toolbar)

        # Create a splitter for the main content
        self.splitter = QSplitter(Qt.Horizontal)
        self.splitter.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        # Setting the plot window with a white background
        self.win = pg.GraphicsLayoutWidget(title='Equilibrium')
        self.win.setBackground('w')
        self.splitter.addWidget(self.win)

        # Create a scroll area for the pf_stack
        self.pfScrollArea = QScrollArea()
        self.pfScrollArea.setWidgetResizable(True)
        self.pfScrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.pfScrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.splitter.addWidget(self.pfScrollArea)

        # Set up the pf_stack
        self.pf_stack = QFrame()
        self.pf_stack_layout = QHBoxLayout(self.pf_stack)
        self.pfScrollArea.setWidget(self.pf_stack)

        # Create a widget for coil controls
        self.coilControlsWidget = QWidget()
        self.coilControlsLayout = QVBoxLayout(self.coilControlsWidget)

        # Setting the virtual circuit toolbar to load a virtual circuit
        self.vc_toolbar = QFrame()
        self.vc_toolbar_layout = QHBoxLayout(self.vc_toolbar)

        self.vc_file_label = QLabel('VC file:')
        self.vc_file_edit = QLineEdit('./')
        self.vc_file_edit.returnPressed.connect(self.set_vc_file)
        self.vc_load_button = QPushButton('Load')
        self.vc_load_button.clicked.connect(self.click_load_vc)

        self.vc_toolbar_layout.addWidget(self.vc_file_label)
        self.vc_toolbar_layout.addWidget(self.vc_file_edit)
        self.vc_toolbar_layout.addWidget(self.vc_load_button)

        # Generating the toolbar for the coil sliders
        self.I_toolbar = QFrame()
        self.I_toolbar_layout = QHBoxLayout(self.I_toolbar)

        self.I_lock_button = QPushButton('Lock')
        self.I_lock_button.clicked.connect(self.lock_coil_ratio)
        self.I_mark_button = QPushButton('Mark')
        self.I_mark_button.clicked.connect(self.mark_equil_state)
        self.I_reset_button = QPushButton('Reset')
        self.I_reset_button.clicked.connect(self.reset_equil)

        self.I_toolbar_layout.addWidget(self.I_lock_button)
        self.I_toolbar_layout.addWidget(self.I_mark_button)
        self.I_toolbar_layout.addWidget(self.I_reset_button)

        self.coilControlsLayout.addWidget(self.vc_toolbar)
        self.coilControlsLayout.addWidget(self.I_toolbar)

        self.pf_name_label = QLabel('Coil currents:')
        self.coilControlsLayout.addWidget(self.pf_name_label)

        # Setting the coil sliders
        self.coil_stack = QFrame()
        self.coil_stack_layout = QVBoxLayout(self.coil_stack)
        self.coil_slider_dict = {}
        for c, coil in enumerate(self.coil_list):
            self.coil_slider_dict[coil] = {}
            self.coil_slider_dict[coil]['label'] = QLabel(f'{coil}: 0kA')
            self.coil_slider_dict[coil]['slider'] = QSlider(Qt.Horizontal)
            self.coil_slider_dict[coil]['toolbar'] = QFrame()
            self.coil_slider_dict[coil]['layout'] = QHBoxLayout(self.coil_slider_dict[coil]['toolbar'])

            if coil == 'pc' or coil == 'p1':
                self.coil_slider_dict[coil]['slider'].setMinimum(0)
                self.coil_slider_dict[coil]['slider'].setMaximum(0)
                self.coil_slider_dict[coil]['slider'].setValue(0)
            else:
                self.coil_slider_dict[coil]['slider'].setMinimum(0)
                self.coil_slider_dict[coil]['slider'].setMaximum(10000)
                self.coil_slider_dict[coil]['slider'].setValue(5000)
            self.coil_slider_dict[coil]['slider'].sliderReleased.connect(self.update_coil_currents_edit)
            self.coil_slider_dict[coil]['slider'].sliderReleased.connect(self.recalculate_equilibrium)

            self.coil_slider_dict[coil]['edit'] = QLineEdit('0')
            self.coil_slider_dict[coil]['edit'].returnPressed.connect(self.update_coil_currents_slider)
            self.coil_slider_dict[coil]['edit'].returnPressed.connect(self.recalculate_equilibrium)

            self.coil_slider_dict[coil]['layout'].addWidget(self.coil_slider_dict[coil]['label'])
            self.coil_slider_dict[coil]['layout'].addWidget(self.coil_slider_dict[coil]['slider'])
            self.coil_slider_dict[coil]['layout'].addWidget(self.coil_slider_dict[coil]['edit'])

            self.coil_stack_layout.addWidget(self.coil_slider_dict[coil]['toolbar'])

        # Adding a slider for VC drive
        self.vc_drive_slider_dict = {}
        self.vc_drive_slider_dict['label'] = QLabel(f'VC drive:')
        self.vc_drive_slider_dict['slider'] = QSlider(Qt.Horizontal)
        self.vc_drive_slider_dict['toolbar'] = QFrame()
        self.vc_drive_slider_dict['layout'] = QHBoxLayout(self.vc_drive_slider_dict['toolbar'])

        self.vc_drive_slider_dict['slider'].setMinimum(0)
        self.vc_drive_slider_dict['slider'].setMaximum(10000)
        self.vc_drive_slider_dict['slider'].setValue(5000)
        self.vc_drive_slider_dict['slider'].sliderReleased.connect(self.drive_vc_currents_slider)
        self.vc_drive_slider_dict['slider'].sliderReleased.connect(self.recalculate_equilibrium)

        self.vc_drive_slider_dict['edit'] = QLineEdit('0')
        self.vc_drive_slider_dict['edit'].returnPressed.connect(self.drive_vc_currents_edit)
        self.vc_drive_slider_dict['edit'].returnPressed.connect(self.recalculate_equilibrium)

        self.vc_drive_slider_dict['layout'].addWidget(self.vc_drive_slider_dict['label'])
        self.vc_drive_slider_dict['layout'].addWidget(self.vc_drive_slider_dict['slider'])
        self.vc_drive_slider_dict['layout'].addWidget(self.vc_drive_slider_dict['edit'])

        # Add VC drive slider to the stack
        self.coil_stack_layout.addWidget(self.vc_drive_slider_dict['toolbar'])

        self.coilControlsLayout.addWidget(self.coil_stack)

        # Create a scroll area for the marks box
        self.marksScrollArea = QScrollArea()
        self.marksScrollArea.setWidgetResizable(True)
        self.marksScrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.marksScrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        # Adding in a side bar to show the mark points
        self.marks_box = QFrame()
        self.marks_box.setStyleSheet('background-color: white')
        self.marks_box_layout = QVBoxLayout(self.marks_box)
        self.marksScrollArea.setWidget(self.marks_box)

        # Add coil controls and marks box to the pf_stack layout
        self.pf_stack_layout.addWidget(self.coilControlsWidget, stretch=7)
        self.pf_stack_layout.addWidget(self.marksScrollArea, stretch=3)

        # Set size policies
        self.win.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.pfScrollArea.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.coilControlsWidget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.marksScrollArea.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)

        # Add the splitter to the layout
        self.verticalLayout.addWidget(self.splitter)
        # Set the stretch factor for the vertical layout
        self.verticalLayout.setStretch(0, 0)  # shot_toolbar
        self.verticalLayout.setStretch(1, 1)  # splitter

        # Set initial sizes for the splitter
        self.splitter.setSizes([self.width() // 2, self.width() // 2])

        self.equil_orig = None
        self.equil_mod = None
        self.equil_lim = None
        self.equil_lock = None
        self.mark_dict = {}
        self.mark_labels = [None for i in self.efit_times if True]  # Creating a none list to hold all the labels
        self.iso1 = None
        self.iso2 = None
        self.iso3 = None
        self.p1 = None

        if time is not None:
            self.update_shot_time()
            self.load_equilibrium()
            self.update_coil_currents_label()


    ################################
    from .freegsnke_util import (set_env,
                                recalculate_equilibrium,
                                load_equilibrium,
                                plot_equilibrium,
                                load_efit_times_and_status,
                                load_static_solver_inputs,
                                )

    from .gui_util import (update_coil_currents_label,
                           update_coil_currents_slider,
                           update_coil_currents_edit,
                           update_shot_time,
                           reset_equil,
                           click_load_vc,
                           showOpenFileDialog,
                           lock_coil_ratio,
                           drive_vc_currents_edit,
                           drive_vc_currents_slider,
                           set_vc_file,
                           mark_equil_state,
                           update_shot_time_slider,
                           )


    set_env(config=None, symmetric_machine=False)
    from freegsnke import build_machine
    ############################### Imports that I wanted to separate for neatness
