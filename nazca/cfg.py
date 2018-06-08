
#-----------------------------------------------------------------------
# This file is part of Nazca.
#
# Nazca is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# Nazca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Nazca.  If not, see <http://www.gnu.org/licenses/>.
#
# 2017 (c)  Ronald Broeke
#-----------------------------------------------------------------------

from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd


#==============================================================================
# plotting related
#==============================================================================
def formatplot():
     #update matplotlib plotting styles
     lines = {'linewidth':3}
     plt.rc('lines', **lines)
     font = {'size':14}
     plt.rc('font', **font)


matplotlib_font = {
    #'family': 'normal',
    'style' : 'normal',
    'weight': 'normal', #'light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'
    'size'  : 12}


#==============================================================================
# figsize for mode field plots
#==============================================================================
modeplotsize = (8,5)


#==============================================================================
# DRC related variables
#==============================================================================
drc_xs = 500
drc_ring_xs = (12, 1)


#==============================================================================
# Cell related variables
#==============================================================================
defaultcellname = 'nazca'
cells = [] # open cell stack. 'active cell' is cells[-1]
self = None # active cell reference
cellnames = dict() # lsit of all cells {cellname: Cell}
basenames = dict() # {basenames" function_id}
topcells = [] # Cells that have been collected to be parsed for GDS export
validate_basename = False # perform validate basenames on substrings if True.
solve_direct = True

def active_cell():
    """Get the active cell.

    The 'active cell' is the last opened Cell object that has not been closed yet.

    Returns:
        Cell: active cell
    """
    if cells:
        return cells[-1]
    else:
        return None

cp = None

xsall = dict() #foundry layer table
xs_layers = dict() #foundry layer table

default_xs_name = 'nazca'
default_xs = {
    'name': ['nazca'],
    'layer': [1111],
    'datatype': [0],
    'accuracy' : [0.005],
    'growx': [0],
    'growy': [0],
    'origin': ['nazca']
}
default_xs_width = 1.0
default_xs_radius = 20.0

default_xserror_name = 'error'
default_xserror = {
    'name': ['error'],
    'layer': [1111],
    'datatype': [1],
    'accuracy' : [0.005],
    'growx': [0],
    'growy': [0],
    'origin': ['nazca']
}
default_xserror_width = 1.0
default_xserror_radius = 20.0


#==============================================================================
# default mask layers
#==============================================================================
default_dump_layer = (1111, 0)
default_error_layer = (1111, 1)
default_layer_AnnotationPin = default_dump_layer
default_layer_FiberIO = default_dump_layer


layerdict = dict()

xsmap = None #map scriptnames of xs to technology names of xs.

share = set() # set of cellname not to prefix.
#some cells may enter the design through different sub-gds files.
#There may be cases where the can be considered to be the same cell.

stubmap = {}

default_df_xs = pd.DataFrame(default_xs)
mask_layers = {default_xs_name: default_df_xs}
XSdict = dict()
#Need to add
#XSdict[default_xs_name] = default_xs_name

# =============================================================================
# pin layers
# =============================================================================
default_layers = {
    'bb_pin':            (1001, 0),
    'bb_pin_text':       (1002, 0),
    'bb_parameter_text': (1003, 0),

    'package_pin':       (1011, 0),
    'package_pin_text':  (1012, 0),

    'bbox':              (1021, 0),
    'bbox_name':         (1022, 0),
    'bbox_pin':          (1023, 0),
    'bbox_pin_text':     (1024, 0),

    'pin_text':          (1004, 0)
    }

cellname_max_height = 50
cellname_scaling = 0.5

stub_default_annotation_layer = 'bb_pin_text'
bbox_stubs = True # add stubs to the cell bounding box
store_pins = False # store pins in cell annotation to read gds as BB with pin reconstruction.

# =============================================================================
# pin shapes
# =============================================================================
# shapes are normalized to 1
pinshapes = {
    'arrow_full': [(0, 0), (-0.5, -0.5), (-0.5, -0.25), (-0.7, -0.25),
        (-0.7, 0.25), (-0.5, 0.25), (-0.5, 0.5)],
    'arrow_lefthalf': [(0, 0), (-0.5, 0.5), (-0.5, 0.25), (-0.7, 0.25), (-0.7, 0)],
    'arrow_righthalf': [(0, 0), (-0.5, -0.5), (-0.5, -0.25), (-0.7, -0.25), (-0.7, 0)],
    'circle': [(0.5, 0), (0.353, 0.353), (0, 0.5), (-0.353, 0.353), (-0.5, 0),
        (-0.353, -0.353), (0, -0.5), (0.353, -0.353)]
}


def reset_pin_settings():
    global pin_settings
    pin_settings = {
        'bb_pin_shape': 'arrow_full',
        'bb_pin_size': 1.0,
        'bb_pin_layer':  default_layers['bb_pin'],
        'bb_pin_annotation_layer': default_layers['bb_pin_text'],
        'bb_pin_annotation_move': (-0.3, 0),
        'bbox_pin_size': 1.0,
        'bbox_pin_layer': default_layers['bbox_pin'],
        'stub_length': 2.0
        }
reset_pin_settings()
overrule_pdk_pinstyle = None #set before loading a foundry to update pin representation
documentation_pin_layer = (60000, 0)

# =============================================================================
# pin adb stub properties
# =============================================================================
pin_instantiate = False
stub_instantiate = False


#==============================================================================
# Font default
#==============================================================================
default_font = 'nazca'

#==============================================================================
# Matplotlib layout
#==============================================================================
plt_cmap_name = 'Set2' #colormap used in mask layout in Matplotlib
plt_cmap = [] # list of the colors in the cmap
plt_alpha = 0.3 # tranparency
plt_figsize = 8
plt_fontsize = 8
plt_background_inside = '#FFFFFF' #background color inside the axes
plt_background_outside = '#EEEEEE' #background color outside the axes


#==============================================================================
# spt output
#==============================================================================
spt = False

def gds_cellname_cleanup(name):
    return name


