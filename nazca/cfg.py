
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
# Matplotlib layout settings
#==============================================================================
plt_cmap_name = 'Set2' #colormap used in mask layout in Matplotlib
plt_alpha = 0.3 # tranparency
plt_figsize = 8 # size of layout
plt_fontsize = 10 #fontsize of annotations in the layout
plt_background_inside = '#FFFFFF' #background color inside the axes
plt_background_outside = '#EEEEEE' #background color outside the axes
plt_cmap = [] # list of the colors in the cmap


# font setting for matplotlib generated mask output
matplotlib_font = {
    #'family': 'normal',
    'style' : 'normal',
    'weight': 'normal', #'light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'
    'size'  : plt_fontsize}


def formatplot():
     """update matplotlib plotting style for graphs."""
     lines = {'linewidth':3}
     plt.rc('lines', **lines)
     font = {'size': plt_fontsize}
     plt.rc('font', **font)


# figsize for mode field plots
modeplotsize = (8, 5)


#==============================================================================
# Layout Settings
#==============================================================================
# The name of the cell that is generated automatically as the first topcell.
defaultcellname = 'nazca'


# add stubs to the cell bounding box
bbox_stubs = True


# store pins in cell annotation to read gds as BB with pin reconstruction.
store_pins = False


# store the cell name as cell annotation
store_bbname = True


# Redirect an unknown layer to the 'dump' layer.
# Set True for PDKs where only known layers should be allowed.
# Default = False works best for tutorials and examples
# and create layers by using them.
redirect_unknown_layers = False


# Dictionary for default layers for special purpose objects:
# They are set in a dict so they can be adjusted during runtime
# and are easily accessible by 'key'.
default_layers = {
    'bb_pin':            (1001, 0),
    'bb_pin_text':       (1002, 0),
    'bb_parameter_text': (1003, 0),
    'bb_name':           (1004, 0),

    'package_pin':       (1011, 0),
    'package_pin_text':  (1012, 0),

    'bbox':              (1021, 0),
    'bbox_name':         (1022, 0),
    'bbox_pin':          (1023, 0),
    'bbox_pin_text':     (1024, 0),

    'pin_text':          (1005, 0),
    'dump':              (1111, 0), # redirect content for unknown-layers here
    'error':             (1111, 1), # error layer
    'docu_pin':          (6000, 0), # Layer to store pin shapes for documentation purposes
    'hull':              (1111, 2) # layer to (optionally) show the convex hull in.
    }


# default xs to create when no xsection is set:
default_xs_name = 'nazca'
default_xs = {
    'name': [default_xs_name],
    'layer': [default_layers['dump'][0]],
    'datatype': [default_layers['dump'][1]],
    'accuracy' : [0.005],
    'growx': [0],
    'growy': [0],
    'origin': ['nazca']
}
default_xs_width = 1.0
default_xs_radius = 20.0


# default error layer to place errornous structures, e.g. impossible interconnects:
default_xserror_name = 'error'
default_xserror = {
    'name': [default_xserror_name],
    'layer': [default_layers['error'][0]],
    'datatype': [default_layers['error'][1]],
    'accuracy' : [0.005],
    'growx': [0],
    'growy': [0],
    'origin': ['nazca']
}
default_xserror_width = 1.0
default_xserror_radius = 20.0


# instantiate pins-shapes as separate cells. default=False
pin_instantiate = False


# instantiate stub cell as separate cells. default=False
stub_instantiate = False


# Font default
default_font = 'nazca'


# define pin shape dict. Shapes are normalized to 1
pinshapes = {
    'arrow_full': [(0, 0), (-0.5, -0.5), (-0.5, -0.25), (-0.7, -0.25),
        (-0.7, 0.25), (-0.5, 0.25), (-0.5, 0.5)],
    'arrow_lefthalf': [(0, 0), (-0.5, 0.5), (-0.5, 0.25), (-0.7, 0.25), (-0.7, 0)],
    'arrow_righthalf': [(0, 0), (-0.5, -0.5), (-0.5, -0.25), (-0.7, -0.25), (-0.7, 0)],
    'circle': [(0.5, 0), (0.353, 0.353), (0, 0.5), (-0.353, 0.353), (-0.5, 0),
        (-0.353, -0.353), (0, -0.5), (0.353, -0.353)]
}


# maximum height of the cell name
cellname_max_height = 50


# size of the cell_name in the cell w.r.t. the size of the cell
cellname_scaling = 0.5


# perform validate basenames on substrings if True.
validate_basename = False


# Solve the pin xya as soon as it is added to the layout.
solve_direct = True


# Clear the mask after an export or not, default = True
# If you want to export multiple independent layouts from a cell, set to True.
# In a Jupyter notebook a step by step is many times more intuitive, in that
# case set to False
export_clear = True


# add convex hull to mask export.
export_hull = False


# obsolete. generate spt output. default=False
spt = False


#==============================================================================
# DRC visualisation related variables
#==============================================================================
drc_xs = 500
drc_ring_xs = (12, 1)


# =============================================================================
# initialize global variables into existence
# =============================================================================
cp = None
self = None # active cell reference
cells = [] # store open cells
cellnames = dict() # lsit of all cells {cellname: Cell}
basenames = dict() # {basenames" function_id}
topcells = [] # Cells that have been collected to be parsed for GDS export
xsall = dict() #foundry layer table
xs_layers = dict() #foundry layer table
layerdict = dict()
xsmap = None #map scriptnames of xs to technology names of xs.
share = set() # set of cellname not to prefix.
stubmap = {}
default_df_xs = pd.DataFrame(default_xs)
mask_layers = {default_xs_name: default_df_xs}
XSdict = dict()


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


#TODO: clean up with try except
def gds_cellname_cleanup(name):
    """Placeholder function name that can be overruled in a foundry"""
    return name


