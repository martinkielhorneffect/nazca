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
cellnames = dict()
cells = [] # list of all cells that have defined (cells not instances)
self = None # active cell reference
topcells = [] # Cells that have been collected to be parsed for GDS export

def active_cell():
    """Get the active cell."""
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
    'origin': ['user']
}
default_xs_width = 1.0
default_xs_radius = 20.0

#==============================================================================
# default mask layers
#==============================================================================
default_dump_layer = 1001
default_layer_AnnotationPin = default_dump_layer
default_layer_FiberIO = default_dump_layer



layerdict = dict()

xsmap = None #map scriptnames of xs to technology names of xs.

share = set() # set of cellname not to prefix.
#some cells may enter the design through different sub-gds files.
#There may be case where the can be considered to be the same cell.

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
    }


bb_pin_shape = 'arrow_full'
bb_pin_size = 1.0

bbox_pin_size = 0.75
bbox_stubs = True # add stubs to the cell bounding box

cellname_max_height = 50
cellname_scaling = 0.5


# =============================================================================
# pin shapes
# =============================================================================
#normalized to 1
pinshapes = {
    'arrow_full': [(0, 0), (-0.5, -0.5), (-0.5, -0.25), (-0.7, -0.25), (-0.7, 0.25), (-0.5, 0.25), (-0.5, 0.5)],
    'arrow_lefthalf': [(0, 0), (-0.5, 0.5), (-0.5, 0.25), (-0.7, 0.25), (-0.7, 0)],
    'arrow_righthalf': [(0, 0), (-0.5, -0.5), (-0.5, -0.25), (-0.7, -0.25), (-0.7, 0)],
    'circle': [(0.5, 0), (0.353, 0.353), (0, 0.5), (-0.353, 0.353), (-0.5, 0),
               (-0.353, -0.353), (0, -0.5), (0.353, -0.353)]
}


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


# =============================================================================
# Tracing elements
# =============================================================================
_trace_id = -1
_trace_ids = set()
traces = defaultdict(list)

def trace_start(trace_id=None):
    global _trace_id, _trace_ids
    _trace_id += 1
    if trace_id is None:
        trace_id = _trace_id
    if trace_id not in _trace_ids:
        _trace_ids.add(trace_id)
    else:
        print("Warning: starting already active trace_id = {} in trace_start.".\
           format(trace_id))
    if trace_id in traces.keys():
        print("Warning: reusing trace_id = {} in trace_start\n"\
              "  This will extend the trace, not start a new one!".format(trace_id))

def trace_stop(trace_id=None):
    global _trace_id, _trace_ids
    if trace_id is None:
         trace_id = _trace_id
    if trace_id in _trace_ids:
        _trace_ids.remove(trace_id)
    else:
        print("Warning: Trying to stop already stopped trace id={}.".\
            format(trace_id))

def trace_length(trace_id=None):
    if trace_id is None:
        trace_id = _trace_id
    length = 0
    for elm in traces[trace_id]:
        try:
            L = elm.cnode.cell.length_geo
            length += L
            #print(elm, L)
        except:
            pass
    return length


