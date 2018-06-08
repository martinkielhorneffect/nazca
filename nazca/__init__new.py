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
#-----------------------------------------------------------------------
#
# Main initialization


"""Nazca module for photonic integrated circuit design."""

#__all__ = ['Cell', 'Pin', 'Polygon', 'Polyline', 'Annotation', 'load_gds',
#    'export_gds', 'export_plt', 'add_layer', 'get_layer', 'show_layers',
#    'show_xsections', 'layout_open', 'cell_open', 'cell_close', 'layout_close',
#    'gds_filter', 'replaceCells', 'print_structure', 'image',
#    'strt', 'bend', 'taper', 'ptaper', 'cp', 'add_xs', 'text', 'linelength',
#    'textheight', 'textfont', '__version__', 'modulepath']

import pandas as pd

from . import cfg
from . import simglobal
from . import slabmode
from .slabsolver import *
from .xsection import *
from .plotting import *
from .gds import *
from .lyp import *
from .util import *
from .convert_AWGtxt2nazca import *
from .gds_import import *
from .gds_base import *
from .colormap import *

from . import netlist
from .layout import *

from . import geometries as geom
from . import mask_layers
from . import mask_elements
from .gds_pcellreplace import *
from .text import *
from .textfont import textfont

from . import interconnects as ic
from . import cp

from .pdk_template import *
from ._version import __version__
import os
# Nazca module path
modulepath = os.path.dirname(os.path.abspath(__file__))
