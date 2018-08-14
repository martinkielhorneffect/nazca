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

# You should have received a copy of the GNU Affero General Public License
# along with Nazca.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------
#==============================================================================
# (c) 2016-2017 Ronald Broeke, Katarzyna Lawniczuk
#==============================================================================

"""Module with a set Buiding Block templates for faciliating PDK creation."""

import os
from collections import OrderedDict
import inspect
from functools import wraps
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
from pprint import pprint
from PIL import Image

import hashlib
import nazca as nd
from nazca import cfg
#from nazca.util import parameters_to_string
import nazca.geometries as geom


#put here to avoid loading utils
def parameters_to_string(param):
    """Create a string from a parameter dictionary.

    param (dict): (parameter_name, value)

    Format:

    "Parameters:
    <parameter> = <value>
    <parameter> = <value>
    ..."

    Returns:
        str: parameters as a string
    """
    plist = ['Parameters:']
    for key, value in param.items():
        plist.append("{} = {}".format(key, value))
    return '\n'.join(plist)


#==============================================================================
# hash
#==============================================================================
_hash_id = None
_hash_params = {}
_hash_name = 'NONAME'
cfg._hashme = []


def hashme(*params):
    """Parametrize the decorator to send base name and parameter names.

    Returns:
        decorator
    """

    hash_class_warning =\
"""Warning: @hashme decorator on Class method:

By default Nazca will guarantee a unique cell name each time a cell is created.
It does so by suffixing cell names with an ordinal counter.
This may result in multiple copies of the same cell, though,
where only one would do.
Decorating the cell generating function with @hashme avoids cell copies of
identical cells by returning an existing cell if the function has been
called with the same parameters earlier.

The @hashme 'state checker' can only be used safely if the state of the
function it decorates is *not* dependent on a state *external* to that
function, i.e. if the variations of the cell it generates *only* depend
on the function parameters. In contrast, Class methods typically depend
on attributes stored at class level, hence, by their very nature Class methods
do not fit well with the @hashme concept to guarantee unique cells for
unique names.

Example of a cell function where @hashme should NOT be used, i.e. in a class method:

    class myElement():
        # stuff
        def make_cell(a, b, c):
            # stuff
            return cell_object

    elm = myElement()

Now compare two ways place a cell multiple times.

* method 1:
    elm.make_cell(a, b, c).put()
    elm.make_cell(a, b, c).put()

a) creates two copies of the cell and suffixes the name with an ordinal counter.
b) if @hashme is used it returns an existing cell if the call profile
   of make_cell has been used earlier.

* method 2:
A secure way in all cases to reuse cells is to call a cell generation function
only once and assign its returns value -the Cell object- to a variable.
Use this variable to put instances of the cells:

    E = elm.cell(a, b, c)
    E.put()
    E.put()
"""

    #print('params:' , params)
    Npar = len(params)
    if Npar == 0:
        print('warning: no name given to @hashme')
        name = 'NONAME'
    else:
        name = params[0]
        params = params[1:]

    def decorator(cellfunc):
        """Decorator for grabbing and hashing a function's full parameter list.

        While executing the wrapper a global hashid and parameter list are set.
        Before leaving the wrapper the hashid is reset to None.

        Returns:
            function: wrapper function
        """
        @wraps(cellfunc)
        def wrapper(*args, **kwargs):
            global _hash_id
            global _hash_params
            global _hash_name
            nonlocal name

            hash_length = 4
            name_long = name
            getargs = inspect.getargspec(cellfunc)
            funcargs = getargs.args
            funcdefaults = getargs.defaults
            hashstr = ''
            _hash_params = OrderedDict()

            if funcargs:
                # also make it work for class methods by removing 'self':
                skip = 0
                if funcargs[0] == 'self':
                    skip = 1
                    funcargs = funcargs[1:]
                    print(hash_class_warning)

                for i, (a, v) in enumerate(zip(funcargs, funcdefaults)):
                    try: # set to kwargs
                        _hash_params[a] = kwargs[a]
                    except:
                        try: # set to args
                            if a != 'self':
                                _hash_params[a] = args[i+skip]
                        except: # set to default
                            _hash_params[a] = v

                for p in sorted(_hash_params.keys()):
                    hashstr += '{}_{}'.format(p, _hash_params[p])
                for p in params:
                    if p in _hash_params.keys():
                        if isinstance(_hash_params[p], float):
                            #avoid formats like 0.0000000001 in names.
                            name_long += "_{:.5f}".format(_hash_params[p]).rstrip('0').rstrip('.')
                        else:
                            name_long += "_{}".format(_hash_params[p])

                _hash_id = hashlib.md5(hashstr.encode()).hexdigest()[:hash_length]
            else:
                _hash_id = ''

            if _hash_id is not '':
                _hash_name = "{}_${}".format(name_long, _hash_id)
            else:
                _hash_name = name_long

            # add cfg to avoid explicit import of this module for hashme attributes
            cfg.hash_id = _hash_id
            cfg.hash_params = _hash_params
            basename = cfg.gds_cellname_cleanup(name)
            cfg.hash_basename = basename                               # base
            cfg.hash_paramsname = cfg.gds_cellname_cleanup(name_long)  # base + params
            cfg.hash_name = cfg.gds_cellname_cleanup(_hash_name)       # base + params + hash
            cfg.hash_func_id = id(cellfunc)

            # check if basename comes from a unique function:
            try:
                funcid = cfg.basenames[basename]
                if funcid != cfg.hash_func_id:
                    print("ERROR: Reusing a basename across functions: '{}'".\
                        format(basename))
                    #raise
            except:
                if cfg.validate_basename:
                    nd.validate_basename(basename)
                # add basename:
                cfg.basenames[basename] = cfg.hash_func_id

            cell_exists = _hash_name in cfg.cellnames.keys()

            if cell_exists:
                return cfg.cellnames[_hash_name]
            else:
                cell_new = cellfunc(*args, **kwargs)

            #reset:
            _hash_id = None
            _hash_params = {}
            _hash_name = ''
            cfg.hash_id  = ''
            cfg.hash_params = ''
            cfg.hash_basename = ''
            cfg.hash_paramsname = ''
            cfg.hash_name = ''
            cfg.hash_func_id = None
            return cell_new
        return wrapper
    return decorator


def get_Cell(name):
    """Get Cell object by providing cell name.

    Args:
        name (str): cell name

    Returns:
        Cell: cell having cellname <name>
    """

    if name in cfg.cellnames.keys():
        return cfg.cellnames[name]
    else:
        raise Exception("Error: Requested cell with name '{}' not existing."\
            "Availabe are: {}.".format(name, sorted(cfg.cellnames.keys())))


def rangecheck(allowed_values):
    """Check if parameter values are in range.

    A dictionary with the allowed values for parameters is sent to rangecheck.
    If any parameter is out of range, a ValueError will be raised and relevant
    error info is provided to solve the issue.

    Args:
        allowed_values (dict): {'<var_name>': (low, <var_name>, high)}

    Raises:
        ValueError if out of range.

    Returns:
        None

    Example:

        Add a check in function `func` on 0 <= a <= 10 and -5 <= b <= 5.
        The call `func(a=100, b=100)` will raise a Value error::

            def func(a, b):
                nd.rangecheck({'a': (0, a, 10), 'b': (-5, b, 5)})

            func(a=100, b=100)

            # output:
            # ValueError: a=100 too large. Allowed values: 0<=a<=10
            # b=100 too large. Allowed values: -5<=b<=5
    """
    err = ''
    for var, (low, value, high) in allowed_values.items():
        if low is None:
            low = float('-inf')
        if high is None:
            high = float('inf')
        if low <= value <= high:
            continue
        elif value < low:
            err += "{0}={3} too small. Allowed values: {1}<={0}<={2}\n".\
                 format(var, low, high, value)
        elif value > high:
            err += "{0}={3} too large. Allowed values: {1}<={0}<={2}\n".\
                 format(var, low, high, value)
    if err != '':
        raise ValueError(err)
    return None

#==============================================================================
#
#==============================================================================
stubs = dict() # dict of all stubs. {stub_name: stub_cell_obj}
stubmap = dict() # dict of xsection_name to its stub's xsection_name. {xs_name: stub_xs_name}

class Functional_group():
    """Class to group building blocks syntactically together.

    Example:
        Create bb_splitters group with various splitter components
        such that the mmi can be found under bb_splitters.<tab>::

            bb_splitters = Functional_group()
            bb_splitters.mmi1x2 = mmi1x2
            bb_splitters.mmi2x2 = mmi2x2
            bb_splitters.mmi3x3 = mmi3x3
    """

    def __init__(self, name=''):
        """Create a Functional group object."""
        self.name = name
BB_Group = Functional_group


def parameters(parameters=None):
    """Put a parameter list as annotation in a building block.

    Args:
        parameters (dict): {parameter_name: value}

    Returns:
        None
    """
    cell = cfg.cells[-1]
    if isinstance(parameters, dict):
        cell.parameters = parameters
        anno = nd.Annotation(text=parameters_to_string(parameters),
            layer='bb_parameter_text')
    elif cell.hashme:
        anno = nd.Annotation(text=parameters_to_string(cell.parameters),
            layer='bb_parameter_text')
    else:
        anno = nd.Annotation(text='', layer='bb_parameter_text')
    return anno


def cellname(cellname=None, length=0, width=None, align='lc'):
    """Create the cellname as a text cell.

    Args:
        cellname (str): name of the cell
        length (float): length available for the BB name in um
        width (float):
        align (str): text alignment (see nazca.text, default = 'lc')

    Returns:
        Cell: text with cellname
    """
    cell = cfg.cells[-1]
    if cellname is None:
       cellname = cell.cell_paramsname

    if width is not None:
        maxheight = min(width*0.8, cfg.cellname_max_height)
    else:
        maxheight = cfg.cellname_max_height

    length = length * cfg.cellname_scaling
    texth = min(maxheight, abs(length) / nd.linelength(cellname, 1))
    txt = nd.text(cellname, texth, layer='bbox_name', align=align)
    return txt


cfg.__mapping_errors = []
def addBBmap(name, params=(), trans=(0, 0, 0)):
    """Apply building block mapping for Phoenix software.

    Returns:
        None

    Exceptions:
        ValueError
    """
    if not cfg.spt:
        return None

    if not isinstance(params, tuple):
        params = (params,)

    #get the right row for <name>:
    mapping = cfg.bbmap.loc[cfg.bbmap['Nazca'] == name]

    #TODO: check also for series with len>1, i.e. multiple matches
    try:
        od_name = mapping['OD_BBname'].iloc[0]
        od_port = mapping['OD_port'].iloc[0]
        odx, ody, oda = mapping[['odx', 'ody', 'oda']].iloc[0]
        x, y, a = mapping[['x', 'y', 'a']].iloc[0]
    except:
        od_name = 'NO_NAME:nazca_name:{})'.format(name)
        x, y, a = 0, 0, 0
        odx, ody, oda = 0, 0, 0
        od_port = 'NO_PORT'
        cfg.__mapping_errors.append("No valid mapping found for '{}'.".format(name))

    #store mapping-info in active cell
    cfg.cells[-1].BBmap = (od_name, od_port, (odx, ody, oda), (x, y, a), trans, params)


#==============================================================================
# Stubs
#==============================================================================
#@hashme('arrow', 'layer', 'pinshape', 'pinshape')
# Do NOT @hashme this function,
# because parameter values are finalized inside the function.
def make_pincell(layer=None, shape=None, size=None):
    """Create a cell to indicate a pin position.

    The cell contains a shape, e.g. an arrow, to point out a location in the layout.
    Available pin shapes are in a dictionary in cfg.pinshapes: {name: polygon}.
    The predefined shapes have been normalized to unit size.

    Args:
        layer (float): layer number to place the pin symbol/shape
        shape (str): name (dict key) of the pin symbol/shape.
        size (float): scaling factor of the pin symbol/shape

    Returns:
        Cell: cell with pin symbol
    """
    pinshape =  cfg.pin_settings['bb_pin_shape']
    if layer is None:
        layer = cfg.pin_settings['bb_pin_layer']
    layer = nd.get_layer(layer)

    if shape is None:
        shape = pinshape

    if size is None:
        size = cfg.pin_settings['bb_pin_size']

    if shape not in cfg.pinshapes.keys():
        print("Warning: arrowtype '{}' not recognized.".format(shape))
        print("Available options are: {}.".format(cfg.pinshapes.keys()))
        print("Fall back type '{}' will be used.".format(pinshape))
        shape = pinshape

    name = "{}_{}_{}_{}".format('arrow', layer, shape, size)
    name = cfg.gds_cellname_cleanup(name)
    if name in cfg.cellnames.keys():
        return cfg.cellnames[name]
    with nd.Cell(name, instantiate=cfg.pin_instantiate, store_pins=False) as arrow:
        outline = [(x*size, y*size) for x, y in cfg.pinshapes[shape]]
        nd.Polygon(layer=layer, points=outline).put(0)
    return arrow


def stubname(xs, width, thick, stubshape=None, pinshape=None, pinsize=None, pinlayer=None):
    """Construct a stub name.

    Args:
        xs (str); xsection name
        width (float): stub width
        thick (float): thickness of stub into cell (length)

    Returns:
        str: stub name
    """
    if width is None:
        return 'stub_{}'.format(xs)
    else:
        name = 'stub_{}_w{}_t{}_{}_{}_s{}_l{}'.\
            format(xs, width, thick, stubshape, pinshape, pinsize, pinlayer)
        return cfg.gds_cellname_cleanup(name)


missing_xs = []
def _makestub(xs_guide=None, width=0, length=2.0, shape=None, pinshape=None,
        pinsize=None, pinlayer=None, cell=None):
    """Create a stub in the logical layers.

    A stub is the stub of a xsection shape around a pin to visualize a connection.
    A pincell is added to the stub to indicate the pin position inside the stub.
    The new stub is added to the stubs dictionary: {name: stubcell}

    Args:
        xs_guide (str): name of xsection
        width (float): stub width
        thick (float): thickness of stub into cell (length)
        shape (str): shape of the stub: 'box' | 'circ' (default = 'box')
        pinshape (string): pinshape used in the stub
        pinsize (float): scaling factor of the pinshape (default = 1)
        cell (Cell): use the provided cell as stub instead of creating a new stub cell

    Returns:
        str: name of the stub
    """
    try:
        xs_logic = cfg.stubmap[xs_guide]
    except:
        # if no stub defined, use the xs as its own stub.
        if xs_guide is None:
            arrow = make_pincell(layer=pinlayer, shape=pinshape, size=pinsize)
            return arrow

        if xs_guide not in cfg.XSdict.keys():
            if xs_guide not in missing_xs:
                missing_xs.append(xs_guide)
                if xs_guide != cfg.default_xs_name:
                    print("Can not make a stub in undefined xsection '{0}'.\n"\
                       "  Possible causes: '{0}' is misspelled or not yet defined.\n"\
                       "  Will use xsection '{2}' instead and continue.\n"
                       "  To define a new xsection:\n"\
                       "      add_xsection(name='{0}')\n"\
                       "  or with layers info and adding a custom stub:\n"\
                       "      add_xsection(name='{0}', layer=1)\n"\
                       "      add_xsection(name='{1}', layer=2)\n"\
                       "      add_stub(xsection='{0}', stub='{1}')".\
                           format(xs_guide, 'stubname', cfg.default_xs_name))
                xs_guide = cfg.default_xs_name

        cfg.stubmap[xs_guide] = xs_guide
        xs_logic = xs_guide

    name = stubname(xs_guide, width, length, shape, pinshape, pinsize, pinlayer)
    if name in stubs.keys():
        return name

    # make a new stub:
    stubshapes = ['box', 'circle']
    arrow = make_pincell(layer=pinlayer, shape=pinshape, size=pinsize)

    if width is None:
        width = 0

    with nd.Cell(name, instantiate=cfg.stub_instantiate, store_pins=False) as C:
        arrow.put(0)

        if cell is not None:
            cell.put()
        else:
            for lay, growx, acc in nd.layeriter(xs_logic):
                if shape is 'circle':
                    outline = geom.circle(radius=0.5*length, N=32)
                else:
                    outline = geom.box(width=width+2*growx, length=length)
                if width != 0:
                    nd.Polygon(layer=lay, points=outline).put(0, 0, 180)

                if shape not in stubshapes:
                    print("Warning: stub shape '{}' not recognized, possible options are {}.".
                        format(shape, stubshapes))
    stubs[name] = C
    return name #name can have change to default_xs


def put_stub(pinname=None, length=None, shape='box', pinshape=None, pinsize=None,
        pinlayer=None, pinshow=True, annotation_layer=None):
    """Add a xsection stub to a pin or multiple pins.

    The xsection and width of the stub is obtained from its pin attributes.
    If no attributes are set, the stub layers are empty and the stubcell
    only contains the pincell.

    Args:
        pinname (str | list of str | None): name(s) of the pin(s) (default = None)
            For pinname=None (default) all chained pins will obtain a stub.
        length (float): length of the stub (thickness into cell)
        shape (string): shape of the stub 'box' | 'cirlce' (default = 'box')
        pinshape (string): pinshape used in the stub
        pinsize (float): scaling factor of the pinshape (default = 1)
        layer (str | int | (int, int)): pin layer
        annotation_layer (str | int | (int, int)): annotation layer

    Returns:
        None
    """
    if isinstance(pinname, str):
        pinname = [pinname]

    if isinstance(pinname, nd.Node):
        pinname = [pinname.name]

    if length is None:
        length = cfg.pin_settings['stub_length']

    if isinstance(shape, nd.Cell):
        stubcell = shape
        shape = shape.cell_name
    else:
        stubcell = None

    if pinname is None or pinname == []:
        #raise ValueError('pinname requries a string name of the pin or list of string names.')
        pinname = [name for name, pin in cfg.cells[-1].pin.items()
            if pin.chain == 1]

    if pinshape is None:
        pinshape = cfg.pin_settings['bb_pin_shape']

    try:
        float(length)
    except:
        raise ValueError("length in 'put_stub' needs to be a number not '{}' of type {}".format(length, type(length)))

    if annotation_layer is None:
        annotation_layer =  cfg.pin_settings['bb_pin_annotation_layer']

    for p in pinname:
        cell = cfg.cells[-1]
        try:
            node = cell.pin[p]
            node.show = pinshow #display pin name in layout and fingerprint
            if stubcell is not None:
                node.show = pinshow #display pin name in layout and fingerprint
                name = "stub_{}".format(shape)
                if name not in stubs.keys():
                    stubs[name] = stubcell
                stubs[name].put(node)
            else:
                if node.xs is None:
                    xs = None
                    width = 0
                else:
                    xs = node.xs
                    try:
                        width = float(node.width)
                    except:
                        width = None

                node.show = True #display pin name in layout and fingerprint
                if xs is None: # pincell only
                    arrow = make_pincell(layer=pinlayer, shape=pinshape, size=pinsize)
                    arrow.put(node)
                else:
                    name = _makestub(xs_guide=xs, width=width, length=length,
                        shape=shape, pinshape=pinshape, pinsize=pinsize,
                        pinlayer=pinlayer)
                    stubs[name].put(node)

            # Move the annotation for readability
            nd.Annotation(layer=annotation_layer, text=p).\
                put(node.move(*cfg.pin_settings['bb_pin_annotation_move']))
        except Exception as E:
            print("warning (core): can't add stub in cell '{}' on pin '{}'. ".\
                format(cell.cell_name, p))
            raise

def set_pdk_pinstyle(pin_settings):
    """Set a custom pin style for a technology.

    Start with the cfg default settings and overrule these for settings provided
    in <pin_settings>.

    Args:
        pin_settings (dict):

    Returns:
        None
    """
    keys = nd.cfg.pin_settings.keys()

    if pin_settings is not None:
        for item, setting in pin_settings.items():
            if item in keys:
                nd.cfg.pin_settings[item] = setting
            else:
                print("Warning: key '{}' not in pin_settings. Available keys: {}".\
                    format(item, keys))

    overrule = nd.cfg.overrule_pdk_pinstyle
    if overrule is not None:
        for item, setting in overrule.items():
            if item in keys:
                nd.cfg.pin_settings[item] = setting
            else:
                print("Warning: key '{}' not in pin_settings. Available keys: {}".\
                    format(item, keys))
        #reset after use:
        nd.cfg.overrule_foundry_pinstyle = None


def export_bb2png(cellcalls, path='', close=True, bbox_pins=True):
    """Create png image for all cells in <cellcalls>.

    Only the Cell objects are used from <cellcalls>, but it does accept as well
    the dict as needed in method 'bb_fingerprint'.

    cellcalls (dict | list of Cells): All the cells {function_call_str: cell_from_call} | list of Cells
    close (bool): close cell immediately after drawing (default = True)
    path (str): path for saving png files
    bbox_pins (bool): add pins to the bbox (default = True)

    Returns:
        None
    """
    export = True

    if isinstance(cellcalls, dict):
        cells = cellcalls.values()
    else:
        cells = cellcalls

    N = len(cells)

    for i, cell in enumerate(cells):
        print('{}/{}'.format(i+1, N))
        bbox_stubs_store = cfg.bbox_stubs
        cfg.bbox_stubs = bbox_pins
        if export:
            output = 'file'
        else:
            output = None
        nd.export_plt(cell, output=output, path=path,
            title=cell.cell_paramsname)
        cfg.bbox_stubs = bbox_stubs_store
        plt.close("all")
    return None


def bb_fingerprint(cellcalls, save=False, filename='fingerprint.json', infolevel=0):
    """Generate a dict with parameter list and pin list of a list off BBs.

    Args:
        cellcalls (dict): All the cells {function_call_str: cell_from_call}

    Returns:
        dict: Dictionary with building block info, Structure is as follows::

            {'<cellname>':
                {'pins':
                    { '<pinname>':
                        {'width': <value>,
                         'xs'   : <value>
                         'xya'  : (<x>, <y>, <a>)
                         'show' : <values>
                        }
                        , ...
                    },
                 'parameters':
                     {'<keywword>': <value>}
                     , ...
                 'paramsname'   : <value>,
                 'basename'     : <value>,
                 'groupname'    : <value>,
                 'function'     : <value>,
                 'length'       : <value>,
                 'width'        : <value>
                }
            }
    """
    fingerprint = {}
    cellcallsorted = sorted(cellcalls.items()) #not effective for Python with scrambled dict.
    for callstr, cell in cellcallsorted:
        try:
            groupname = cell.groupname
        except:
            groupname = ''

        pinsetting = {}
        if infolevel > 0:
            print(cell.cell_basename)
        for name, pin in sorted(cell.pin.items()):
            if name != 'org':
                pinsetting[name] = {
                    'width': pin.width,
                    'xsection': pin.xs,
                    'pointer': [float('{:.5f}'.format(xya)) for xya in pin.xya()],
                    'show': pin.show,
                    'type': pin.type,
                    'remark': pin.remark
                }
                if infolevel > 0:
                    print(name,  pinsetting[name])

        paramsetting = {}
        try: # 'parameters' may be empty
            for p in cell.parameters:
                paramsetting[p] = cell.parameters[p]
        except:
            pass

        #bbox
        try:
            length = cell.length
            width = cell.width
        except:
            length = 'None'
            width = 'None'
        fingerprint[cell.cell_name] = {
            'pins': pinsetting,
            'parameters': paramsetting,
            'paramsname': cell.cell_paramsname,
            'basename': cell.cell_basename,
            'groupname': groupname,
            'function': callstr,
            'length': length,
            'width': width,
            'connect': (cell.default_in, cell.default_out)}

    if save:
        with open(filename, 'w') as fp:
            json.dump(fingerprint, fp, indent=2, sort_keys=True)
        print("Saved fingerprint to {}".format(filename))
    if infolevel > 1:
        pprint(fingerprint)
    return fingerprint


def validate_black_to_white_mapping(black2whiteMap, allBBcells, infolevel=0):
    """Validate if all white and black boxes are mapped.

    By increasing the infolevel more information will be displayed on where
    black <-> white mappings are missing.

    Args:
        black2whiteMap (dict): {blackbox-basename: whitebox-function}
        allBBcells (list of Cell): all black box cells.
        infolevel: amount of feedback to stdout (default = 0)

    Returns:
        bool: True if mapping is successful.
    """

    if infolevel > 0:
        print('Validate black <-> white box mappings:')

    okay = True
    mesg = []
    #print("\nStatus of white-box implementation:")
    mappings = ""
    header = "   {0:30}{1}".format("BASENAME", "FUNCTION")
    black_in_map = list(black2whiteMap.keys())
    basenames = {cell.cell_basename: cell for cell in allBBcells}
    black_count = len(basenames)
    white_found = 0
    for basename, cell in sorted(basenames.items()):
        func = ""
        defined = '?'
        if basename in black_in_map:
            defined = '*'
            white_found += 1
            func = black2whiteMap[basename]
        mappings += "{0}  {1:30}{2}\n".format(defined, basename, func)

    black_orphan_count = black_count - white_found
    text = ("Found {} black orphans".format(black_orphan_count))
    print("{}".format(text))
    if black_orphan_count > 0:
        okay = False
        mesg.append(text)
        if infolevel > 0:
            print(header)
            print(mappings)

    white_orphans = set(black2whiteMap) - set(basenames.keys())
    white_orphan_count = len(white_orphans)
    text = "Found {} white orphans".format(white_orphan_count)
    print("{}".format(text))
    if white_orphan_count > 0:
        okay = False
        mesg.append("{}{}".format(text, '.'))
        if infolevel > 0:
            print(header)
            for name in white_orphans:
                print("?  {:30}{}".format(name, black2whiteMap[name]))

    if not okay:
        print("\nError(s) found in white <-> black box mappings.")
        for m in mesg:
            print("- {}{}".format(m, '.'))
    else:
        print("\nSuccesful black <-> white mapping.")


bbox_pinnames = [
    'lb', 'lc', 'lt',
    'tl', 'tc', 'tr',
    'rt', 'rc', 'rb',
    'br', 'bc', 'bl',
    'cc']

def put_boundingbox(pinname, length, width, raise_pins=True, outline=True,
        align='lc', name=True, params=True, move=(0, 0, 0)):
    """Create bounding box (bbox) cell inside the active cell.

    This function places a bbox cell and raises by default the bbox pins into
    the active cell. The bbox displays a bbox outline (can be switched off).
    By default it also adds the active cellname and parameters.

    Args:
        pin (str): pin to place bbox on (center left of bbox)
        length (float): length of the bbox
        with (float): width of the bbox
        raise_pins (bool): raise bbox pins into active cell (default = True)
        outline (bool): draw bbox outline (default = True)
        align (str): align the bbox on the specified bbox pin <pinname> (default = 'lc')
        name (bool): display the (active) cell name in the bbox. (default = True)
        params (bool): add parameter annotation to the bbox
        move (tuple): move the bbox placement by (float, float, float)

    Returns:
        None
    """
    cell = cfg.cells[-1]
    _paramsname = cell.cell_paramsname
    _parameters = cell.parameters
    with nd.Cell(name='bbox', instantiate=False) as bbox:
        outline = geom.box(length, width)
        nd.Polygon(layer='bbox', points=outline).put(0)

        stbs = cfg.bbox_stubs and outline
        nd.Pin('lb', type='bbox', show=stbs).put(0, -0.5*width, 180)
        nd.Pin('lc', type='bbox', show=stbs).put(0, 0, 180)
        nd.Pin('lt', type='bbox', show=stbs).put(0, 0.5*width, 180)

        nd.Pin('tl', type='bbox', show=stbs).put(0, 0.5*width, 90)
        nd.Pin('tc', type='bbox', show=stbs).put(0.5*length, 0.5*width, 90)
        nd.Pin('tr', type='bbox', show=stbs).put(length, 0.5*width, 90)

        nd.Pin('rt', type='bbox', show=stbs).put(length, 0.5*width, 0)
        nd.Pin('rc', type='bbox', show=stbs).put(length, 0, 0)
        nd.Pin('rb', type='bbox', show=stbs).put(length, -0.5*width, 0)

        nd.Pin('br', type='bbox', show=stbs).put(length, -0.5*width, -90)
        nd.Pin('bc', type='bbox', show=stbs).put(0.5*length, -0.5*width, -90)
        nd.Pin('bl', type='bbox', show=stbs).put(0, -0.5*width, -90)

        nd.Pin('cc', type='bbox', show=cfg.bbox_stubs).put(0.5*length, 0, 0)

        if stbs:
            nd.put_stub(['lb', 'tl', 'rt', 'br'],
                pinshape='arrow_righthalf',
                pinsize=cfg.pin_settings['bbox_pin_size'],
                pinlayer=cfg.pin_settings['bbox_pin_layer'],
                annotation_layer='bbox_pin_text')
            nd.put_stub(['lt', 'tr', 'rb', 'bl'],
                pinshape='arrow_lefthalf',
                pinsize=cfg.pin_settings['bbox_pin_size'],
                pinlayer=cfg.pin_settings['bbox_pin_layer'],
                annotation_layer='bbox_pin_text')
            nd.put_stub(['lc', 'tc', 'rc', 'bc'],
                pinshape='arrow_full',
                pinsize=cfg.pin_settings['bbox_pin_size'],
                pinlayer=cfg.pin_settings['bbox_pin_layer'],
                annotation_layer='bbox_pin_text')

        if name:
            cellname(cellname=_paramsname, length=length, width=width, align='cc').\
                put(bbox.pin['cc'])
        if params:
            parameters(parameters=_parameters).put(bbox.pin['cc'])

    # options to align the bbox in this cell w.r.t. to its pins:
    align_shift = {
        'lb': (0, 0.5*width),
        'lc': (0, 0),
        'lt': (0, -0.5*width),
        'rb': (-length, 0.5*width),
        'rc': (-length, 0),
        'rt': (-length, -0.5*width),
        'bc': (-0.5*length, 0.5*width),
        'cc': (-0.5*length, 0),
        'tc': (-0.5*length, -0.5*width)
    }
    dx = align_shift[align][0]
    dy = align_shift[align][1]

    box = bbox.put(cfg.cells[-1].pin[pinname].move(dx, dy).move(*move))
    cell.length = length
    cell.width = width

    if raise_pins:
        box.raise_pins(bbox_pinnames, show=False)

    return None


#TODO: add flop
def load_gdsBB(gdsfilename, cellname, pinfilename=None, newcellname=None,
        layermap=None, cellmap=None, flip=False, scale=1.0, stubs=True, native=True,
        bbox=False, bboxbuf=0, prefix='', instantiate=None):
    """Load a gds cell and the corresponding pin info from file.

    This function uses method 'load_gds' for loading the gds file
    and adds pins and stubs according to the pinfile.
    The combination creates building block cell with connectivity.
    Without pinfile method 'load_gds" can be used instead.

    If stubs == True then both the stubs and the loaded GDS are instantiated
    inside the returned Cell.
    This to ensure that for this cell the gds and stubs remaind aligned under flipping.

    Args:
        gdsfilename (str): gds filename to import from
        cellname (str): cellname to import
        pinfilename (str): optional name of csv file containing pin info for stubs.
        newcellname (str): new cell name
        layermap (dict): mapping of layer {number_in: number_out}
        flip (bool): mirror the gds (default = False)
        scale (float): scaling factor (default = 1.0). Use with great care
            as scaled building blocks will not make much sense from a functional
            prespective.
            pin positions will be scaled along with the gds, however,
            xs information, like 'width' will *not* be scaled
        stubs (bool): place stubs (default = True)
        native (bool):
        bbox (bool):
        bboxbuf (float):

    Returns:
        Cell: Nazca Cell with the loaded gds and (if provided) pins and stubs.
    """
    if stubs:
        instantiate = True
        bboxgds = False
    else:
        instantiate = False
        bboxgds = bbox
    #TODO: add instantiate override

    with nd.Cell(cellname+'_stubs', instantiate=instantiate) as C:
        #TODO: if instantiate is true the instance and gds file may get the same name
        #      causing a topology error in GDS.
        nd.load_gds(
            filename=gdsfilename,
            cellname=cellname,
            newcellname=newcellname,
            layermap=layermap,
            cellmap=cellmap,
            scale=scale,
            native=native,
            bbox=bboxgds,
            bboxbuf=bboxbuf,
            prefix=prefix).\
                put(0, flip=flip)
        if flip:
            s = -1.0
        else:
            s = 1.0
        if pinfilename is not None:
            df = pd.read_csv(pinfilename, delimiter= ',', skiprows=1, header=None,
                names = ['io', 'x', 'y', 'a', 'xs', 'w'])
            df.fillna(0.0, inplace=True)
            for row in df.itertuples():
                i, io, x, y, a, xs, w = row
                try:
                    xs = xs.strip()
                    if xs.upper() == 'NONE':
                        xs = None
                except:
                    xs = None
                if not isinstance(w, float):
                    w = 0.0

                nd.Pin(name=io, xs=xs, width=w).put(scale*x, s*scale*y, a)
                if stubs: #and xs is not None :
                    put_stub(io)
        if instantiate and bbox:
            C.autobbox = True
            C.bboxbuf = bboxbuf
    return C


def _PIL2array(img):
    return np.array(img.getdata(), np.bool).reshape(img.size[1], img.size[0])


def image(name, layer=1, size=256, pixelsize=1, threshold=0.5,
        cellname=None, invert=False, align='cc', box_layer=None, box_buf=0):
    """Read an image file and return a nazca cell with the image.

    Image format can be png, jpg, gif, bpm, eps and others, as supported by Pillow.
    Note that the output resolution (size) does not exceed the image resolution.
    Increase <pixelsize> in this case to obtain a larger logo in gds.

    A rectangular box can be added around the logo by providing a <box_layer>.
    This box can be enlarged beyond the original image size by setting <box_buf> > 0.

    Args:
        name (str): name of the image file
        layer (int): layer number that the image will be written to (default 1)
        size (int): maximum bounding box size in pixels (default 256)
        pixelsize (float): pixel size in micron (default 1)
        threshold (float): black/white threshold (default 0.5)
        cellname (str): Nazca cell name (default image filename)
        invert (bool): flag to invert black & white (default False)
        align (str): two character string for image alignment (default 'cc')
            allowed:
            lt, ct, rt,
            lc, cc, rc,
            lb, cb, rb

        box_layer (str | int | tuple): layer reference to generate a rectangular
            box behind the text, e.g. for tiling exclusion areas (NOFILL)
            (default = None)
        box_buf (float): extra buffer for the box_layer in um

    Returns:
        Cell: cell with image

    Examples:
        Load a logo in a cell and put and/or export it to gds::

            import nazca as nd

            logo = nd.image('mylogo.png', align='lb') # left/bottom alignment
            logo.put(0)
            nd.export_gds()

        or::

            import nazca as nd

            nd.export_gds(logo, filename='mylogo.gds')
    """
    if cellname is None:
        cellname = os.path.basename(name)
    p = pixelsize
    threshold = int(threshold * 256)
    a = {'lb', 'cb', 'rb', 'lc', 'cc', 'rc', 'lt', 'ct', 'rt'}
    if align not in a:
        print("Invalid alignment specification '{}' for image '{}'.".format(align, name))
        print("Allowed values are {}.".format(a))
        print("Using default value 'cc'")
        align = 'cc'
    halign = {'l': 0, 'c': -0.5, 'r': -1}
    valign = {'b': 0, 'c': -0.5, 't': -1}

    im = Image.open(name)
    gray = im.convert('L')
    # resize keep aspect, only if smaller
    gray.thumbnail((size, size), Image.ANTIALIAS)
    bw = gray.point(lambda x: 0 if x<threshold else 255, '1')
    pix = _PIL2array(bw)
    width, height = bw.size
    width_tot = width*p + 2*box_buf
    height_tot = height*p + 2*box_buf
    print('Generating {}x{} pixels image of {:.0f}x{:.0f} um2, edge is {} um.'.\
         format(width, height, width*p, height*p, box_buf))
    h0 = halign[align[0]] * width_tot
    v0 = valign[align[1]] * height_tot

    with nd.Cell(cellname) as C:
        for line in range(height):
            x1 = x0 = lb = lw = 0
            y0 = (height - line) * p
            for pixel in pix[line]:
                if pixel == invert:
                    lb += 1
                    if lw > 0:
                        x0 += lw * p
                        lw = 0
                else:
                    lw += 1
                    if lb > 0:
                        x1 = x0 + lb * p
                        xy = [(x0, y0), (x1, y0), (x1, y0-p), (x0, y0-p)]
                        nd.Polygon(layer=layer, points=xy).\
                            put(h0+box_buf, v0+box_buf, 0)
                        x0 = x1
                        lb = 0
            if lb > 0:
                x1 = x0 + lb * p
                xy = [(x0, y0), (x1, y0), (x1, y0-p), (x0, y0-p)]
                nd.Polygon(layer=layer, points=xy).\
                    put(h0+box_buf, v0+box_buf, 0)

        if box_layer is not None:
            nd.Polygon(layer=box_layer, points=[(0, 0), (0, height_tot),
                (width_tot, height_tot), (width_tot, 0)]).put(h0, v0, 0)
    return C


# =============================================================================
# predefined example cell for documentation examples
# =============================================================================
def example_cell():
    """Create a Nazca example cell.

    Returns:
        cell
    """
    with nd.Cell('example_cell') as example:
        s = nd.strt(length=20).put(0)
        b = nd.bend(angle=45).put()
        nd.Pin('a0', pin=s.pin['a0']).put()
        nd.Pin('b0', pin=b.pin['b0']).put()
    return example

dir_path = os.path.dirname(os.path.abspath(__file__))

def nazca_logo(layers={'ring':1, 'box':2, 'bird':3}, cell_name='nazca_logo',
        scale=1.0):
    """Return the nazca logo as a cell.

    Args:
        layers (dict): gds layers of the logo. default = {'ring':1, 'box':2, 'bird':3}
        cell_name (str): default = 'nazca_logo'
        scale (float): scaling factor (default = 1)
    """
    return load_gdsBB(
	    gdsfilename=os.path.join(dir_path, 'nazca_logo.gds'),
		cellname='nazca',
		newcellname=cell_name,
        scale=1.0,
	    layermap={1:layers['ring'], 2:layers['box'], 3:layers['bird']})

