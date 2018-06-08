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
import operator
from functools import wraps
import hashlib
import matplotlib.pyplot as plt

import pandas as pd
import nazca as nd
from nazca import cfg
from nazca.util import parameters_to_string
import nazca.geometries as geom


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
                    print(
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
""")

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
                            #avoid things like 0.0000000001 in names.
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

            # add cfg version to avoid importing template_core for hashme attributes
            cfg.hash_id = _hash_id
            cfg.hash_params = _hash_params
            cfg.hash_basename = cfg.gds_cellname_cleanup(name)         # base
            cfg.hash_paramsname = cfg.gds_cellname_cleanup(name_long)  # base + params
            cfg.hash_name = cfg.gds_cellname_cleanup(_hash_name)       # base + params + hash

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
            return cell_new
        return wrapper
    return decorator


def rangecheck(allowed_values):
    """Check if parameter values are in range.

    Args:
        allowed values (dict): {str: tuple} as {var_name: (low, var_value, high)}

    Raises:
        ValueError if out of range.

    Returns:
        None

    Example:

        check 0 <= a <= 10 and -5 <= b <= 5::

            def func(a, b):
                nd.rangecheck({
                    'a', (0, a, 10),
                    'b', (-5, a, 5)
                })
    """
    err = ''
    for var, (low, value, high) in allowed_values.items():
        if low <= value <= high:
            pass
        elif value < low:
            err += "{0}={3} too small. Allowed values: {1}<={0}<={2}\n".\
                 format(var, low, high, value)
        elif value > high:
            err += "{0}={3} too large. Allowed values: {1}<={0}<={2}\n".\
                 format(var, low, high, value)
    if err != '':
        raise ValueError(err)
    return err

#==============================================================================
#
#==============================================================================
stubs = dict() # All stubs. {stub_name: stub_cell_obj}
stubmap = dict() # Map xs to its stub_xs. {xs_name: stub_xs_name}


class Functional_group():
    """
    Class to group building blocks syntactically together.

    e.g. a bb_splitters group with various splitter components:

        bb_splitters.mmi1x2

        bb_splitters.mmi2x2

        bb_splitters.mmi3x3
    """

    def __init__(self, name=''):
        self.name = name
BB_Group = Functional_group


def parameters(parameters=None):
    """Put parameter list as annotation in a building block.

    Args:
        parameters (dict): {parameter_name: value}

    Returns:
        None
    """

    #store parameters as a dict in the cell for PDK checks

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
    """Create text of the cellname and position current pin to place it.

    Args:
        cellname (str): name of the cell
        length (float): length available for the BB name in um
        align (str): text alignment (see nazca.text, default = 'lc')

    Returns:
        Cell: text with cellname formated
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
arrow = None

#@hashme('arrow', 'layer', 'pinshape', 'pinshape')
# Do not @hashme this function (autoname) this function
# as parameter values are finalized inside the function.
def make_pincell(layer=None, pinshape=None, pinsize=None):
    """Create a cell with an arrow shape to indicate pin positions.

    Args:
        layer (float): layer number to place the arrow
        arrowtype (str): geometry type of the arrow: full | lefthalf | righthalf

    Returns:
        None
    """
    pinshape_default = 'arrow_full'
    if layer is None:
        layer = 'bb_pin'

    if layer not in cfg.layerdict.keys():
        try:
           layer = cfg.default_layer_Annotation
           nd.add_layer(name=layer, layer=layer)
        except:
            pass

    if pinshape is None:
        pinshape = pinshape_default

    if pinsize is None:
        pinsize = cfg.bb_pin_size

    if pinshape not in cfg.pinshapes.keys():
        print("Warning: arrowtype '{}' not recognized.".format(pinshape))
        print("Available options are {}.".format(cfg.pinshapes))
        print("I will use arrowtype = '{}'.".format(pinshape_default))
        pinshape = pinshape_default

    name = "{}_{}_{}_{}".format('arrow', layer, pinshape, pinsize)
    if name in cfg.cellnames:
        return cfg.cellnames[name]
    with nd.Cell(name) as arrow:
        outline = [(x*pinsize, y*pinsize) for x, y in cfg.pinshapes[pinshape]]
        nd.Polygon(layer=layer, points=outline).put(0)
    return arrow


def stubname(xs, width, thick, shape=None, size=None):
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
        name = 'stub_{}_w{}_t{}_{}'.format(xs, width, thick, shape, size)
        return cfg.gds_cellname_cleanup(name)


missing_xs = []
def makestub(xs_guide, width, length=2.0, shape=None, pinshape=None, pinsize=None, cell=None):
    """Factory for creating stubs in the logical layers.

    Args:
        xs_guide (str): name of xsection
        width (float): stub width
        thick (float): thickness of stub into cell (length)
        shape (str): shape of the stub: circ | rect (default = 'rect')

    Returns:
        None
    """
    arrow = make_pincell(pinshape=pinshape, pinsize=pinsize)
    try:
        xs_logic = cfg.stubmap[xs_guide]
    except:
        # if no stub defined, use the xs as its own stub.
        if xs_guide is None:
            return arrow

        if xs_guide not in cfg.XSdict.keys():
            if xs_guide not in missing_xs:
                missing_xs.append(xs_guide)
                print("Can not make a stub in an undefined xsection '{0}'.\n"\
                   "  Possible causes: '{0}' is misspelled or not yet defined.\n"\
                   "  Will use xsection '{2}' instead and continue.\n"
                   "  To define a new xsection:\n"\
                   "  add_xsection(name='{0}')\n"\
                   "  or with layers info and adding a custom stub:\n"\
                   "  add_xsection(name='{0}', layer=1)\n"\
                   "  add_xsection(name='{1}', layer=2)\n"\
                   "  add_stub(xsection='{0}', stub='{1}')".\
                       format(xs_guide, 'stubname', cfg.default_xs_name))
                xs_guide = cfg.default_xs_name

        cfg.stubmap[xs_guide] = xs_guide
        xs_logic = xs_guide

    stubshapes = ['box', 'circle']
    name = stubname(xs_guide, width, length, shape, pinsize)

    if width is None:
        width = 0

    with nd.Cell(name) as C:
        arrow.put(0)

        if cell is not None:
            cell.put()
        else:
            for lay, growx, acc in nd.layeriter(xs_logic):
                if shape is 'circle':
                    outline = geom.circle(radius=0.5*length, N=32)
                else:
                    outline = geom.box(width=width+2*growx, length=length)
                nd.Polygon(layer=lay, points=outline).put(0, 0, 180)

                if shape not in stubshapes:
                    print("Warning: stub shape '{}' not recognized, possible options are {}.".
                        format(shape, stubshapes))
    stubs[name] = C
    return name #name can have change to default_xs


def put_stub(pinname=None, length=2.0, shape='box', pinshape=None, pinsize=1.0,
        layer=None, annotation_layer=None):
    """Add a xsection stub to a pin.

    The xsection and width of the stub is obtained from the pin attributes.

    Args:
        pinname (str | list of str): name(s) of the pin(s) (default = None)
            If pinname is None all the chained pins will obtain a stub.
        length (float): length of the stub (thickness into cell)

    Returns:
        None
    """

    if isinstance(pinname, str):
        pinname = [pinname]

    if isinstance(pinname, nd.Node):
        pinname = [pinname.name]

    if isinstance(shape, nd.Cell):
        stubcell = shape
        shape = shape.cell_name
    else:
        stubcell = None

    if pinname is None or pinname == []:
        #raise ValueError('pinname requries a string name of the pin or list of string names.')
        pinname = [name for name, pin in cfg.cells[-1].pin.items() if pin.chain == 1]

    if pinshape is None:
        pinshape = cfg.bb_pin_shape

    if annotation_layer is None:
        annotation_layer = 'bb_pin_text'

    for p in pinname:
        cell = cfg.cells[-1]
        try:
            node = cell.pin[p]
            if node.xs is None:
                xs = None
                width = 0
            else:
                xs = node.xs
                width = node.width

            node.show = True #display pin name in layout
            if xs is None:
                arrow = make_pincell(pinshape=pinshape, pinsize=pinsize, layer=layer)
                arrow.put(node)
            else:
                name = stubname(xs, width, length, shape, pinsize)
                if name not in stubs.keys():
                    name = makestub(xs_guide=xs, width=width, length=length,
                        shape=shape, pinshape=pinshape, cell=stubcell)
                stubs[name].put(node)

                # Move the annotation a bit so it is inside the arrow and the labels of
                # connected pins can be distinguished.
                #TODO: 0.3 should be parametrized.

            nd.Annotation(layer=annotation_layer, text=p).put(node.move(-0.3))
        except:
            print("warning (core): can't add stub in cell '{}': pin '{}' not existing".\
                format(cell.cell_name, p))


def make_pdk_lists(functioncalls):
    """Create a list of Cell and a dictionary of function calls strings.

    Returns:
        list, Dict: list of Cell, dict of {Cell, functioncall_ string}
    """
    BBcells = [] # list with all BBs to show in the manual
    BBdict = {} # dictionary for bbfingerprint.
    for F in functioncalls:
        BBcells.append(eval(F))
        #BBfunc_truncs.append(F[:F.index('(')])
        BBdict[BBcells[-1]] = F
    return BBcells, BBdict


def bb_fingerprint(cellcall, export_png=False, path='', \
        bbox_stubs=True, close_plt=True):
    """Generate a dict with parameter list and pin list of a list off BBs.

    Args:
        functioncall (list of str): function calls to fingerprint
        export_png (bool): generate png output and save in <path>
        path (str): path for saveing png files
        bbox_stubs (bool): add stubs to the bbox (default = True)
            Only relevant with export_png True

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
    BBsettings = {}
    cellcallsorted = sorted(cellcall.items(), key=operator.itemgetter(1))
    #print('sorted:', cellcallsorted)
    for cell, callstr in cellcallsorted:
        try:
            groupname = cell.groupname
        except:
            groupname = ''

        if export_png:
            bbox_stubs_store = cfg.bbox_stubs
            cfg.bbox_stubs = bbox_stubs
            nd.export_plt(cell, output='file', path=path, title=cell.cell_paramsname)
            cfg.bbox_stubs = bbox_stubs_store

        pinsetting = {}
        for name, pin in sorted(cell.pin.items()):
            if name != 'org':
                pinsetting[name] = {
                    'width':pin.width,
                    'xsection':pin.xs,
                    'pointer':pin.xya(),
                    'show':pin.show,
                    'type':pin.type
                }

        paramsetting = {}
        try: # parameters may be empty
            for p in cell.parameters:
                paramsetting[p] = cell.parameters[p]
        except:
            pass
        try:
            length = cell.length
            width = cell.width
        except:
            length = 'None'
            width = 'None'
        BBsettings[cell.cell_name] = {
            'pins': pinsetting,
            'parameters': paramsetting,
            'paramsname': cell.cell_paramsname,
            'basename': cell.cell_basename,
            'groupname': groupname,
            'function': callstr,
            'length': length,
            'width': width}

        plt.close()
    return BBsettings


bbox_pinnames = [
    'lb', 'lc', 'lt',
    'tl', 'tc', 'tr',
    'rt', 'rc', 'rb',
    'br', 'bc', 'bl',
    'cc']

def put_boundingbox(pinname, length, width, raise_pins=True, outline=True,
        align='lc', name=True, params=True, move=(0, 0, 0)):
    """Create  bounding box (bbox) cell.

    Args:
        pin (str): pin to place bbox on (center left of bbox)
        length (float): length of the bbox
        with (float): width of the bbox
        raise_pins (bool): raise bbox pins into active cell (default = True)
        outline (bool): draw bbox outline (default = True)
        align (str): align the bbox on pin <pinname> (default = 'lc')

    Return:
        None
    """
    cell = cfg.cells[-1]
    _paramsname = cell.cell_paramsname
    _parameters = cell.parameters
    with nd.Cell(name='bbox', instantiate=False) as bbox:
        outline = geom.box(length, width)
        nd.Polygon(layer='bbox', points=outline).put(0)

        stbs = cfg.bbox_stubs
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
                pinshape='arrow_righthalf', pinsize=cfg.bbox_pin_size,
                layer='bbox_pin', annotation_layer='bbox_pin_text')
            nd.put_stub(['lt', 'tr', 'rb', 'bl'],
                pinshape='arrow_lefthalf', pinsize=cfg.bbox_pin_size,
                layer='bbox_pin', annotation_layer='bbox_pin_text')
            nd.put_stub(['lc', 'tc', 'rc', 'bc'],
                pinshape='arrow_full', pinsize=cfg.bbox_pin_size,
                layer='bbox_pin', annotation_layer='bbox_pin_text')

        if name:
            cellname(cellname=_paramsname, length=length, width=width, align='cc').\
                put(bbox.pin['cc'])
        if params:
            parameters(parameters=_parameters).put(bbox.pin['cc'])

    # options to align the bbox in this cell w.r.t. to its pins:
    align_options = [ 'lb', 'lc', 'lt', 'rb', 'rc', 'rt', 'cc']

    align_translate = [
        (0, 0.5*width), (0, 0), (0, -0.5*width),
        (-length, 0.5*width), (-length, 0), (-length, -0.5*width),
        (-0.5*length, 0)
    ]

    dx = align_translate[align_options.index(align)][0]
    dy = align_translate[align_options.index(align)][1]

    box = bbox.put(cfg.cells[-1].pin[pinname].move(dx, dy).move(*move))
    cell.length = length
    cell.width = width

    if raise_pins:
        box.raise_pins(bbox_pinnames, show=False)

    return None


def load_gdsBB(gdsfilename, cellname, pinfilename=None, newcellname=None,
        layermap=None, flip=False, scale=1, stubs=True):
    """Load a GDS building block and its corresponding pin file.

    If stubs == True then both the stubs and the loaded GDS are instantiated
    inside the Cell returned by load _gdsBB to keep them in sync under flipping.

    Args:
        gdsfilename (str): gds filename to import from
        cellname (str): cellname to import
        pinfilename (str): optional name of csv file containing pin info for stubs.
        newcellname (str): new cell name
        layermap (dict): mapping of layer {number_in: number_out}
        flip (bool): mirror the gds (default = False)
        scale (float): scaling factor (default = 1), not implemented
        stubs (bool): place stubs (default = True)

    Returns:
        Cell
    """
    if stubs:
        instantiate = True
    else:
        instantiate = False

    with nd.Cell(cellname+'_stubs', instantiate=instantiate) as C:
        #TODO: if this is true the instance and gds file may get the same name
        #      causing a topology error in Klayout
        nd.load_gds(
            filename=gdsfilename,
            cellname=cellname,
            newcellname=newcellname,
            layermap=layermap,
            scale=scale
            ).put(0, flip=flip)
        #TODO: if the GDS gets flipped as original the pins need to follow
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

                nd.Pin(name=io, xs=xs, width=w).put(x, y, a)
                if stubs and xs is not None :
                    put_stub(io)
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

def nazca_logo(layers={'ring':1, 'box':2, 'bird':3}, cell_name='nazca_logo'):
	return load_gdsBB(
		gdsfilename=os.path.join(dir_path, 'nazca_logo.gds'),
		cellname='nazca',
		newcellname=cell_name,
		layermap={1:layers['ring'], 2:layers['box'], 3:layers['bird']})

