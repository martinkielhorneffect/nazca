#!/usr/bin/env python3
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
# @author: Ronald Broeke (c) 2016-2017
# @email: ronald.broeke@brightphotonics.eu
#-----------------------------------------------------------------------
# -*- coding: utf-8 -*-
"""
This module defines
1- layers in DataFrame 'layer_table'
2- xsections in DataFrame 'xsection_table.
3- xsection layers in DataFrame 'xsection_layer table', adding layers to each 'xs'
4- The 'layer_table' and 'xsection_layer table' are joined using the
   'layer2xsection table'.

The result is stored in a dictionary containing the xsection-layer information
per 'xs' needed in layout: 'xsection_layers'

Layers and xsection can be added in two ways
1 - Loading from csv file: load_layers, load_xs
2 - Functions add_layer, add_xsection
When adding layers to a 'xs' the xs_map is built automatically.

It is also possible to first load layers and xs, and then extend them via
'add_layer' and/or 'add_layer2xsection'.

Ronald Broeke 2017(c)
"""


import pandas as pd
import matplotlib as mpl
from matplotlib.colors import rgb2hex
from . import cfg
from . import xsection
from pprint import pprint

# keep track of layers used by the designer that have not been defined:
unknown_layers = set()
unknown_xsections = set()
cfg.layerset = set()

cfg.xsection_table = pd.DataFrame()
xsection_table_attr = ['xsection', 'xsection_foundry', 'origin', 'stub']

#layer-table:
cfg.layer_table = pd.DataFrame()
layer_table_attr_csv  = ['layer', 'datatype', 'layer_name', 'layer_name_foundry', 'accuracy', 'origin',       'remark']
#rename attributes for the xs-layer join to avoid name clashes:
layer_table_attr_join = ['layer', 'datatype', 'layer_name', 'layer_name_foundry', 'accuracy', 'origin_layer', 'remark']


#layercolor table:
cfg.colors = pd.DataFrame()


#xsection-layer-map table, mapping layers to xsections:
cfg.xsection_layer_map = pd.DataFrame()
xsection_layer_map_attr_csv  = ['layer', 'datatype', 'xsection', 'growx', 'growy', 'origin']
#rename attributes for the xs-layer join to avoid name clashes:
xsection_layer_map_attr_join = ['layer', 'datatype', 'xsection', 'growx', 'growy', 'origin_xs']

xsection_table_attr

#final set of export attributes for joined xsection-layers table:
mask_layers_attr = ['xsection', 'layer', 'datatype', 'growx', 'growy', 'accuracy']


def add_layer(
    name=None,
    layer=None,
    accuracy=None,
    fab_name=None,
    origin=None,
    remark=None,
    # next color info:
    frame_color=None,
    fill_color=None,
    frame_brightness=None,
    fill_brightness=None,
    dither_pattern=None,
    valid=None,
    visible=None,
    transparent=None,
    width=None,
    marked=None,
    animation=None,
    alpha=None,
    unknown=False):

    """Create a new mask layer.

    Args:
        name (str): layer name
        layer (int | tuple): layer number or (layer, datatype)

    Returns:
        DataFrame: Table of mask layers
    """

    #if layer is None:
    if unknown: #redirected from get_layer
        _layer, _datatype = layer
    else:
        _layer, _datatype = get_layer(layer)
    #layer = (cfg.default_xs['layer'][0], 0)
    #if datatype is None:
    #    datatype = 0
    if accuracy is None:
        accuracy = cfg.default_xs['accuracy'][0]
    if origin is None:
        origin = cfg.default_xs['origin'][0]
    if name is None:
        name = 'autoname_{}_{}'.format(_layer, _datatype)

    #rewrite layer:
    if not cfg.layer_table.empty and name in cfg.layer_table['layer_name'].values:
        #print('Warning: overwriting already existing layer name \'{}\'.'.format(name))
        select = cfg.layer_table ['layer_name']==str(name)
        columns = ['layer', 'datatype', 'accuracy']
        cfg.layer_table.loc[select, columns] = _layer, _datatype, accuracy
    #new layer:
    else:
        newlayer = {
            'layer': [_layer],
            'datatype': [_datatype],
            'accuracy': [accuracy],
            'layer_name': str(name),
            'layer_name_foundry': fab_name,
            'origin': [origin],
            'remark': remark}
        df = pd.DataFrame(newlayer)
        cfg.layer_table = pd.concat([cfg.layer_table, df], ignore_index=True)
        cfg.layer_table.index.name = 'id'

        #create a dict for speeding up get_layer
        tab = cfg.layer_table.set_index('layer_name')
        tab2 = tab.loc[:, ['layer', 'datatype']].stack().unstack(0)
        cfg.layerdict = tab2.to_dict(orient='list')

    set_layercolor(layer,
        frame_color, fill_color, frame_brightness, fill_brightness,
        dither_pattern, valid, visible, transparent, width, marked, animation,
        alpha)

    merge_xsection_layers_with_layers()
    return cfg.layer_table


def add_XSdict(name, name_foundry=None, origin='user', stub=None):
    """Create a new Xsection object named <name>.

    If a Xsection with <name> already exists, the existing Xsection is returned.

    Args:
        name (str): xsection name.
        name_foundry (str): optional xsection name used by the foundry
        origin (str): source/creator of the layer, e.g. 'user', 'foundry' or any string.
        stub: xsection name to be used for the stub of this xsection

    Returns:
        Xsection: Xsection object with name <name>
    """

    if name not in cfg.XSdict.keys():
        cfg.XSdict[name] = xsection.Xsection(name=name)

        # add new info to xsection_layers:
        D = {
            'xsection': name,
            'xsection_foundry': name_foundry,
            'origin': origin,
            'stub': stub
        }
        cfg.xsection_table = cfg.xsection_table.append(D, ignore_index=True)

        #update mask_layers:
        merge_xsection_layers_with_layers()
    return cfg.XSdict[name]

add_xsection = add_XSdict


def add_layer2xsection(xsection='name', layer=None, growx=None,
        growy=None, accuracy=0.001, fab_name=None, origin=None, remark=None):
    """Add a layer to the Xsection object with name <xsection>.

    If <xsection> does not exist yet it will be created.
    The layer entry is described by (layer, datatype).
    If the (layer, datatype) does not exist yet it will be created.
    The keyword 'accuracy' is only processed in case of a new layer.

    Args:
        xsection (str): xsection name for use in nazca
        layer (int | tuple): layer number or (layer, datatype)
        growx (float): growth in x-direction of polygon in the (layer, datatype)
            (default = 0.0)
        growy (float): growth in y-direction of polygon in the (layer, datatype)
            (default = 0.0)
        accurary (float): accuracy of the grid in um of the (layer, datatype)
            (default = 0.001)
        fab_name (str): technology xsection name, e.g. the foundry name for the layer
        origin (str): source/creator of the layer, e.g. 'user', 'foundry'
        remark (str): extra info about the layer

    Returns:
        Xsection: Xsection object having name <xsection>
    """

    _layer, _datatype = get_layer(layer)
    if growx is None:
        growx = cfg.default_xs['growx'][0]
    if growy is None:
        growy = cfg.default_xs['growy'][0]
    if origin is None:
        origin = cfg.default_xs['origin'][0]

    if cfg.layer_table.empty:
        add_layer(layer=layer, accuracy=accuracy)
    else:
        mask = (cfg.layer_table['layer'] == _layer) &\
            (cfg.layer_table['datatype'] == _datatype)
        hasLayerDatatype = cfg.layer_table[mask]
        if hasLayerDatatype.empty:
            add_layer(layer=layer, accuracy=accuracy)

    if xsection not in cfg.XSdict.keys():
        add_xsection(name=xsection)

    #check on layer, datatype instead of name for doubles:
    if not cfg.xsection_layer_map.empty:
        select = (cfg.xsection_layer_map['xsection']==xsection) &\
            (cfg.xsection_layer_map['layer']==_layer) &\
            (cfg.xsection_layer_map['datatype']==_datatype)
        LayerDatatype = cfg.xsection_layer_map[select]
        if not LayerDatatype.empty:
            cfg.xsection_layer_map.loc[select, ['growx', 'growy']] = growx, growy
            cfg.xs_list = cfg.xsection_layer_map['xsection'].unique()
            merge_xsection_layers_with_layers()
            return cfg.XSdict[xsection]
    else:
        cfg.default_xs_name = xsection # set default xs to first user defined xs.

    #new layer entry
    D = {
        'layer': int(_layer),
        'datatype': int(_datatype),
        'xsection': xsection,
        'xsection_foundry': fab_name,
        'growx': growx,
        'growy': growy,
        'accuracy': accuracy,
        'origin': origin,
        'remark': remark
    }
    cfg.xsection_layer_map = cfg.xsection_layer_map.append(D, ignore_index=True)

    merge_xsection_layers_with_layers()
    return cfg.XSdict[xsection]

add_xs = add_layer2xsection
add_xsection_layer = add_layer2xsection


def load_layers(filename):
    """Load layer definitions from csv file.

    Args:
        filename (str): name of layer definition file in csv format

    Returns:
        DataFrame: table with layers
    """

    #TODO: add to existing table, not replace
    #TODO: check if layers are unique when adding tables.
    cfg.layer_table = pd.read_csv(filename, delimiter=',')
    cfg.layer_table.dropna(how='all', inplace=True)
    cfg.layer_table.index.name = 'id'
    #delete unknown_layers when overwriting cfg.layer_table
    global unknown_layers
    unknown_layers = set()

    #TODO: add filename to column to track origin of content:
    #cfg.layer_table['filename'] = filename

    #create a dict for speeding up get_layer function
    tab = cfg.layer_table.set_index('layer_name')
    tab2 = tab.loc[:, ['layer', 'datatype']].stack().unstack(0)
    cfg.layerdict = tab2.to_dict(orient='list')

    merge_xsection_layers_with_layers()
    return cfg.layer_table


def load_xsection_layer_map(filename):
    """Load the assignment of layers to xsections from <filename>.

    Layers or xsections that do not yet exist are created automatically.

    Args:
        filename (str): xsection file in csv format

    Returns:
        DataFrame: merged xsections and layers.
    """

    cfg.xsection_layer_map = pd.read_csv(filename, delimiter=',')
    cfg.xsection_layer_map.dropna(inplace=True)

    #add xsections in the map to the xsection_table if missing:
    for name in cfg.xsection_layer_map['xsection']:
        #if name not in cfg.xsection_table['xsection'].values:
        add_XSdict(name)

    return merge_xsection_layers_with_layers()


def load_xsections(filename):
    """Load list of xsection from file with stub mapping.

    In addition create a stub map dictionary of a xsection to its stub-xsection.

    Args:
        filename (str): xsection map filename in csv format

    Returns:
        DataFrame: table with loaded xsections
    """

    cfg.xsection_table = pd.read_csv(filename, delimiter=',')
    for name in cfg.xsection_table['xsection']:
        add_XSdict(name)

    temptable = cfg.xsection_table.set_index('xsection')
    cfg.stubmap = temptable['stub'].dropna().to_dict()

    #TODO: check naming consistency of stubs
    #TODO: check for unique names
    #TODO: fill missing nazca_name with fab_name?
    #TODO: check if all xs are present in the xs definition.

    #update mask_layers:
    merge_xsection_layers_with_layers()
    return cfg.xsection_table


def add_xsection_stub(xsection, stub):
    """Assign a stub xsection to a xsection.

    The default stub is the xsection itself.
    """
    cfg.stubmap[xsection] = stub


def load_bbmap(filename):
    """Load the mapping table for nazca to spt.

    The spt format is used for third party OptoDesigner mask assembly.

    Args:
        filename (str): name of file containing bb-mapping in csv format

    Returns:
        None
    """
    bbmaptable = pd.read_csv(filename, delimiter=',')
    bbmaptable.dropna(how='all', inplace=True)
    columns = ['Nazca', 'OD_BBname', 'OD_port', 'odx', 'ody', 'oda', 'x', 'y', 'a']
    cfg.bbmap = bbmaptable[columns]
    cfg.spt = True
    return None


def merge_xsection_layers_with_layers():
    """Create a dictionary containing tables of xsections joined with layers.

    The xsection names are the dictionary keys.
    The dictionary contains all attributes for mask export.

    Left-join xsection table with layer table on (layer, datatype) and
    chop it up into a dictionay with xsection as key and the join result as
    values.

    Returns:
        None
    """
    if cfg.layer_table.empty or\
        cfg.xsection_layer_map.empty or\
        cfg.xsection_table.empty:
        return None

    #TODO: check if column names exists in cfg.xsection_layer_map.empty
    #TODO: deal with NaN entries, i.e. no (layer, datatype) match.

    right = cfg.layer_table[layer_table_attr_csv].rename(
        columns=dict(zip(layer_table_attr_csv, layer_table_attr_join)))
    left = cfg.xsection_layer_map[xsection_layer_map_attr_csv].rename(
        columns=dict(zip(xsection_layer_map_attr_csv, xsection_layer_map_attr_join)))
    Merged = pd.merge(left=left, right=right, how='left',
        on=['layer', 'datatype'])

    #reset the table when merge is called for a fresh build.
    cfg.xs_list = cfg.xsection_table['xsection'].unique()
    #cfg.mask_layers = dict()

    for xs in cfg.xs_list:
        rows = Merged['xsection'] == xs
        add_xsection(xs).mask_layers = Merged[rows][mask_layers_attr]

    return None


def load_masklayers(layer_file=None, xsection_layer_file=None):
    """Load layer and xsection files and merge them.

    This function combines
    1. load_layers()
    2. load_xsection_layer_map()
    3. merge()

    Args:
        layer_file (str): layer file name of csv file
        xs_file (str): xsection file name of csv file

    Returns:
        dict: {xsection name: DataFrame with mask layer}
    """
    load_layers(filename=layer_file)
    load_xsection_layer_map(filename=xsection_layer_file)
    return merge_xsection_layers_with_layers()


def get_xsection(name):
    """Return the Xsection object corresponding to <name>.

    Args:
        name (str | Node): xsection name or a pin object

    Returns:
        Xsection: Xsection object with name <name>
    """

    try: # if name is a pin (Node) get its xsection name
        name = name.xs
    except:
        pass
    if not name in cfg.XSdict.keys():
        msg = "\nNo xsection object existing under xsection name '{}'.".format(name)
        msg += ' Available xsections are:'
        for xname in cfg.XSdict.keys():
            msg += "\n  '{}'".format(xname)
        msg += "\nAlternatively, add a new xsection as:"
        msg += "\n  add_xsection(name='{}')".format(name)
        raise Exception(msg)

    return cfg.XSdict[name]

#getXS = get_xsection


def get_layer(layer):
    """Get the layer number for specific layer reference.

    If the <layername> does not exist a default layer number is returned.

    Args:
        layer (str| tuple | int): layer reference by name | (layer, datatype) | layer

    Returns:
        tuple: (layer, datatype)
    """

    #print(layername)
    #TODO: check if int or tuple are in layer definition, if not: warn once
    layID = None
    if isinstance(layer, int):
        layID = (layer, 0)
    elif isinstance(layer, tuple):
        if len(layer) == 2:
            layID = layer
        else:
            print("Error: invalid tuple format for layer {}. Length must be 2.".\
                format(layer))
    else:
        try:
            l, d = cfg.layerdict[layer]
            layID = (int(l), int(d))
        except:
            pass
    if layID is None: #layer not (yet) in unknown_layers:
        text = ''
        if layer in cfg.default_layers:
            layID = cfg.default_layers[layer]
            text = 'internal Nazca '
            warning = False #return layID # avoid warnings
        else:
            if isinstance(cfg.default_dump_layer, int):
                layID = (int(cfg.default_dump_layer), 0)
            elif isinstance(cfg.default_dump_layer, tuple):
                if len(cfg.default_dump_layer) == 2:
                    layID = cfg.default_dump_layer
            warning = True
        if warning and layer not in unknown_layers:
            if isinstance(layer, str):
                print("Warning: {0}layer '{1}' not set. "\
                      "Setting it now and redirecting output to layer {2}. "\
                      "You can explicitly add a new layer by: "\
                      "add_layer(name='{1}', layer=<num>)\n"
                      "Available layers are:\n{3}".\
                      format(text, layer, layID, sorted(list(cfg.layerdict.keys()))))
            else:
                print("Warning: layer {0} not set. "\
                      "Redirecting output to collecting layer {1}. "\
                      "You can explicitly add a new layer by: "\
                      "add_layer(name=<name>, layer={0})\n"\
                      "Available layers are:\n{2}".\
                     format(layer, layID, sorted(list(cfg.layerdict.keys()))))
            unknown_layers.add(layer)
            add_layer(name=str(layer), layer=cfg.default_dump_layer)

    if layID not in cfg.layerset:
        cfg.layerset.add(layID)
        add_layer(name=str(layID), layer=layID, unknown=True)
    return layID


def clear_layers():
    """Drop all rows from layers and xsection tables"""
    cfg.layer_table.drop(cfg.layer_table.index, inplace=True)
    cfg.xsection_layer_map.drop(cfg.xsection_layer_map.index, inplace=True)
    cfg.colors.drop(cfg.colors.index, inplace=True)
    for xs in cfg.xs_list:
        get_xsection(xs).mask_layers = dict()
    cfg.xs_list = []


def clear_xsections():
    """Drop all rows from layers and xs tables"""
    cfg.xsection_layer_map.drop(cfg.xsection_layer_map.index, inplace=True)
    #cfg.mask_layers = dict()
    for xs in cfg.xs_list:
        get_xsection(xs).mask_layers = dict()
    cfg.xs_list = []


def empty_xsections():
    """Delete info on all xs. Keep layer info."""
    cfg.xsection_layer_map = pd.DataFrame()
    merge_xsection_layers_with_layers()


def load_parameters(filename):
    cfg.dfp = pd.read_csv(filename, delimiter=',')


def get_parameter(name):
    """get parameter value for specific name"""
    return float(cfg.dfp[cfg.dfp['name'] == name]['value'])


#==============================================================================
# plt
#==============================================================================
def set_plt_properties(figsize=14, cmap=None, N=32, alpha=0.3):
    """
    Set the default colormap to use in Matplotlib output for mask layout.

    Args:
        figsize (float): size of the matplotlib figure
        cmap (str): name of the colormap
        N (int): the number of colors in the map in case of a 'linearSegmentedColormap'
        alpha (float): transparency of the cmap to see through layers.

    Returns:
        None
    """

    if cmap is None:
        cmap = cfg.plt_cmap_name

    mp = mpl.cm.get_cmap(cmap)
    if isinstance(mp, mpl.colors.LinearSegmentedColormap):
        cfg.plt_cmap = []
        for n in range(N):
            cfg.plt_cmap.append(mp(n/N))
    else:
        cfg.plt_cmap = mp.colors

    cfg.plt_cmap_size = len(cfg.plt_cmap)
    cfg.plt_figsize = figsize
    cfg.plt_alpha = alpha

#make sure colors are set:
set_plt_properties()


#==============================================================================
# colors
#==============================================================================
color_defaults = {
    'frame_brightness': 0,
    'fill_brightness': 0,
    'fill_color': 0,
    'frame_color': 0,
    'dither_pattern': 'I0',
    'valid': True,
    'visible': True,
    'transparent': True,
    'width': 0,
    'marked': False,
    'animation': 1,
    'alpha' : 0.3
}


def load_layercolors(filename):
    """Read colormap from a .csv Nazca color table file.

    Returns:
        None
    """

    #For Klayout tabs with groups the layer column contains a '*'
    #and Pandas will create type->str
    #For tabs without groups the layer is read as type->int
    df1 = pd.read_csv(filename)
    df1['layer'] = df1['layer'].astype(str)
    #df1.loc[df1['width'] == 'na', 'width'] = 0
    df1.replace('na', '', inplace=True)

    #TODO: creates a brand new table and forgets present settings
    #  needs an update to handle mutliple layer/color sets.
    cfg.colors = df1[df1['layer'] != '*'].dropna(subset=['layer'])
    cfg.colors.reset_index(inplace=True)
    try:
        cfg.colors['layer'] = cfg.colors['layer'].astype(str)
        cfg.colors['datatype'] = cfg.colors['datatype'].astype(str)
    except:
        print('Warning: Check load_layercolors.')
        #print('Group in load_layercolors:', cfg.colors['layer'], '\n')
    #TODO: the try will trigger on any layer entry with text.
    # as a resuls all the layers and datatype remain a string.

    #cfg.colors.set_index(['layer', 'datatype'], inplace=True)
    cfg.colors['alpha'] = 0.3
    return None


# note: fill_color and frame_color will come from default colormap.
def set_layercolor(
    layer=None,
    frame_color=None,
    fill_color=None,
    frame_brightness=None,
    fill_brightness=None,
    dither_pattern=None,
    valid=None,
    visible=None,
    transparent=None,
    width=None,
    marked=None,
    animation=None,
    alpha=None):
    """Set layer color information.

    For missing values for fill_color and frame_color, the default_colors
    as specified in colormap cfg.plt_cmap will be applied. This colormap can
    be adjusted by the user.

    Args:
        layer (int | tuple): layer number or (layer, datatype)
        ...

    Returns:
        None
    """

    if layer is None:
        raise ValueError('Need to provide a layer number, or name.')
    else:
        layer, datatype = get_layer(layer)

    colors = {
        'layer': str(layer),
        'datatype': str(datatype),
        'fill_color': fill_color,
        'frame_color': frame_color,
        'frame_brightness': frame_brightness,
        'fill_brightness': fill_brightness,
        'dither_pattern': dither_pattern,
        'valid': valid,
        'visible': visible,
        'transparent': transparent,
        'width': width,
        'marked': marked,
        'animation': animation,
        'alpha' : alpha
    }

    if not cfg.colors.empty:
        mask = cfg.colors['layer'] == str(layer)
        mask2 = cfg.colors['datatype'] == str(datatype)
        select = cfg.colors.loc[mask & mask2]
        items = len(select)
    else:
        items = 0

    if items == 0:
        #TODO: not here or the colors can't be changed later on.
        #col = mpl.colors.to_hex(cfg.plt_cmap[layer % cfg.plt_cmap_size])
        col = rgb2hex(cfg.plt_cmap[int(layer % cfg.plt_cmap_size)])

        if fill_color is None:
            colors['fill_color'] = col
        if frame_color is None:
            colors['frame_color'] = col
        for attr in colors.keys():
            if colors[attr] is None:
                if attr in color_defaults.keys():
                    colors[attr] = color_defaults[attr]

        df = pd.DataFrame(colors, index=[0])
        cfg.colors = pd.concat([cfg.colors, df], ignore_index=True)
        cfg.colors['layer'] = cfg.colors['layer'].astype(str)
        cfg.colors['datatype'] = cfg.colors['datatype'].astype(str)
    elif items == 1:
        for attr in colors.keys(): #use existing values for undefined attributes
            if colors[attr] is None:
                #TODO: Better use a try to catch color_defaults are missing:
                if attr in color_defaults.keys():
                    colors[attr] = select[attr].iloc[0]
        df = pd.DataFrame(colors, index=select.index) #.set_index(idx)

        # TODO: this could be a new (extra) layername with an old layer number
        # now it's always overwritten when add_layer() is called:
        cfg.colors.loc[mask & mask2] = df
    else:
        print('Multiple color entries for ML({}, {}). Skipping...'.format(layer, datatype))


#==============================================================================
# Show or get tables
#==============================================================================
def get_layers():
    df = cfg.layer_table[['nazca_name', 'layer', 'datatype', 'accuracy']]
    return df

def show_xsections():
    """Print the xsections table.

    Returns:
        None
    """
    #print('-------------------------------')
    print('xsections:')
    if not cfg.xsection_table.empty:
        pprint(cfg.xsection_table[xsection_table_attr])
    else:
        print('No xsections defined.')

def show_layers():
    """Print the layer table."""
    #print('-------------------------------')
    #print("DataFrame 'cfg.layer_table' with layer definitions:")
    print("layers:")
    try:
        df = cfg.layer_table[['layer_name', 'layer', 'datatype', 'accuracy']]
        pprint(df)
    except:
        print('  ...no layers defined.')

def show_layercolors():
    """Print the layercolor table."""
    #print('-------------------------------')
    #print("DataFrame 'cfg.layer_table' with layer definitions:")
    print("layercolors:")
    df = cfg.colors[['fill_color', 'frame_color', 'width']]
    pprint(df)


def show_xsection_layer_map():
    """Print the xsection_layer_map table to stdout.

    Returns:
        None
    """
    #print('-------------------------------')
    print("xsection-layer_map:")
    if not cfg.xsection_layer_map.empty:
        df = cfg.xsection_layer_map[['xsection', 'layer', 'datatype', 'growx', 'growy']]
        pprint(df)
        #return df
    else:
        print('No xs defined.')


def show_mask_layers():
    """Print mask_layers dictionary with layer export definitions to stdout.

    Returns:
        None
    """
    #print('-------------------------------')
    #print("Dictionary of DataFrames 'cfg.xs_layers' with the description of each xs:")
    columns = ['layer', 'datatype', 'growx', 'growy', 'accuracy']
    print("mask_layers:")
    for xs in cfg.xs_list:
        #for sorted(cfg.mask_layers.keys()):
        under = '_'*(29-len(xs))
        print("xsection: '{}' {}".format(xs, under))
        pprint(get_xsection(xs).mask_layers.reset_index()[columns])
        print('')




