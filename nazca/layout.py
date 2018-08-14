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
# (c) 2016-2017 Ronald Broeke
#

"""Generate layout to gds or screen."""

from collections import defaultdict
from math import sin, cos, radians
import pandas as pd
import os.path
from IPython.core.display import display, HTML
from pprint import pprint

import matplotlib.pyplot as mplt
from matplotlib.patches import Polygon as matPolygon
from matplotlib.collections import PatchCollection

from . import gds as gdsmod
from . import gds_base as gbase
from . import gds_import as gstream
from .netlist import Netlist, Cell, Pointer, Polygon, Polyline, Annotation
from .mask_layers import add_layer, get_layer
#from .pdk_template_core import strip_groupname
from . import cfg
from . import util

try:
    import svgwrite
    cfg.SVGWRITE = True
except:
    cfg.SVGWRITE = False


# =============================================================================
# Nazca cell
# =============================================================================
class ClsNazca():
    """Helper class to handle Nazca cell tree transformations.

    Note that this class does not reconstruct the netlist in the
    cell it returns, but manipulates it on element and layer level.

    Use cases:
    flatten a celltree, change and/or delete or filter out layers,
    change cell names.
    """
    def __init__(self, instantiate=True, cellmap=None, layermap=None,
            flat=None, infolevel=0):
        self.cellsopen = []
        self.flat = flat
        self.instantiate = instantiate
        self.layermap = layermap
        self.cellmap = cellmap
        self.oldcellrefs = {} # link old cell to new cell
        self.newcellrefs = {} # link new cell to old cell
        self.topcell = None
        self.infolevel = infolevel


    def open(self, params):
        """Open new Cell object based on namedtuple <cellinfo>.

        Args:
            cell (Cell): original cell to be copied/filtered
            level (int): cell level in hierarchy (0 is topcell)

        """
        self.level = params.level
        if params.cell_name is None:
            cell_name = params.cell.cell_name + '_N'
        else:
            cell_name = params.cell_name
        if self.infolevel > 0:
            print("{}.{}ClsNazca: open  '{}' source:'{}'".\
                format(self.level, '  '*self.level, cell_name, params.cell.cell_name))
        newcell = Cell(name=cell_name, instantiate=self.instantiate)
        self.cellsopen.append(newcell)
        self.oldcellrefs[newcell] = params.cell
        self.newcellrefs[params.cell] = newcell
        if params.level == 0:
            self.topcell = newcell
        return newcell

    def close(self):
        """Close Cell object."""
        newcell = self.cellsopen.pop()
        newcell.close()
        if self.infolevel > 0:
            print("{}.{}ClsNazca: close '{}'".\
                format(self.level, '  '*self.level, newcell.cell_name))
        return newcell

    def add_polygon(self, layer, xy):
        """"""
        Polygon(points=xy, layer=layer).put(0)
        #print("ClsNazca: add_polygon to", cfg.cells[-1].cell_name, xy[:3], layer)


    def add_polyline(self, layer, xy):
        """"""
        Polyline(points=xy, layer=layer).put(0)
        #print("ClsNazca: add_polylines")


    def add_annotation(self, text, layer, pos):
        """"""
        #print("ClsNazca: add_annotation")
        Annotation(text=text, layer=layer).put(*pos)


    def add_instance(self, inode, xya, flip):
        """"""
        if self.infolevel > 0:
            print("{}.{}ClsNazca: add_instance '{}' xya:{}, flip:{}".\
                format(self.level, '  '*self.level, inode.cell.cell_name, xya, flip))
        cell = inode.cell
        if cell.instantiate: #and not flat:
            #TODO, flatten array's here if not instantiated?
            self.newcellrefs[cell].put(*xya, flip=inode.flip, array=inode.array)


    def add_gds(self):
        """"""
        if self.infolevel > 0:
            print("ClsNazca: add_gds")


# =============================================================================
# GDS
# =============================================================================
manual_nazca_version = None #override auto nazca versioning if not Nonne.
class ClsGDS():
    """Helper class to handle gdsii compatible export of masks."""

    def __init__(self):
        pass

    def open(self, filebasename):
        self.filename = filebasename+'.gds'
        self.outfile = open(self.filename, 'wb')
        self.outfile.write(gdsmod.layout_open(name=manual_nazca_version))
        self.content = dict() # collect cell content per level.

    def write(self, level):
        self.outfile.write(b''.join(self.content[level]))

GDS = ClsGDS()

#==============================================================================
# svg
#==============================================================================
class ClsSVG():
    """Helper class to handle svg compatible export of masks."""

    def __init__(self):
        #if filename is not None:
        #    filename = os.path.splitext(filename)[0] + '.svg'
        self.minx = 1e6
        self.miny = 1e6
        self.maxx = -1e6
        self.maxy = -1e6

    def open(self, filebasename='nazca_export'):
        """Open SVG.

        Open after layer colors have been loaded.

        Returns:
            None
        """
        self.filename = filebasename+'.svg'
        # create a map of defined layer colors to speed up lookup
        self.drawing = svgwrite.Drawing(filename=self.filename, size=('100%', '100%'),
            profile='tiny', debug=False)
        self.layermap = defaultdict(list)
        if not cfg.colors.empty:
            for i, row in cfg.colors.loc[:, ['layer', 'datatype']].iterrows():
                #int('1.0') doesn't work, int(float('1.0')) does
                ML = (int(float(row[0])), int(float(row[1])))
                self.layermap[ML].append(i)
        self.cmap = cfg.plt_cmap
        self.cmaplen = len(self.cmap)
        self.alpha = cfg.plt_alpha


    def add_polygon(self, layer, points, bbox):
        """Add a polygon to the SVG drawing."""

        lay = (int(layer[0]), int(layer[1]))
        for colorentry in self.layermap[lay]:
            colors = cfg.colors.loc[colorentry]
            if isinstance(colors, pd.DataFrame):
                colors = colors.iloc[0] #TODO: make loop , not only first color
            if not colors['visible']:
                continue
            edgecolor = colors['frame_color']
            facecolor = colors['fill_color']
            alpha = colors['alpha']
            lw = colors['width']
            if colors['dither_pattern'] == 'I1':
                fill_opacity = 0
            else:
                fill_opacity = alpha
            self.g = self.drawing.g(
                stroke=edgecolor,
                stroke_opacity=alpha,
                fill=facecolor,
                fill_opacity=fill_opacity,
                stroke_width=lw)
            self.mask = self.drawing.add(self.g)

            p = self.drawing.polygon(points=points, transform="scale(1,-1)")
            self.mask.add(p)

            if self.minx > bbox[0]:
                self.minx = bbox[0]
            if self.maxx < bbox[2]:
                self.maxx = bbox[2]
            if self.miny > bbox[1]:
                self.miny = bbox[1]
            if self.maxy < bbox[3]:
                self.maxy = bbox[3]

    def close(self):
        x, y, w, h = self.minx, self.miny, self.maxx-self.minx, self.maxy-self.miny
        self.drawing.viewbox(minx=x, miny=-y-h, width=w, height=h)
        self.drawing.save()

SVG = ClsSVG()


#==============================================================================
# matplotlib export
#==============================================================================
class ClsMatplotlib():
    """Helper class to handle matplotlib export of masks."""

    def __init__(self):
        """Construct an object handling Matplotlib exports."""

      #  try:
      #      font = cfg.matplotlib_font
      #
      #  except:
      #  font = {
      #      #'family'  : 'normal',
      #      'style'  : 'normal',
      #      'weight' : 'light', #'light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'
      #      'size'   : cfg.plt_fontsize}
      #  mplt.rc('font', **font)


    def open(self, cell=None, title=None):
        """Inititialize Matplotlib mask output."""
        font = {
            #'family'  : 'normal',
            'style'  : 'normal',
            'weight' : 'light', #'light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'
            'size'   : cfg.plt_fontsize}
        mplt.rc('font', **font)



        self.cell = cell

        #print('TITLE', title)
        if title is not None:
            #self.title = '{} - {}'.format(title, cell.cell_name)
            self.title = '{}'.format(title)
        else:
            self.title = cell.cell_name
            self.title = 'cell: {}'.format(self.title)

        self.cellname = cell.cell_name
        self.patches = []
        self.colors = []
        self.minx = 1e6
        self.miny = 1e6
        self.maxx = -1e6
        self.maxy = -1e6
        self.figsize = cfg.plt_figsize

        # create a map of type list of defined layer colors to speed up lookup
        self.layermap = defaultdict(list)
        if not cfg.colors.empty:
            for i, row in cfg.colors.loc[:, ['layer', 'datatype']].iterrows():
                ML = (int(float(row[0])), int(float(row[1])))
                self.layermap[ML].append(i)

        #print('self.layermap\n')
        #pprint(self.layermap)
        # backup colormap for layers without color setting
        self.cmap = cfg.plt_cmap
        self.cmaplen = len(self.cmap)
        self.alpha = cfg.plt_alpha


    def add_polygon(self, layer, points, bbox):
        """Add a polygon to the plot.

        Note that using linewidth > 0 significantly slows down drawing in Matplotlib.

        Args:
            layer (int, int): mask layer
            points (list of (float, float)): polygon points
            bbox (tuple of floats): bounding box of polygon (x1, y1, x2, y2)

        Returns:
            None
        """
        if layer == cfg.default_layers['docu_pin']:
            return None
        polygon = None
        lay = (int(layer[0]), int(layer[1]))
        if lay not in self.layermap:
            # use default colors if no explicit layer color has been set
            col = self.cmap[lay[0] % self.cmaplen]
            edgecolor = col
            facecolor = col
            alpha = self.alpha
            lw = 0
            polygon = matPolygon(
                points,
                closed=True,
                edgecolor=edgecolor,
                facecolor=facecolor,
                fill=True,
                lw=lw,
                alpha=alpha)

        else:
            # use preset layer colors
            # like in Klayout there can be multiple colors per layer.
            for colorentry in self.layermap[lay]:
                colors = cfg.colors.loc[colorentry]
                if isinstance(colors, pd.DataFrame):
                    colors = colors.iloc[0] #TODO: make loop , not only first color
                if not colors['visible']:
                    continue
                edgecolor = colors['frame_color']
                facecolor = colors['fill_color']
                alpha = colors['alpha']
                lw = colors['width']
                try: #needed because lw == '' occurs
                    lw = float(lw)
                except:
                    lw = 0
                if colors['dither_pattern'] == 'I1':
                    fill = False
                else:
                    fill = True
                #lw = 0 #matplotlib very slow with lw != 0

                polygon = matPolygon(
                    points,
                    closed=True,
                    edgecolor=edgecolor,
                    facecolor=facecolor,
                    fill=fill,
                    lw=lw,
                    alpha=alpha)

        if polygon is not None:
            self.patches.append(polygon)
            if self.minx > bbox[0]:
                self.minx = bbox[0]
            if self.maxx < bbox[2]:
                self.maxx = bbox[2]
            if self.miny > bbox[1]:
                self.miny = bbox[1]
            if self.maxy < bbox[3]:
                self.maxy = bbox[3]

            addbbox = False
            if addbbox:
                polygon = matPolygon(
                        [(bbox[0], bbox[1]), (bbox[2], bbox[1]), (bbox[2], bbox[3]), (bbox[0], bbox[3])],
                        closed=True,
                        edgecolor='#00FF00',
                        facecolor=facecolor,
                        fill=False,
                        lw=1,
                        alpha=alpha)
                self.patches.append(polygon)


    def add_polyline(self, layer, points, width, bbox):
        """Add a polyline to the plot."""


    def close(self):
        """Close plt and show plot."""

        if not self.patches:
            print('Nothing to draw.')
            return None

        dx = self.maxx-self.minx
        dy = self.maxy-self.miny
        bufx = 0.05*dx
        bufy = 0.05*dy
        buf = max(bufx, bufy)
        domain = max(dx, dy)
        xtot = dx+2*buf
        ytot = dy+2*buf
        aspect = ytot/xtot
        factor = 1.0
        if aspect > factor:
            ysize = self.figsize
            xsize = ysize/aspect
        else:
            xsize = self.figsize
            ysize = 1.2*xsize*aspect

        fig, ax = mplt.subplots(figsize=(xsize, ysize),\
            facecolor=cfg.plt_background_outside)
        p = PatchCollection(self.patches, match_original=True)

        ax.add_collection(p)
        ax.set(xlim=[self.minx-buf, self.minx+dx+buf],
             ylim=[self.miny-buf, self.miny+dy+buf])

        # format spines:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.spines['left'].set_visible(True)
        ax.spines['left'].set_smart_bounds(True)

        ax.spines['bottom'].set_visible(True)
        ax.spines['bottom'].set_smart_bounds(True)

        ax.yaxis.set_ticks_position('left')
        #ax.yaxis.set_ticklabels('')
        ax.xaxis.set_ticks_position('bottom')

        ax.set_title(self.title)

        # cell info available under self.cell:
        try:
            pin_ignore = cfg.pin_settings['bb_docu_pin_ignore']
        except:
            pin_ignore = []
        if cfg.pin_settings['bb_pin_layer'] != cfg.default_layers['docu_pin']:
            dt = 0.02*domain
            for name, pin in self.cell.pin.items():
                if name == 'org' or\
                    (pin.type == 'bbox' and  not cfg.bbox_stubs) or\
                    name in pin_ignore:
                       continue
                at = pin.pointer.a
                xt = pin.pointer.x + dt*cos(radians(at))
                yt = pin.pointer.y + dt*sin(radians(at))
                ax.text(xt, yt, name, ha='center', va='center', rotation=0, wrap=True)
        else: #docu-pins
            try:
                pin_scale = cfg.pin_settings['bb_docu_pin_scale']
            except:
                pin_scale = 1/15
                print("No 'bb_docu_pin_scale' set. Using {}".format(pin_scale))
            try:
                pin_shape = cfg.pin_settings['bb_docu_pin_shape']
            except:
                pin_shape = 'arrow_full'
                print("No 'bb_docu_pin_shape' set. Using '{}'".format(pin_shape))

            docu_patches = []
            d_text = 0.40
            fontsize = 20 * 15*pin_scale
            scale = domain*pin_scale
            shape = cfg.pinshapes[pin_shape]
            dt = -d_text*pin_scale*domain
            for name, pin in self.cell.pin.items():
                #not pin.show or
                if name == 'org' or\
                    (pin.type == 'bbox' and  not cfg.bbox_stubs) or\
                    name in pin_ignore:
                        continue
                at = pin.pointer.a
                xt = pin.pointer.x + dt*cos(radians(at))
                yt = pin.pointer.y + dt*sin(radians(at))

                ax.text(xt, yt, name, fontsize=fontsize, rotation=0.,
                    ha="center", va="center", weight='bold'
                    #bbox=dict(boxstyle="round", ec=(1, 1, 1), fc=(1, 1, 1))
                )
                points = util.transform_polygon(points=shape,
                    da=at, scale=scale, x=pin.pointer.x, y=pin.pointer.y)
                polygon = matPolygon(
                    xy=points,
                    closed=True,
                    edgecolor='#000000',
                    facecolor='#FFFFFF',
                    fill=True,
                    lw=3,
                    alpha=0.5)
                docu_patches.append(polygon)
            if docu_patches:
                p = PatchCollection(docu_patches, match_original=True)
                ax.add_collection(p)

        mplt.axis('scaled')
        mplt.tight_layout()

        #ax.set_facecolor(cfg.plt_background_inside)
#==============================================================================
#         #Watermark:
#         ax.set_facecolor((0.97, 0.97, 0.97))
#         import matplotlib.image as image
#         datafile = '/nazca/nazca-logo_illustration_beeldmerk.png'
#         im = image.imread(datafile)
#         mplt.imshow(im, extent=[self.minx-bufx, self.minx+dx+bufx, self.miny-bufy, self.miny+dy+bufy],
#             alpha=0.05, zorder=3)
#
#==============================================================================
PLT = ClsMatplotlib()


#==============================================================================
# spt
#==============================================================================
class ClsSPT():
    """Helper class to handle spt compatible export of masks."""

    def __init__(self):
        self.Fspt = None
        self.bblock_cnt = 0

    def open(self, filename=None):
        """Initialize spt file export."""
        self.Fspt = None # spt-file handle
        self.spt_filename = os.path.splitext(self.filename)[0] + '.spt'
        self.Fspt = open(self.spt_filename, 'w')
        self.Fspt.write(
"""
// -----------------------------------------------------------------
// Building block references for mask assembly with Phoenix tooling.
// This file may be required by the foundry.
// * Exported by Nazca Design *
// -----------------------------------------------------------------

""")

    def add_line(self, cell, x, y, a, flip):
        """Add building block line to spt file."""
        self.bblock_cnt += 1

        try:
            BBname = cell.BBmap[0]
            port = cell.BBmap[1]
            tx, ty, ta = cell.BBmap[3]
            fx, fy, fa = cell.BBmap[4] #parametrized translation
            params = cell.BBmap[5]
        except:
            message = "ERROR: no mapping for cell {}.".format(cell.cell_name)
            print(message)
            self.Fspt.write(message)
            return None

        p = Pointer(x, y, a)
        if flip:
            x, y, a = p.move(tx+fx, -ty-fy, -ta-fa).xya()
            port = port+'->box@flipPortY(org)'
        else:
            x, y, a = p.move(tx+fx, ty+fy, ta+fa).xya()
            port = port+'->box@org'

        varlist = []
        od_exceptions = ['box', 'fname_final'] # no quotes for these in spt
        for param in params:
            if isinstance(param, str) and param not in od_exceptions:
                c = '\"{}\"'.format(param)
                varlist.append(c)
            else:
                varlist.append(str(param))

        self.Fspt.write("ml::{0}({1}+[{2:.5f}, {3:.5f}, {4:.5f}] : {5}) CMP{6};\n".\
            format(BBname, port, x, y, a,', '.join(varlist), self.bblock_cnt))

    def close(self):
        """Close spt output routine."""
        if self.Fspt is not None:
            self.Fspt.close()
            print("...Wrote file '{}'".format(self.spt_filename))

SPT = ClsSPT()


# =============================================================================
# Nazca building block
# =============================================================================
class BBlock():
    """Helper class to export a layout as a Nazca GDS building block.
    """

    def __init__(self, filebasename):
        """Initialize BB file: open and add header.

        Args:
            filename (str): basename of the output file (without extension).

        return:
            None
        """
        self.filebasename = filebasename
        pyfilename = os.path.join(self.filebasename)+'.gds.py'
        self.bblockfile = open(pyfilename, 'w')
        bbfile_header = \
"""# Nazca building block(s) file.
#
# Note1: Do *not* edit coordinates in the cell definition(s) below.
# Note2: Always use this file with the intended, accompanying gds library.
#
# In order to use the BB(s) in this file,
# include the following lines in your design file:
#
# import {0}
# {0}.<bb_name>.put()
#
# where '{0}' is the name of this file without .py extension, and
# where <bb_name> needs to be replaced with any cell variable name in this file.
#
# Executing this file stand-alone exports a gds with all BBs in this file.

from nazca import Cell, Pin, load_gds, export_gds
""".format(pyfilename[:-3], pyfilename)

        self.bblockfile.write(bbfile_header)
        self.xpos = 0
        self.foot = ''


    def body(self, topcell):
        """Write BB file body element.

        Returns:
            None
        """
        indent = '    '
        gdsfilename = os.path.join(self.filebasename+'.gds')
        cellname = topcell.cell_paramsname
        cellvar = cellname.replace('.', '_')
        self.bblockfile.write(
            "\nwith Cell(name='BB', autobbox=False, instantiate=False) as {0}:\n".\
            format(cellvar))
        self.bblockfile.write("{}load_gds(filename='{}', cellname='{}').put(0)\n".\
            format(indent, gdsfilename, cellname))
        self.foot += ("{}{}.put('org', {:.1f})\n".\
            format(indent, cellvar, self.xpos))

        try:
            length = topcell.length
        except:
            length = 0
        if length is not None:
            self.xpos += 1.1*length
        for name, pin in sorted(topcell.pin.items()):
            if name == 'org':
                continue
            x, y, a = pin.pointer.xya()
            if pin.xs is None:
                xs = None
            else:
                xs = "'{}'".format(pin.xs)
            self.bblockfile.write(
                "{}Pin(name='{}', width={}, xs={}).put({:.5f}, {:.5f}, {:.5f})\n".\
                format(indent, name, pin.width, xs, x, y, a))
        return None


    def footer(self):
        """Write BB file footer.

        Returns:
            None
        """
        self.bblockfile.write("\nif __name__ == '__main__':\n")
        self.bblockfile.write(self.foot)
        self.bblockfile.write("    export_gds()\n")
        return None



tab = '  '
class Export_layout():
    """Class to export a mask layout.
    """

    def __init__(self, layermap=None, layermapmode='all'):
        self.gds_files = dict() # existing gds files used in the layout
        self.layermapmodes = ['none', 'all']
        self.reset()

        #make sure hull layer is known
        if cfg.export_hull:
            if 'hull' in cfg.default_layers.keys():
                get_layer(cfg.default_layers['hull'])

        self.setlayermap(layermap=layermap, mode=layermapmode)
        return None


    def reset(self):
        """Reset export settings.

        Returns:
            None
        """
        self.topcells = None
        self.path = ''
        self._filename='nazca_export'
        self.bblock = False
        self.infolevel = 0
        self.show_cells = False
        self.clear = cfg.export_clear
        self.title = None
        self.output = None
        self._ascii = False #gds as ascii
        self._gds = False
        self._flat = False
        self._spt = False
        self._plt = False
        self._svg = False
        self.nazca = False
        self.layermap = {}
        self.layermapmode = 'all'
        self.setlayermap()
        return None


    def newtree(self):
        """Create a Nazca export object."""
        self._newtree = ClsNazca()
        return self._newtree


    def setlayermap(self, layermap=None, mode='all'):
        """Create the layermap for export.

        Layermap is set internally.

        Returns:
            dict: {layer: export_layer}, layermap as reference for the user
        """
        if layermap is None:
            layermap = {}

        if mode is None:
            mode = 'all'
        elif mode.lower() in self.layermapmodes:
            self.layermapmode = mode.lower()
        else:
            print("Warning: Trying to set a non-existing layermapmode '{}'."\
                " Valid options as {}".format(mode, self.layermapmodes))
            self.layermapmode = 'all'

        if self.layermapmode == 'all':
            for layer in cfg.layerset:
                if layer is not None:
                    self.layermap[layer] = layer
        elif self.layermapmode == 'none':
            for layer in cfg.layerset:
                if layer is not None:
                    self.layermap[layer] = None

        if layermap is not None:
            for layerold, layernew in layermap.items():
                layerold = get_layer(layerold)
                if layerold in cfg.layerset:
                    layer_new = get_layer(layernew)
                    if layer_new not in cfg.layerset:
                        print("Warning (rebuild): mapping to non-exiting layer {0}."\
                              " Adding it for you, but honestly, you should add it:"
                              " add_layer(layer={0})".format(layernew))
                        add_layer(name=str(layernew), layer=layernew)
                    self.layermap[layerold] = layernew

        return self.layermap

    @property
    def filebasename(self):
        return self._filebasename

    @property
    def filename(self):
        return self._filename
    @filename.setter
    def filename(self, name):
        if name is None:
            name = self._filename
        if isinstance(name, str):
            last = name[-4:].lower()
            base = name[:-4]
            if last == '.gds':
                self.gds = True
                self._filebasename = base
            elif last == '.svg':
                self.svg = True
                self._filebasename = base
            elif last == '.plt':
                self.plt = True
                self._filebasename = base
            else:
                self._filebasename = name
        else:
            print("WARNING: Filename provided is not a string but a {}."\
                 " Using '{}' instead.".\
                     format(type(name), self._filebasename))
        return self._filebasename


    @property
    def spt(self):
        return self._spt
    @spt.setter
    def spt(self, value):
        if value:
            self.flat = True
        self._spt = value
        return self._spt

    @property
    def svg(self):
        return self._svg
    @svg.setter
    def svg(self, value):
        if value:
            self.flat = True
        self._svg = value
        return self._svg

    @property
    def plt(self):
        return self._plt
    @plt.setter
    def plt(self, value):
        if value:
            self.flat = True
        self._plt = value
        return self._plt

    @property
    def flat(self):
        return self._flat
    @flat.setter
    def flat(self, value):
        if not value:
            self._spt = False
            self._svg = False
            self._plt = False
        self._flat = value
        return self._flat

    @property
    def ascii(self):
        return self._ascii
    @ascii.setter
    def ascii(self, value):
        if value:
            self._gds = True
        self._ascii = value
        return self._ascii

    @property
    def gds(self):
        return self._gds
    @gds.setter
    def gds(self, value):
        if not value:
            self._ascii = False
        self._gds = value
        return self._gds


    def add_polygons(self, pgon_iter=None, level=None, params=None):
        """Add polygons content to cell.

        Args:
            params (namedtuple): replaces all other parameters if not None

        Returns:
            None
        """
        if params is not None:
            pgon_iter = params.iters['polygon']
            level = params.level

        for pgon, xy, bbox in pgon_iter:
            layer = self.layermap[pgon.layer]
            if layer is None:
                continue
            if self.gds:
                GDS.content[level].append(
                    gbase.gds_polygon(xy, lay=int(layer[0]),
                        datatype=int(layer[1])))
            if self.plt:
                PLT.add_polygon(layer, xy, bbox)
            if self.svg:
                SVG.add_polygon(layer, xy, bbox)
            if self.nazca:
                self.CELL.add_polygon(layer, xy)
        return None


    def add_polylines(self, pline_iter=None, level=None, params=None):
        """Add polylines content to cell.

        In Matplotlib and svg output the polylines are converted to polygons.

        Args:
            params (namedtuple): replaces all other parameters if not None

        Returns:
            None
        """
        if params is not None:
            pline_iter = params.iters['polyline']
            level = params.level

        for pline, xy, bbox in pline_iter:
            if self.gds:
                GDS.content[level].append(
                    gbase.gds_polyline(xy, pline.width, lay=int(pline.layer[0]),
                        datatype=int(pline.layer[1]), pathtype=pline.pathtype))
            if self._plt or self._svg:
                polygon_xy = util.polyline2polygon(xy, width=pline.width, miter=0.5)
                if self.plt:
                    PLT.add_polygon(pline.layer, polygon_xy, bbox)
                if self.svg:
                    SVG.add_polygon(pline.layer, polygon_xy, bbox)
                if self.nazca:
                    self.CELL.add_polylines(pline.layer, polygon_xy, bbox)
        return None


    def add_annotations(self, anno_iter=None, level=None, create=False, cell=None, params=None):
        """Add annotation content to cell.

        Args:
            instantiate (bool): the instantiation level of the cell where
               the annotation are from: To check for black boxes, which
               can't be flattened.
            params (namedtuple): replaces all other parameters if not None

        Returns:
            None
        """
        if params is not None:
            anno_iter = params.iters['annotation']
            level = params.level
            create = params.cell_create
            cell = params.cell

        for anno, xy in anno_iter:
            layer = int(anno.layer[0])
            datatype = int(anno.layer[1])
            if (layer, datatype) == cfg.default_layers['bb_name'] :
                if not create:
                    pass
                    #print("ERROR: it's not allowed to flatten a black-box cell '{}'".format(cell.cell_name))
                name = anno.text.split('\n')
                if not cell.cell_name.startswith(name[0]):
                    pass
                    #print("ERROR: bb_name and cell_name of a black-box must be the same: '{}' != '{}'".\
                    #    format(name, cell.cell_name))
            if self.gds:
                GDS.content[level].append(
                    gbase.gds_annotation(xy=xy, string=anno.text,
                        lay=layer, datatype=datatype))
            if self.nazca:
                self.CELL.add_annotation(anno.text, anno.layer, xy)
        return None


    def add_instances(self, instance_iter=None, level=None, infolevel=0, params=None):
        """Add instances to cell.

        Args:
            params (namedtuple): replaces all other parameters if not None

        Returns:
            None
        """
        if params is not None:
            instance_iter = params.iters['instance']
            level = params.level

        for inode, [x, y, a], flip in instance_iter:
            cell = inode.cell
            if self.spt:
                try:
                    cell.BBmap[3]
                except:
                    continue
                SPT.add_line(cell, x, y, a, flip)
                continue

            if cell.instantiate and not self.flat:
                #only add references for instantiated cells
                if self.gds:
                    GDS.content[level].append(
                        gdsmod.cell_reference([x, y], cell.cell_name, a,
                            flip=flip ^ inode.flip, array=inode.array,
                            mag=inode.scale))
                    #TODO: check: for external gds files in cells named "load_gds",
                    # the array setting needs to be transferred to the gds in this instance
                    # if instantiate is False on "load_gds".
                    if self.infolevel > 2:
                        print("{}  add instance '{}' @ ({:.3f}, {:.3f}, {:.3f})".\
                            format(level*tab, cell.cell_name, x, y, a))
                        print('{}    instance flip = {}'.format(level*tab, inode.flip))
                        print('{}    instance array = {}'.format(level*tab, inode.array))

                if self.nazca:
                    self.CELL.add_instance(inode, [x, y, a], flip)
        return None


    def add_gdsfiles(self, gdsfile_iter, level):
        """Add gds instances to cell.

        Returns:
            None
        """
        for gdsinfo, [x, y, a], flip in gdsfile_iter:
            if self.gds:
                filename, cellname, newcellname, layermap,\
                    cellmap, scale, strm = gdsinfo
                GDS.content[level].append(gdsmod.cell_reference([x, y],
                    newcellname, a, mag=scale, flip=flip))
                if newcellname not in self.gds_files:
                    self.gds_files[newcellname] = \
                        (filename, cellname, layermap, cellmap, scale, strm)
            if self.nazca:
                self.CELL.add_gds()
        return None


    def export_topcells(self, topcells):
        """Export topcells.

        Loop over each topcell in the topcells list.

        Returns:
            None
        """
        # initialize export formats with single output file across topcells.
        bbstr = ''
        if self.bblock:
            bbstr = '.lib'
            bblock = BBlock(filebasename=self.filebasename+bbstr)
        if self.gds:
            GDS.open(filebasename=self.filebasename+bbstr)
        if self.output == 'file':
            mplt.ioff()
            #print('Switched of Matplotlib interactive output.')
        if self.nazca:
            self.CELL = ClsNazca(instantiate=self.instantiate)

        cells_visited = set()
        for topcell in topcells:

            # initialize formats with a new output file per topcell
            if self.plt: #separate plot for each topcell.
                PLT.open(topcell, self.title)
            if self.svg: #separate plot for each topcell.
                SVG.open(self.filename)
            if self.spt:
                SPT.open(self.filename)

            NL = Netlist()
            cell_iter = NL.celltree_iter2(topcell, flat=self.flat,
                cells_visited=cells_visited, infolevel=self.infolevel)
            for params in cell_iter:
                if params.cell_start:
                    if params.cell_create: #take care of formats with cell hierarchy opening
                        if self.gds:
                            GDS.content[params.level] = []
                            GDS.content[params.level].append(
                                gdsmod.cell_open(params.cell.cell_name))
                        if self.nazca:
                            self.CELL.open(params)
                    self.add_polygons(params.iters['polygon'], params.parent_level)
                    self.add_polylines(params.iters['polyline'], params.parent_level)
                    self.add_annotations(params.iters['annotation'], params.parent_level,
                        create=params.cell_create, cell=params.cell)
                    self.add_instances(params.iters['instance'], params.parent_level)
                    self.add_gdsfiles(params.iters['gdsfile'], params.parent_level)
                else:
                    if params.cell_close: # take care of formats with cell hierarchy closing
                        if self.gds:
                            GDS.content[params.level].append(gdsmod.cell_close())
                            GDS.write(params.level)
                            del GDS.content[params.level]
                        if self.nazca:
                            self.CELL.close()

            cells_visited |= NL.cells_visited

            # Close exports per topcell
            if self.plt:
                PLT.close()
                if self.info:
                    print('...plotting plt')
                    if self.output != 'file':
                        mplt.show()
            if self.svg:
                print("...Wrote file '{}'".format(SVG.filename))
                SVG.close()
            if self.spt:
                SPT.close()
            if self.output is not None:
                try:
                    file = os.path.join(self.path, topcell.cell_paramsname) + '.png'
                except:
                    print("Warning: BB has no basename '{}'".format(topcell.cell_name))
                    file = os.path.join(self.path, topcell.cell_name) + '.png'
                print("export BB file: '{}'".format(file))
                mplt.savefig(file)
            if self.bblock:
                bblock.body(topcell)
            if self.show_cells:
                print('Cells processed:', cells_visited)

        # Close exports common across topcells
        if self.bblock:
            bblock.footer()
        if self.gds: # save all external GDS file based instances
            if self.infolevel > 0:
                print('----\nsave external gds instances')
            for newcellname, mapping in self.gds_files.items():
                filename, cellname, layermap, cellmap, scale, strm = mapping
                g = gstream.GDSII_stream(filename, cellmap=cellmap, layermap=layermap)
                stream = g.GDSII_stream_cell(newcellname)
                GDS.outfile.write(stream)
                if self.infolevel > 0:
                    print("{}'{}'".format(tab*2, filename))
            GDS.outfile.write(gdsmod.layout_close())
            GDS.outfile.close()
            print("...Wrote file '{}'".format(GDS.filename))
            if self.ascii:
                ga = gstream.GDSII_stream(GDS.filename)
                ga.ASCII_write(GDS.filename+'.asc')
            if util.isnotebook():
                display(HTML('<pre>...<a href="{0}" target="_blank">{0}</a></pre>'.\
                    format(GDS.filename)))


    def generate_layout(self, topcells=None):
        """Internal wrapper function before exporting the layout.

        Create final topcells list and set final flags for export.

        Args:
            topcells (Cell | list of Cells): Cell(s) to export
                (default = None, which, exports the 'nazca' default gds cell)
            filename (str): gds output filename (default = 'nazca_export.gds')
            ascii (bool): export ascii version of gds (default = False)
            show_cells (bool): print exported cell names to stdout (default = False)
            gds (bool): export gds (default = True)
            flat (bool): export flat gds, i.e. no hierarchy (default = False)
            spt (bool): export spt file (default = False)
            plt (bool): generate matplotlib based layout (default = False)
            clear (bool): clear mask layout between consecutive exports (default = True)
            title (str): title for the layout if the format allows for a title
            output (str): type of output stream (screen, file ...)
            path (str): output dir for saved Matplotlib plots (default = '')

        Returns:
            None
        """
        if isinstance(topcells, str):
            print('WARNING: You are trying to export a string instead of a cell.')
            print('Did you mean: export_gds(filename=\'' + topcells + '\')?')
            return 0

        if self.infolevel > 0:
            print('gds:{}, plt:{}, spt:{}, flat:{}'.format(
                self.gds, self.plt, self.spt, self.flat))
            print('Generate layout:')

        # close cells in reverse order.
        if not self.clear: #keep default topcell open
            for cell in cfg.cells[1:-1:-1]:
                cell.close()
                if not cfg.solve_direct:
                    cfg.cells[0]._solve() #solve (still open) default topcell
        else:
            for cell in cfg.cells[::-1]:
                cell.close()

        # construct a list of topcells
        if topcells is None:
            topcells = []
        if not topcells: # 0 topcells
            topcells = [cfg.defaultcell]
        elif isinstance(topcells, Cell): # 1 topcell
            topcells = [topcells]
        else: # >1 topcells
            topcells = set(topcells + cfg.topcells)

        self.export_topcells(topcells)

        if self.clear: #Create a new default cell.
            clear_layout()
            if self.infolevel > 0:
                print("Recreated topcell '{}'.".format(cfg.defaultcellname))

        return None


def clear_layout():
    """Remove all cell references to start a brand new layout.

    A new topcell 'nazca' will be created.

    Returns:
        None
    """
    for name, cell in cfg.cellnames.items():
        del cell
    cfg.cellnames = {}

    #cfg.cellnames.pop(cfg.defaultcellname) #avoid triggering cell reuse warning
    cfg.defaultcell = Cell(name=cfg.defaultcellname)
    cfg.cp = cfg.defaultcell.org
    return None


def export_clear():
    """Clear the default topcell.

    Note that export_clear does not clear or delete any other cells.

    Returns:
        None
    """
    if cfg.cells:
        for cell in cfg.cells[::-1]:
            cell.close()
        cfg.cellnames.pop(cfg.defaultcellname)
        cfg.defaultcell = Cell(name=cfg.defaultcellname)
        cfg.cp = cfg.defaultcell.org
        #if self.infolevel:
        #    print("Recreated topcell '{}'.".format(cfg.defaultcellname))
    return None


#==============================================================================
#
#==============================================================================
def verify_topcells(topcells):
    #return None
    if topcells is None:
        return None
    elif isinstance(topcells, Cell):
        return None
    elif isinstance(topcells, list):
        for cell in topcells:
            if not isinstance(cell, Cell):
                raise ValueError('You need to provide a Cell or list of Cells.')
        return None
    raise ValueError("You need to provide a Cell object or list of Cell objects"\
        " for keyword 'topcells'. Instead got an object of type {}.".format(type(topcells)))


def rebuild(cell, instantiate=True, flat=False, layermap=None, layermapmode=None):
    """Flatten, rename, relayer, reshape and/or filter a Nazca Cell object.

    The original cell remains unchanged.

    Args:
        cell (Cell): input cell(tree)

    Returns:
        Cell: rebuild input cell
    """
    export = Export_layout()
    export.nazca = True
    export.flat = flat
    export.setlayermap(layermap, layermapmode)
    export.instantiate = instantiate
    export.generate_layout(cell)
    return export.CELL.topcell

def celltree_iter(cell, level=0, position=None, flat=False,
         cells_visited=None, infolevel=0):
    return Netlist().celltree_iter2(cell=cell, position=position, flat=flat,
         cells_visited=cells_visited, infolevel=infolevel)

def export(topcells=None,
        filename=None,
        gds=False,
        ascii=False,
        plt=False,
        svg=False,
        spt=False,
        flat=False,
        infolevel=0,
        show_cells=False,
        info=True, #progress info
        clear=None,
        title=None,
        output=None,
        path='',
        bb=False):
    """Export layout to gds file for all cells in <topcells>.

    Args:
        topcells (Cell | list of Cells): Cell(s) to export
            (default = None, which, exports the 'nazca' default gds cell)
        filename (str): gds output filename (default = 'nazca_export.gds')
        clear (bool): clear mask layout between consecutive exports (default = True)
        gds (bool): export gds (default = True)
        ascii (bool): export ascii version of gds (default = False)
        svg (bool): export spt file (default = False)
        plt (bool): generate matplotlib based layout (default = False)
        spt (bool): export spt file (default = False)
        flat (bool): export flat gds, i.e. no hierarchy (default = False)
        infolevel (int): amount of debug info to stdout (default = 0)
        show_cells (bool): (default = False)
        info (bool): (default = True)
        title (str): title for the layout if the format allows for a title
        output (str): type of output stream (screen, file ...)
        path (str): output dir for saved Matplotlib plots (default = '')
        bb (bool): export as a library bb (default = False)

    Returns:
        None
    """
    verify_topcells(topcells)
    #TODO: implement layermap filter.
    export = Export_layout()
    export.filename = filename
    export.ascii = ascii
    export.infolevel = infolevel
    export.show_cells = show_cells
    export.gds = gds
    export.flat = flat
    export.spt = spt
    export.plt = plt
    if svg and not cfg.SVGWRITE:
        export.svg = False
        print("Warning: could not load module 'svgwrite'. Skipping svg export.")
    else:
        export.svg = svg
    export.info = info
    export.title = title
    export.output = output
    export.path = path
    export.bblock = bb
    if clear is None:
        export.clear = cfg.export_clear
    else:
        export.clear = clear
    if export.info:
           print('Starting layout export...')
    if export.gds:
        if info:
            print('...gds generation')
    if export.spt:
        print('...spt generation')
        if cfg.__mapping_errors:
            print("Warning: bb mappings missing:")
            for text in cfg.__mapping_errors:
                print('  ', text)
    if export.plt:
        if export.info:
            print('...matplotlib generation')
        mplt.show()
    if export.svg:
        if info:
            print('...svg generation')

    export.generate_layout(topcells)
    #print('done') do not want done in the notebook
    return None


def export_svg(topcells=None, title=None, path='', **kwargs):
    export(topcells, plt=False, gds=False, svg=True, info=False,
        title=title, **kwargs)


def export_plt(topcells=None, clear=None, title=None, output=None, path='', **kwargs):
    """Export layout with Matplotlib for all cells in <topcells>.

    Args:
        topcells (Cell | list of Cells): Cell(s) to export
            (default = None, which exports the 'nazca' default gds cell)
        clear (bool): clear mask layout between consecutive exports (default = True)
        title (str): title for the layout if the format allows for a title
        output (str): Matplotlib output stream (screen, file ...) (default = None -> screen)
        path (str): output dir for saved Matplotlib plots (default = '')

    Returns:
        None
    """
    export(topcells, plt=True, gds=False, info=False, clear=clear, title=title,
        output=output, path=path, **kwargs)


def export_spt(topcells=None, filename=None):
    export(topcells, filename=filename,
        plt=False, gds=False, spt=True, info=True)


def export_gds(topcells=None, filename=None, flat=False, spt=False, clear=None,
        bb=False, **kwargs):
    """Export layout to gds file for all cells in <topcells>.

    Args:
        topcells (Cell | list of Cells): Cell(s) to export
            (default = None, which exports the 'nazca' default gds cell)
        filename (str): gds output filename (default = 'nazca_export.gds')
        flat (bool): export flat gds, i.e. no hierarchy (default = False)
        clear (bool): clear mask layout between consecutive exports (default = True)
        bb (bool): Export design as a building block (default = False)'

    Returns:
        None
    """
    verify_topcells(topcells)
    export(topcells, filename=filename, plt=False, gds=True, spt=spt,
        info=True, flat=flat, clear=clear, bb=bb, **kwargs)


def layout(layermap=None, layermapmode='all'):
    """Create a layout object for rebuilding cells.'

    Returns:
        Export_layout: layout object
    """
    ly = Export_layout(layermap=layermap, layermapmode=layermapmode)
    ly.nazca = True
    ly.CELL = ClsNazca(infolevel=0)
    return ly



