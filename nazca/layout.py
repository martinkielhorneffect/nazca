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
import numpy as np
import pandas as pd
import os.path

import matplotlib.pyplot as mplt
from matplotlib.patches import Polygon as matPolygon
from matplotlib.collections import PatchCollection

import svgwrite

from . import gds
from . import gds_base as gbase
from . import gds_import as gstream
from .netlist import Netlist, Cell, Pointer, Polygon
#from .pdk_template_core import strip_groupname
from . import cfg
from .util import isnotebook
from IPython.core.display import display, HTML
from pprint import pprint
from . import util


#==============================================================================
# spt
#==============================================================================
class SPT():
    """Handle spt compatible export of masks."""

    def __init__(self):
        self.Fspt = None
        self.bb_cnt = 0

    def open(self, filename=None):
        """Initialize spt file export."""
        self.Fspt = None # spt-file handle
        self.spt_filename = os.path.splitext(filename)[0] + '.spt'
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
        self.bb_cnt += 1

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
            format(BBname, port, x, y, a,', '.join(varlist), self.bb_cnt))

    def close(self):
        """Close spt output routine."""
        if self.Fspt is not None:
            self.Fspt.close()
            print("...Wrote file '{}'".format(self.spt_filename))

Spt = SPT()


#==============================================================================
# svg
#==============================================================================
class SVG():
    """Handle svg compatible export of masks."""

    def __init__(self, filename='nazca_export.svg'):
        #if filename is not None:
        #    filename = os.path.splitext(filename)[0] + '.svg'
        self.drawing = svgwrite.Drawing(filename=filename, size=('100%', '100%'),
            profile='tiny', debug=False)
        self.minx = 1e6
        self.miny = 1e6
        self.maxx = -1e6
        self.maxy = -1e6

    def open(self):
        """Open svg.

        Open after layer colors have been loaded.

        Returns:
            None
        """
        # create a map of defined layer colors to speed up lookup
        self.layermap = defaultdict(list)
        if not cfg.colors.empty:
            for i, row in cfg.colors.loc[:, ['layer', 'datatype']].iterrows():
                ML = (int(row[0]), int(row[1]))
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

Svg = SVG()


#==============================================================================
# matplotlib export
#==============================================================================
class Matplot():
    """Handle matplotlib export of masks."""

    def __init__(self):
        """Construct an object handling Matplotlib exports."""

        try:
            font = cfg.matplotlib_font
        except:
            font = {
                #'family'  : 'normal',
                'style'   : 'normal',
                'weight' : 'light', #'light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'
                'size'   : 24}
        mplt.rc('font', **font)


    def open(self, cell=None, title=None):
        """Inititialize Matplotlib mask output."""
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

        # create a map of defined layer colors to speed up lookup
        self.layermap = defaultdict(list)
        if not cfg.colors.empty:
            for i, row in cfg.colors.loc[:, ['layer', 'datatype']].iterrows():
                ML = (int(row[0]), int(row[1]))
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
        if layer == cfg.documentation_pin_layer:
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
        if cfg.pin_settings['bb_pin_layer'] != cfg.documentation_pin_layer:
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
matplot = Matplot()


tab = '  '
class Export_layout():

    def __init__(self):
        self.gds_files = dict() # existing gds files used in the layout

    def add_polygons(self, netlist, cnode, trans, flip, infolevel=0):
        """Add polygons content to cell.

        Returns:
            bytestring: gds stream
        """
        content = []
        s = 1
        if flip:
            s = -1
        for org, pgon in netlist.polygon_iter(cnode, infolevel=infolevel):
            if flip:
                org = org.copy()
                org.flip()
            [x0, y0, a0] = org.multiply_ptr(trans).get_xya()
            a = s*np.radians(a0)
            xy = [(x0+cos(a)*u-sin(a)*v, y0+s*(sin(a)*u+cos(a)*v))
                for u, v in pgon.points]
            if self.gds:
                content.append(gbase.gds_polygon(xy, lay=int(pgon.layer[0]),
                    datatype=int(pgon.layer[1])))
            if self.plt or self.svg:
                x1, y1, x2, y2 = pgon.bbox
                bbox = ( x0+cos(a)*x1-sin(a)*y1, y0+s*(sin(a)*x1+cos(a)*y1),
                         x0+cos(a)*x2-sin(a)*y2, y0+s*(sin(a)*x2+cos(a)*y2) )
                if self.plt:
                    matplot.add_polygon(pgon.layer, xy, bbox)
                if self.svg:
                    Svg.add_polygon(pgon.layer, xy, bbox)

        return b''.join(content)


    def add_polylines(self, netlist, cnode, trans, flip, infolevel=0):
        """Add polylines content to cell.

        Returns:
            bytestring: gds stream
        """
        content = []
        s = 1
        if flip:
            s = -1
        for org, pline in netlist.polyline_iter(cnode, infolevel=infolevel):
            if flip:
                org = org.copy()
                org.flip()
            [x0, y0, a0] = org.multiply_ptr(trans).get_xya()
            a = s*np.radians(a0)
            xy = [(x0+cos(a)*u-sin(a)*v, y0+s*(sin(a)*u+cos(a)*v))
                for u, v in pline.points]
            if self.gds:
                content.append(gbase.gds_polyline(xy, pline.width,
                    lay=int(pline.layer[0]), datatype=int(pline.layer[1]),
                    pathtype=pline.pathtype))
            if self.plt or self.svg:
                x1, y1, x2, y2 = pline.bbox
                bbox = ( x0+cos(a)*x1-sin(a)*y1, y0+s*(sin(a)*x1+cos(a)*y1),
                         x0+cos(a)*x2-sin(a)*y2, y0+s*(sin(a)*x2+cos(a)*y2) )
                polygon_xy = util.polyline2polygon(xy, width=pline.width, miter=0.5)
                if self.plt:
                    matplot.add_polygon(pline.layer, polygon_xy, bbox)
                    #matplot.add_polyline(pline.layer, xy, pline.width, bbox)
                if self.svg:
                    Svg.add_polygon(pline.layer, polygon_xy, bbox)
        return b''.join(content)


    def add_annotations(self, netlist, cnode, trans, flip, infolevel=0):
        """Add annotation content to cell.

        Returns:
            bytestring: gds stream
        """
        content = []
        for org, anno in netlist.annotation_iter(cnode, infolevel=infolevel):
            if flip:
                org = org.copy()
                org.flip()
            x0, y0 = org.multiply_ptr(trans).get_xy()
            if self.gds:
                content.append(gbase.gds_annotation(xy=[x0, y0], string=anno.text,
                    lay=int(anno.layer[0]), datatype=int(anno.layer[1])))
        return b''.join(content)


    def add_gdsfiles(self, netlist, cnode, trans, flip, infolevel=0):
        """Add gds instances to cell.

        Returns:
            bytestring: gds stream
        """
        content = []
        gdsfiles = netlist.gdsfile_iter(cnode, infolevel=infolevel)
        for org, filename, cellname, newcellname, layermap, cellmap, scale,\
                strm in gdsfiles:
            [x0, y0, a0] = org.multiply_ptr(trans).get_xya()
            if self.gds:
                content.append(gds.cell_reference([x0, y0], newcellname, a0, mag=scale, flip=flip))
                if newcellname not in self.gds_files:
                    self.gds_files[newcellname] = (filename, cellname,
                        layermap, cellmap, scale, strm)
        return b''.join(content)


    def add_instances(self, instance_iter, level, infolevel=0):
        """Add instances to cell.

        Returns:
            bytestring: gds stream
        """
        content = []
        for instance in instance_iter:
            cnode = instance #, cnode_iter, position = instance
            cellname = cnode.cell.cell_name
            org = cnode.pointer.copy()
            if self.flat:
                try:
                    cnode.cell.BBmap[3]
                except:
                    continue
                if self.spt:
                    if self.fliplist[-1]:
                        org.flip()
                    [x0, y0, a0] = org.multiply_ptr(self.translist[-1]).get_xya()
                    Spt.add_line(cnode.cell, x0, y0, a0, self.fliplist[-1])
                    continue

            if self.gds:
                if not cnode.cell.instantiate:
                     continue
                [x0, y0, a0] = org.multiply_ptr(self.translist[-1]).get_xya()
                content.append(gds.cell_reference([x0, y0], cellname, a0,
                    flip=self.fliplist[-1] ^ cnode.flip, array=cnode.array))
                if self.infolevel > 1:
                    print("{}  add instance '{}' @ ({:.3f}, {:.3f}, {:.3f})".format(
                        level*tab, cellname, x0, y0, a0))
                    print('{}    instance flip = {}'.format(level*tab, cnode.flip))
                    print('{}    instance array = {}'.format(level*tab, cnode.array))
        return b''.join(content)


    def close_cells(self, start_level, stop_level):
        """Close gds "export levels" up to the last branching."""
        level = start_level
        while level >= stop_level:
            if self.translist:
                self.translist.pop()
                self.fliplist.pop()

            #close gds levels:
            if level in self.export_levels:
                self.export_levels.pop()
                if self.gds:
                    if self.infolevel > 0:
                        print('{}<< close gds cell'.format(level*tab))
                    self.cell_content[level].append(gds.cell_close())
                    self.gds_outfile.write(b''.join(self.cell_content[level]))
                    del self.cell_content[level]
            level -= 1


    fliplast = False
    def process_cell(self, netlist, cnode, orgtrans, level, instance, flip):
        """Process a cell (via its cnode) and generate content from it.

        Returns:
            None
        """
        global fliplast
        cellname = cnode.cell.cell_name

        if instance or level == 0: #create new export level
            self.export_levels.append(level)
            self.translist.append(Pointer(0, 0, 0))
            self.fliplist.append(False)
            fliplast = flip
            self.cell_content[level] = []
            if self.gds:
                if self.infolevel > 0:
                    print("{}>> open gds cell '{}'".format(level*tab, cellname))
                self.cell_content[level].append(gds.cell_open(cellname))

        else: #flattened cell -> content translates into a parent cell.
            translate = self.translist[-1].copy()
            if self.fliplist[-1]:
                orgtrans = orgtrans.copy()
                orgtrans.flip()
            fliplast = self.fliplist[-1] ^ flip
            self.fliplist.append(fliplast)
            translate.move_ptr(orgtrans)
            self.translist.append(translate)

        #generate cell content:
        parent_level = self.export_levels[-1]
        values = (netlist, cnode, self.translist[-1], self.fliplist[-1], self.infolevel)
        self.cell_content[parent_level].append(self.add_polygons(*values))
        self.cell_content[parent_level].append(self.add_polylines(*values))
        self.cell_content[parent_level].append(self.add_annotations(*values))
        self.cell_content[parent_level].append(self.add_gdsfiles(*values))
        self.cell_content[parent_level].append(self.add_instances(cnode.cnode_nb_iter(), level))

        return None


    def export_topcells(self, topcells, show_cells):
        """Export topcells.

        Loop over each topcell in the topcells list.

        Returns:
            None
        """
        if self.output == 'file':
            mplt.ioff()
            #print('Switched of Matplotlib interactive output.')

        if self.bb:
            #bbfile open and header
            pyfilename = os.path.join(self.fln[:-4])+'.py'
            bbfile = open(pyfilename, 'w')
            bbfile_header =\
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

            bbfile.write(bbfile_header)
            bbfile_footer = ""
            bbfile_x = 0

        cells_visited = set()
        for topcell in topcells:
            if self.plt: #separate plot for each topcell.
                matplot.open(topcell, self.title)
            if self.svg: #separate plot for each topcell.
                Svg.open()
            self.export_levels = [] # netlist levels corrresponding to export levels.
            self.cell_content = dict() # collect cell content per level.
            self.translist = []
            self.fliplist = []
            level_last = 0 # top level is 0
            level_duplicate = 100 # level

            netlist = Netlist()
            cells = netlist.celltree_iter(topcell.cnode, flat=self.flat,
                infolevel=0)
            for cnode, level, org, flip in cells:
                if level > level_duplicate: #skip cells that are copies of previous cells.
                    continue
                else:
                    level_duplicate = 100

                cellname = cnode.cell.cell_name
                self.close_cells(start_level=level_last, stop_level=max(1, level))

                if self.infolevel > 0:
                    info = ("'{}' @ ({:.3f}, {:.3f}, {:.3f}), inst={}, flip={}".\
                        format(cellname, org.x, org.y, org.a,
                            cnode.cell.instantiate, flip))

                if not self.flat:
                    inst = cnode.cell.instantiate
                else:
                    inst = False

                visited = cellname in cells_visited
                if not inst or not visited: #new cell or flat:
                    if self.infolevel > 0:
                        if not visited:
                            print('{}create cell {}'.format('-'*len(tab)*level, info))
                        else:
                            print('{}re-use cell {}'.format('-'*len(tab)*level, info))
                    cells_visited.add(cellname)
                    self.process_cell(netlist, cnode, org, level, inst, flip)
                    level_last = level
                elif inst: #copy of earlier cell:
                    if self.infolevel > 0:
                        print('{}reuse cell {}'.format('-'*len(tab)*level, info))
                    level_duplicate = level
                    level_last = level-1

            self.close_cells(start_level=level, stop_level=0)

            if show_cells:
                print('Cells processed:', cells_visited)

            if self.plt:
                matplot.close()
            if self.output is not None:
                try:
                    file = os.path.join(self.path, topcell.cell_paramsname) + '.png'
                except:
                    print("Warning: BB has no basename '{}'".format(topcell.cell_name))
                    file = os.path.join(self.path, topcell.cell_name) + '.png'
                print("export BB file: '{}'".format(file))
                mplt.savefig(file)

            if self.bb:
                #bbfile body
                gdsfilename = os.path.join(self.fln)
                cellname = topcell.cell_paramsname
                cellvar = cellname.replace('.', '_')
                bbfile.write("\nwith Cell(name='BB', autobbox=False, instantiate=False) as {0}:\n".\
                    format(cellvar))
                bbfile.write("    load_gds(filename='{}', cellname='{}').put(0)\n".\
                    format(gdsfilename, cellname))
                bbfile_footer += ("    {}.put('org', {:.1f})\n".format(cellvar, bbfile_x))

                if topcell.length is not None:
                    bbfile_x += 1.1*topcell.length
                for name, pin in sorted(topcell.pin.items()):
                    if name == 'org':
                        continue
                    x, y, a = pin.pointer.xya()
                    if pin.xs is None:
                        xs = None
                    else:
                        xs = "'{}'".format(pin.xs)
                    bbfile.write("    Pin(name='{}', width={}, xs={}).put({:.5f}, {:.5f}, {:.5f})\n".\
                        format(name, pin.width, xs, x, y, a))

        if self.bb:
            #bbfile footer
            bbfile.write("\nif __name__ == '__main__':\n")
            bbfile.write(bbfile_footer)
            bbfile.write("    export_gds()\n")

        # Save all external GDS file based instances:
        if self.gds:
            if self.infolevel > 0:
                print('----\nsave external gds instances')
            for newcellname, mapping in self.gds_files.items():
                filename, cellname, layermap, cellmap, scale, strm = mapping
                g = gstream.GDSII_stream(filename, cellmap=cellmap, layermap=layermap)
                stream = g.GDSII_stream_cell(newcellname)
                self.gds_outfile.write(stream)
                if self.infolevel > 0:
                    print("{}'{}'".format(tab*2, filename))



    def generate_layout(self,
        topcells=None,
        filename='nazca_export.gds',
        do_ascii=False,
        infolevel=False,
        show_cells=False,
        gds_=True,
        flat=False,
        spt=False,
        plt=False,
        svg=False,
        clear=False,
        title=None,
        output=None,
        path='',
        bb=False):
        """Internal wrapper function before exporting the layout.

        Create final topcells list and set final flags for export.

        Args:
            topcells (Cell | list of Cells): Cell(s) to export
                (default = None, which, exports the 'nazca' default gds cell)
            filename (str): gds output filename (default = 'nazca_export.gds')
            do_ascii (bool): export ascii version of gds (default = False)
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
        self.fln = filename
        self.gds = gds_
        if do_ascii:
            self.gds = True
        self.flat = flat
        self.spt = spt
        self.plt = plt
        self.svg = svg
        self.Fspt = Spt.Fspt
        self.title = title
        self.output = output
        self.path = path
        self.infolevel = infolevel
        self.bb = bb

        if isinstance(topcells, str):
            print('WARNING: You are trying to export a string instead of a cell.')
            print('Did you mean: export_gds(filename=\'' + topcells + '\')?')
            return 0

        if not isinstance(filename, str):
            print("WARNING: Filename provided is not a string but a {}."\
                  " Using filename 'export_nazca.gds' instead.".format(type(self.fln)))
            self.fln = 'export_nazca.gds'

        if self.infolevel > 0:
            print('gds:{}, plt:{}, spt:{}, flat:{}'.format(
                self.gds, self.plt, self.spt, self.flat))
            print('Generate layout:')

        # close cells in reverse order.
        if not clear: #keep default topcell open
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

        if self.gds:
            self.gds_outfile = open(self.fln, 'wb')
            self.gds_outfile.write(gds.layout_open())
            self.export_topcells(topcells, show_cells)
            self.gds_outfile.write(gds.layout_close())
            self.gds_outfile.close()
            print("...Wrote file '{}'".format(self.fln))
            if do_ascii:
                ga = gstream.GDSII_stream(self.fln)
                ga.ASCII_write(self.fln+'.asc')
            if isnotebook():
                display(HTML('<pre>...<a href="{0}" target="_blank">{0}</a></pre>'.format(self.fln)))
        else:
           self.export_topcells(topcells, show_cells)

        if clear: #Create a new default cell.
            clear_layout()
            if self.infolevel > 0:
                print("Recreated topcell '{}'.".format(cfg.defaultcellname))

        return None


def clear_layout():
    """Remove all cell references to start a new layout.

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
    """

    if cfg.cells:
        for cell in cfg.cells[::-1]:
            cell.close()
        cfg.cellnames.pop(cfg.defaultcellname)
        cfg.defaultcell = Cell(name=cfg.defaultcellname)
        cfg.cp = cfg.defaultcell.org
        #if self.infolevel:
        #    print("Recreated topcell '{}'.".format(cfg.defaultcellname))

#==============================================================================
#
#==============================================================================
def verify_topcells(topcells):
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


def export(topcells=None,
           filename=None,
           do_ascii=False,
           infolevel=0,
           show_cells=False,
           gds=True,
           flat=False,
           spt=False,
           plt=False,
           svg=False,
           info=True,
           clear=True,
           title=None,
           output=None,
           path='',
           bb=False):
    """Export layout to gds file for all cells in <topcells>.

    Args:
        topcells (Cell | list of Cells): Cell(s) to export
            (default = None, which, exports the 'nazca' default gds cell)
        filename (str): gds output filename (default = 'nazca_export.gds')
        do_ascii (bool): export ascii version of gds (default = False)
        show_cells (bool): (default = False)
        gds (bool): export gds (default = True)
        flat (bool): export flat gds, i.e. no hierarchy (default = False)
        spt (bool): export spt file (default = False)
        plt (bool): generate matplotlib based layout (default = False)
        info (bool): (default = True)
        clear (bool): clear mask layout between consecutive exports (default = True)
        title (str): title for the layout if the format allows for a title
        output (str): type of output stream (screen, file ...)
        path (str): output dir for saved Matplotlib plots (default = '')

    Returns:
        None
    """

    verify_topcells(topcells)
    #TODO: implement layermap filter.
    export = Export_layout()

    if filename is None:
        filename = 'nazca_export.gds'
    if info:
           print('Starting layout export...')
    if gds:
        if info:
            print('...gds generation')
    if spt:
        flat = True
        print('...spt generation')
        if cfg.__mapping_errors:
            print("Warning: bb mappings missing:")
            for text in cfg.__mapping_errors:
                print('  ', text)
        Spt.open(filename)
    if plt:
        flat = True
        if info:
            print('...matplotlib generation')
    if svg:
        flat = True
        if info:
            print('...svg generation')

    export.generate_layout(topcells, filename, do_ascii, infolevel=infolevel,
        gds_=gds, flat=flat, spt=spt, plt=plt, svg=svg, clear=clear, title=title,
        output=output, path=path, bb=bb)

    if spt:
        Spt.close()
    if plt:
        if info:
            print('...plotting plt')
        if output != 'file':
            mplt.show()
    if svg:
        Svg.close()
    if gds:
        pass
        #print('Done.')
    return None


def export_svg(topcells=None, title=None, path='', **kwargs):
    verify_topcells(topcells)
    export(topcells, plt=False, gds=False, svg=True, info=False,
        title=title, **kwargs)


def export_plt(topcells=None, clear=True, title=None, output=None, path='', **kwargs):
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
    verify_topcells(topcells)
    export(topcells, plt=True, gds=False, info=False, clear=clear, title=title,
        output=output, path=path, **kwargs)


def export_spt(topcells=None, filename=None):
    verify_topcells(topcells)
    export(topcells, filename=filename,
        plt=False, gds=False, spt=True, info=True)


def export_gds(topcells=None, filename=None, flat=False, spt=False, clear=True,
        bb=False, **kwargs):
    """Export layout to gds file for all cells in <topcells>.

    Args:
        topcells (Cell | list of Cells): Cell(s) to export
            (default = None, which exports the 'nazca' default gds cell)
        filename (str): gds output filename (default = 'nazca_export.gds')
        flat (bool): export flat gds, i.e. no hierarchy (default = False)
        clear (bool): clear mask layout between consecutive exports (default = True)

    Returns:
        None
    """
    verify_topcells(topcells)
    export(topcells, filename=filename, plt=False, gds=True, spt=spt,
        info=True, flat=flat, clear=clear, bb=bb, **kwargs)
