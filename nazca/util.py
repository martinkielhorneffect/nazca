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
# Utility routines
#
# (c) 2016-2017  Xaveer Leijtens, Ronald Broeke
#
from math import hypot
from collections import OrderedDict
from . import gds_base as gbase
from .netlist import Cell
from PIL import Image
import os
import numpy as np
import hashlib
from nazca.netlist import Polygon

__all__ = ['parameters_to_string', 'string_to_parameters',
        'get_cell_annotation', 'get_cell_polyline', 'get_cell_polygon',
        'md5', 'isnotebook', 'image']


def parameters_to_string(param):
    """Create a string from a parameter dictionary.

    param (dict): (parameter_name, value)

    Returns:
        str: parameters as a string
    """
    plist = ['Parameters:']
    for key, value in param.items():
        plist.append("{} = {}".format(key, value))
    return '\n'.join(plist)


def string_to_parameters(string):
    """Convert a string to a parameter dictionaty.

    Returns:
        ordered dict: <parameter_name: parameter_value>
    """
    lines = string.split('\n')
    p = OrderedDict()
    if (lines[0] == 'Parameters:'):
        for line in lines[1:]:
            param = line.split('=', 1)
            if len(param) == 2:
                p[param[0].strip()] = param[1].strip()
    return p


def get_cell_annotation(cell, convert=False):
    """Yield the <cell>'s annotations one by one.

    If not convert then return XY as integers, as in GDSII file (nm).
    If convert then return XY as floats (um).

    cell (str): GDS cell name
    convert (bool): default convert=False

    Yields:
        int, tuple(int|float, int|float), str:
            annotation layer, position, text
    """
    if convert:
        conv = gbase.gds_db_user
    else:
        conv = 1
    for e in cell.elements:
        if e.etype == gbase.gds_record.TEXT:
            lay, pos, text = e.annotation
            pos[0] *= conv
            pos[1] *= conv
            yield lay, pos, text


def get_cell_polyline(cell, convert=False):
    """Yield the <cell>'s polylines one by one.

    If not convert then return XY as integers, as in GDSII file (nm).
    If convert then return XY as floats (um).

    cell (str): GDS cell name
    convert (bool): default convert=False

    Yields:
        int, tuple(int, int) or (float, float): layer, XY
    """
    if convert:
        conv = gbase.gds_db_user
    else:
        conv = 1
    for e in cell.elements:
        if e.etype == gbase.gds_record.PATH:
            lay, points = e.polyline
            XY = []
            for i in range(0, len(points), 2):
                XY.append((points[i] * conv, points[i+1] * conv))
            yield lay, XY

def get_cell_polygon(cell, convert=False):
    """Yield the <cell>'s polygones one by one.

    If not convert then return XY as integers, as in GDSII file (nm).
    If convert then return XY as floats (um).

    cell (str): GDS cell name
    convert (bool): default convert=False

    Yields:
        int, tuple(int, int) or (float, float): layer, XY
    """
    if convert:
        conv = gbase.gds_db_user
    else:
        conv = 1
    for e in cell.elements:
        if e.etype == gbase.gds_record.BOUNDARY:
            lay, points = e.polygon
            XY = []
            for i in range(0, len(points), 2):
                XY.append((points[i] * conv, points[i+1] * conv))
            yield lay, XY


def md5(s, N=5):
    """Return first N characters of md5 hash of string s."""
    return hashlib.md5(s.encode()).hexdigest()[:N]


def isnotebook():
    """Check if code is run in a Jupyter notebook.

    Returns:
        bool: wether call is made from a notebook (True)
    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':  # Jupyter notebook or qtconsole?
            return True
        elif shell == 'TerminalInteractiveShell':  # Terminal running IPython?
            return False
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

def PIL2array(img):
    return np.array(img.getdata(),
            np.bool).reshape(img.size[1], img.size[0])

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
            # or
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
    pix = PIL2array(bw)
    width, height = bw.size
    width_tot = width*p + 2*box_buf
    height_tot = height*p + 2*box_buf
    print('Generating {}x{} pixels image of {:.0f}x{:.0f} um2, edge is {} um.'.\
         format(width, height, width*p, height*p, box_buf))
    h0 = halign[align[0]] * width_tot
    v0 = valign[align[1]] * height_tot

    with Cell(cellname) as C:
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
                        Polygon(layer=layer, points=xy).\
                            put(h0+box_buf, v0+box_buf, 0)
                        x0 = x1
                        lb = 0
            if lb > 0:
                x1 = x0 + lb * p
                xy = [(x0, y0), (x1, y0), (x1, y0-p), (x0, y0-p)]
                Polygon(layer=layer, points=xy).\
                    put(h0+box_buf, v0+box_buf, 0)

        if box_layer is not None:
            Polygon(layer=box_layer, points=[(0, 0), (0, height_tot),
                (width_tot, height_tot), (width_tot, 0)]).put(h0, v0, 0)
    return C


def klayout2nazca(string):
    """Convert a copy-pasted path (x,y points) from Klayout into list of points.

    Args:
        string (str): string with points to convert

    Returns:
        list of points (x,y)

    Example::

        klayout2nazca("1.0, 3.5, 2.0, 5.5")

        out: [(1.0, 3.5), (2.0, 5.5)]
    """
    return [[float(x) for x in line.split('\t')]
        for line in string.strip().split('\n')]

def pointdxdy(xy, i, w):
    """Helper function to return the left and right point to the starting
    point of a line segment from (x1, y1) to (x2, y2) with width w. The
    line may contain many points, but only the first two are used.

    Args:
        xy: list of points (x,y)
        i (int): index of point in list
        w (float): width of line segment

    Returns:
        point (dx, dy)
    """
    dx = xy[i+1][1] - xy[i+0][1]
    dy = xy[i+0][0] - xy[i+1][0]
    l = hypot(dx, dy)
    dx *= w / 2 / l
    dy *= w / 2 / l
    return [dx, dy]


def corner(xy, i, dxy):
    """Helper function to return the coordinates of 4 corner points of a line
    segment with a given distance in x and y.

    Args:
        xy: list of points (x,y)
        i (int): index of point in list
        dxy: point (dx, dy)

    Returns:
        tuple of 4 points (x,y) which are the corner points of the line
        segment.
    """
    return ((xy[i][0]   + dxy[0], xy[i][1]   + dxy[1]),
            (xy[i][0]   - dxy[0], xy[i][1]   - dxy[1]),
            (xy[i+1][0] + dxy[0], xy[i+1][1] + dxy[1]),
            (xy[i+1][0] - dxy[0], xy[i+1][1] - dxy[1]),)


# Intersection point (xi, yi) of two lines that go through
# (x0, y0), (x1, y1) and (x2, y2), (x3, y3).
# Fails for parallel lines (divide by zero), but the algorithm below
# ensures this is not the case.
def intersect(xy0, xy1, xy2, xy3): # four points
    """Helper function to intersect two lines. Intersection point (xi, yi)
    of two lines that go through (x0, y0), (x1, y1) and (x2, y2), (x3, y3).
    Fails for parallel lines (divide by zero), but the algorithm in
    polyline2polygon ensures this is not the case.

    Args:
        list of four points (x,y)

    Returns:
        intersection point (xi,yi)
    """
    x0, y0 = xy0
    x1, y1 = xy1
    x2, y2 = xy2
    x3, y3 = xy3
    D  = (x3-x2)*(y1-y0)-(x1-x0)*(y3-y2)
    Dx = (x3-x2)*(y2-y0)-(x2-x0)*(y3-y2)
    xi = x0 + Dx/D * (x1-x0)
    yi = y0 + Dx/D * (y1-y0)
    return (xi, yi)


def polyline2polygon(xy, w):
    """Return a polygon that contains the outline points of a polyline with
    given width.
    Since we have to specify the outline of two or more segments that make
    an angle, we have to know what to do with the gap between those
    segments at the outside of the corner. In order to determine which
    points are on the outside of the corner we use the following algorithm:
    Given a line segment between P0 (x0,y0) and P1 (x1,y1), another point P
    (x,y) has the following relationship to the line segment. Compute
    (y - y0) (x1 - x0) - (x - x0) (y1 - y0). If it is less than 0 then P is
    to the right of the line segment, if greater than 0 it is to the left,
    if equal to 0 then it lies on the line segment.
    The routine fills an array from the start with the anticlockwise points
    and from the end with the clockwise points.

    Args:
        xy (list): list of (x,y) points that hold the polygon
        w (float): width of the polyline

    Returns:
        the polygon
    """

    # the fraction of the width of the line segments that is used to
    # determine if a single point is sufficient to describe the outline, or
    # that two points are needed (miter limit).
    fracwidth = 0.5
    dsqrmax = (fracwidth * w)**2
    n = len(xy)
    if n < 2:
        raise ValueError("Polyline2polygon: need at least 2 points for polyline.")
    # Start with the first two points
    dxy1 = pointdxdy(xy, 0, w)
    cxy1 = corner(xy, 0, dxy1)
    xy_start = [cxy1[0]]
    xy_end = [cxy1[1]]
    for i in range(1, n-1): # loop over the points in the polyline
        dxy0 = dxy1
        # Shift corner points from next to current segment.
        cxy0 = cxy1
        # Get corner points for next segment.
        dxy1 = pointdxdy(xy, i, w)
        cxy1 = corner(xy, i, dxy1)
        # left or right turn
        lrt = (xy[i+1][1]-xy[i-1][1]) * (xy[i][0]-xy[i-1][0]) -\
              (xy[i+1][0]-xy[i-1][0]) * (xy[i][1]-xy[i-1][1])
        # Distance (squared) between the two points at the kink (2->4 == 3->5)
        dsqr = (cxy1[0][0] - cxy0[2][0])**2 + (cxy1[0][1] - cxy0[2][1])**2
        # the inside corner point is on the intersection of the two inside
        # lines.
        if lrt > 0: # Left turn
            # Outside corner: use two points, unless these points are close.
            if dsqr < dsqrmax:
                xy_start.append(intersect(cxy0[0], cxy0[2], cxy1[0], cxy1[2]))
            else:
                xy_start.append((xy[i][0]+dxy0[0], xy[i][1]+dxy0[1]))
                xy_start.append((xy[i][0]+dxy1[0], xy[i][1]+dxy1[1]))
            # Inside corner: always intersect.
            xy_end.append(intersect(cxy0[1], cxy0[3], cxy1[1], cxy1[3]))
        elif lrt < 0: # Right turn
            # Outside corner: use two points, unless these points are close.
            if dsqr < dsqrmax:
                xy_end.append(intersect(cxy0[1], cxy0[3], cxy1[1], cxy1[3]))
            else:
                xy_end.append((xy[i][0]-dxy0[0], xy[i][1]-dxy0[1]))
                xy_end.append((xy[i][0]-dxy1[0], xy[i][1]-dxy1[1]))
            # Inside corner: always intersect.
            xy_start.append(intersect(cxy0[0], cxy0[2], cxy1[0], cxy1[2]))
        else:
            continue # No turn: goto next point
    # Last two points.
    xy_start.append(cxy1[2])
    xy_end.append(cxy1[3])
    return xy_start + list(reversed(xy_end))
