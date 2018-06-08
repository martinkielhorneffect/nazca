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
# (c) 2016-2018  Xaveer Leijtens, Ronald Broeke
#
from math import hypot, sqrt
from collections import OrderedDict
from . import gds_base as gbase
import numpy as np
import hashlib
from math import radians, cos, sin


__all__ = ['parameters_to_string', 'string_to_parameters',
        'get_cell_annotation', 'get_cell_polyline', 'get_cell_polygon',
        'make_iter', 'md5', 'isnotebook']


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


def string_to_parameters(string):
    """Convert a string to a parameter dictionary.

    The returned parameter values are represented as type str.

    Expected format of <string>:

    "parameters:
    <parameter> = <value>
    <parameter> = <value>
    ..."

    Header 'parameters:' is case incensitive and spaces will be stripped.

    Args:
        string (str): parameters

    Returns:
        OrderedDict: {<parameter_name>: <parameter_value>}
    """
    lines = string.split('\n')
    p = OrderedDict()
    if (lines[0].lower() == 'parameters:'):
        for line in lines[1:]:
            param = line.split('=', 1)
            if len(param) == 2:
                p[param[0].strip()] = param[1].strip()
            else:
                print("Warning: string_to_parameter: "
                    "Expected one keyword and one value, but found instead: {}\n"\
                    "Provided string: {}".\
                    format(param, string))
    else:
        print("Error: string_to_parameter: "
            "Expected header 'parameters:'\n"
            "Provided string: {}".format(string))
    return p


def get_cell_annotation(cell, convert=False):
    """Yield the <cell>'s annotations one by one.

    If convert is False then return XY as integers, as in GDSII file (nm).
    If convert is True then return XY as floats (um).

    cell (str): GDS cell name
    convert (bool): default convert=False

    Yields:
        int, (int, int) | (float, float): annotation layer, position, text
    """
    if convert:
        conv = gbase.gds_db_user
    else:
        conv = 1
    for e in cell.elements:
        if e.etype == gbase.GDS_record.TEXT:
            lay, pos, text = e.annotation
            pos[0] *= conv
            pos[1] *= conv
            yield lay, pos, text


def get_cell_polyline(cell, convert=False):
    """Yield the <cell>'s polylines one by one.

    If convert is False then return XY as integers, as in GDSII file (nm).
    If convert is True then return XY as floats (um).

    cell (str): GDS cell name
    convert (bool):  convert polyline's values to float (default = False)

    Yields:
        int, (int, int) | (float, float): layer, XY
    """
    if convert:
        conv = gbase.gds_db_user
    else:
        conv = 1
    for e in cell.elements:
        if e.etype == gbase.GDS_record.PATH:
            lay, points = e.polyline
            XY = []
            for i in range(0, len(points), 2):
                XY.append((points[i] * conv, points[i+1] * conv))
            yield lay, XY


def get_cell_polygon(cell, convert=False):
    """Yield the <cell>'s polygons one by one.

    If convert is False then return XY as integers, as in GDSII file (nm).
    If convert is True then return XY as floats (um).

    cell (str): GDS cell name
    convert (bool): convert polygon's values to float (default = False)

    Yields:
        int, (int, int) | (float, float): layer, XY
    """
    if convert:
        conv = gbase.gds_db_user
    else:
        conv = 1
    for e in cell.elements:
        if e.etype == gbase.GDS_record.BOUNDARY:
            lay, points = e.polygon
            XY = []
            for i in range(0, len(points), 2):
                XY.append((points[i] * conv, points[i+1] * conv))
            yield lay, XY


def make_iter(x):
    """Return x as tuple, if x is not a string and not iterable."""
    if x is None:
        return tuple()
    elif type(x) is str or not hasattr(x, "__iter__"):
        return x,
    else:
        return x


def md5(x, N=5):
    """Return first N characters of md5 hash of argument x.
    It hashes the (default) string representation of the object.
    """
    return hashlib.md5('{}'.format(x).encode()).hexdigest()[:N]


def isnotebook():
    """Check if code is run in a Jupyter notebook.

    Returns:
        bool: True if call is made from a notebook
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


def klayout2nazca(string):
    """Convert a copy-pasted path (x,y points) from Klayout into list of points.

    Args:
        string (str): string with points to convert

    Returns:
        list of (float, float): list of points (x,y)

    Example::

        klayout2nazca("1.0, 3.5, 2.0, 5.5")

        out: [(1.0, 3.5), (2.0, 5.5)]
    """
    return [[float(x) for x in line.split('\t')]
        for line in string.strip().split('\n')]


def pointdxdy(xy, i, w):
    """Helper function to return the left and right point to the starting
    point of a line segment from (x1, y1) to (x2, y2) with width w.

    The line may contain many points, but only the first two are used.

    Args:
        xy (list of (float, float)): list of points (x,y)
        i (int): index of point in list
        w (float): width of line segment

    Returns:
        (float, float): point (dx, dy)
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
        xy (float, float): list of points (x,y)
        i (int): index of point in list
        dxy (float, float): point (dx, dy)

    Returns:
        tuple of 4 points (x,y): corner points of the line segment
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
    """Helper function to intersect two lines.

    Intersection point (xi, yi)
    of two lines that go through (x0, y0), (x1, y1) and (x2, y2), (x3, y3).
    Fails for parallel lines (divide by zero), but the algorithm in
    polyline2polygon ensures this is not the case.

    Args:
        list of four points (x,y)

    Returns:
        (float, float): intersection point (xi,yi)
    """
    x0, y0 = xy0
    x1, y1 = xy1
    x2, y2 = xy2
    x3, y3 = xy3
    D  = (x3-x2)*(y1-y0)-(x1-x0)*(y3-y2)
    Dx = (x3-x2)*(y2-y0)-(x2-x0)*(y3-y2)
    if abs(D) < 1e-16: # Don't divide, just return the point "in between".
        return ((x0+x1+x2+x3)/4, (y0+y1+y2+y3)/4)
    xi = x0 + Dx/D * (x1-x0)
    yi = y0 + Dx/D * (y1-y0)
    return (xi, yi)


def polyline_length(xy):
    """Return the lenght of the polyline, which is the sum of the line
    segments in the polyline.

    Args:
        xy (list): list of (x,y) points that hold the polygon.

    Returns:
        length (float): the length of the polyline.
    """
    length = 0
    for i in range(1, len(xy)):
        length += sqrt((xy[i][0]-xy[i-1][0])**2+(xy[i][1]-xy[i-1][1])**2)
    return length


def polyline2polygon(xy, width=2, miter=0.5):
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
        width (float): width of the polyline (default 2)
        miter (float): maximum fraction of the width before an extra point
            is added in outside corners (default 0.5)

    Returns:
        list of (float, float): the polygon
    """

    # the fraction of the width of the line segments that is used to
    # determine if a single point is sufficient to describe the outline, or
    # that two points are needed (miter limit).
    dsqrmax = (miter * width)**2
    n = len(xy)
    if n < 2:
        raise ValueError("Polyline2polygon: need at least 2 points for polyline.")
    # Start with the first two points
    dxy1 = pointdxdy(xy, 0, width)
    cxy1 = corner(xy, 0, dxy1)
    xy_start = [cxy1[0]]
    xy_end = [cxy1[1]]
    for i in range(1, n-1): # loop over the points in the polyline
        dxy0 = dxy1
        # Shift corner points from next to current segment.
        cxy0 = cxy1
        # Get corner points for next segment.
        dxy1 = pointdxdy(xy, i, width)
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


def transform_polygon(points, dx=0.0, dy=0.0, da=0.0, scale=1.0,
        flipx=False, flipy=False, x=0.0, y=0.0):
    """Transform a polygon by translation, rotation, scaling and/or flipping.

    The transformation first applies (dx, dy) to reposition the origin.
    Subsequently, the scale, rotate and flips are applied, where order does not matter.
    Finally, a (x, y) translation is performed.

    Args:
        polygon (list of (float, float)): points (x, y)
        dx (float): x translation in um (default = 0.0)
        dy (float): y translation in um (default = 0.0)
        da (float): a translation in deg (default = 0.0)
        scale (float): scaling factor (default = 1.0)
        flipx (bool): flip x coordinate x -> -x (default = False)
        flipy (bool): flip y coordinate y -> -y (default = False)
        x (float): final x translation (after other transformations)
        y (float): final y translation (after other transformations)

    Returns:
        (list of (float, float)): transformed polygon points
    """
    fu,fv = 1, 1
    if flipx:
        fu = -1
    if flipy:
        fv = -1
    a = radians(da)
    xy = []
    for u, v in points:
        u = (u+dx)*fu
        v = (v+dy)*fv
        xy.append( ( x + scale*(cos(a)*u - sin(a)*v),
                     y + scale*(sin(a)*u + cos(a)*v) ))
    return xy

