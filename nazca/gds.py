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
# Higher level GDS routines.
#
# (c) 2011 Xaveer Leijtens
# (c) 2017 Ronald Broeke
"""
Module with higher level (cell level) routines to open a GDS stream
and add cells and instances of cells to it.

Example:
    Create a simple gds stream with two empty cells with cell_1 under cell_2::

        import nazca.gds as gds

        strm = gds.layout_open()

        strm += gds.cell_open(name='cell_1')
        strm += gds.cell_close()

        strm += gds.cell_open(name='cell_2')
        strm += gds.cell_reference(xy=(0, 0), name='cell_1')
        strm += gds.cell_close()

        strm += gds.layout_close()

        open(simple_stream.gds, 'bw').write(strm)
"""

from . import gds_base as gbase
from ._version import __version__

# Bit states for the STRANS Bit Array.
GDS_BIT_REFLECT  = 0x8000
GDS_BIT_ABSMAGN  = 0x0004
GDS_BIT_ABSANGLE = 0x0002

__all__ = ["layout_open", "cell_open", "layout_close", "cell_close",
           "cell_reference"]


def layout_open(name='Nazca {}'.format(__version__)):
    """Open a GDS stream.

    Args:
        name (str): name of the layout

    Returns:
        bytestring: stream for opening a layout
    """
    strm = gbase.gds_header() \
         + gbase.gds_bgnlib() \
         + gbase.gds_libname(name) \
         + gbase.gds_units()
    return strm


def cell_open(name):
    """Open a GDS cell.

    Args:
        name (str): cell name

    Returns:
        bytestring: stream for cell opening
    """
    return gbase.gds_bgnstr() + gbase.gds_strname(name)


def cell_close():
    """Close a GDS cell.

     Returns:
        bytestring: stream for cell closure
    """
    return gbase.gds_endstr()


def cell_reference(xy, name, adeg=0, mag=1, flip=False, array=None): # angle in degrees
    """Add a cell reference (instance) to a GDS cell.

    Args:
        xy (): (x, y) position of the instance
        name (str): cell name to instantiate
        adeg (float); cell angle in degrees (default = 0)
        mag (float); scaling factor (default = 1.0)
        flip (bool): flag is cell is mirrored (default = False)
        array (dx, dy): create an array reference if values are given

    Returns:
        bytestring: stream for a cell reference
    """
    if array is None:
        strm = gbase.gds_sref() + gbase.gds_sname(name)
        xy = [xy]
    else:
        strm = gbase.gds_aref() + gbase.gds_sname(name)
        col, row = array[0], array[2]
        xy = [xy] + \
             [[xy[0]+array[1][0], xy[1]+array[1][1]]] + \
             [[xy[0]+array[3][0], xy[1]+array[3][1]]]

    # Insert STRANS record if needed. If this record is omitted, the
    # defaults for the element are no reflection, non-absolute
    # magnification, and non-absolute angle (which is what we want anyway).
    # So insert only if needed for flipping.
    if flip:
        strm += gbase.gds_strans(GDS_BIT_REFLECT)
    if adeg != 0:
        strm += gbase.gds_angle(adeg)
    if mag != 1:
        strm += gbase.gds_mag(mag)
    if array is not None:
        strm += gbase.gds_colrow(col, row)
    strm += gbase.gds_xy(xy, False) + gbase.gds_endel()
    return strm

def layout_close():
    """Close a GDS stream.

    Returns:
        bytestring: stream for closing a layout
    """
    return gbase.gds_endlib()
