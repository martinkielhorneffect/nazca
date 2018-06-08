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
#==============================================================================
# (c) 2016-2017 Ronald Broeke, Katarzyna Lawniczuk
#==============================================================================

"""Module defining black box templates for PDK implementation."""



from math import atan, tan, sin, cos, degrees, radians

import nazca as nd
import nazca.pdk_template_core as pdk



def Tp_PhotoDetector(
        length=100, width=50, buf=20,
        name='BBname', groupname='',
        pinwidth=None, xs=None,
        icon=None):
    """Template for a photodetector with option to draw a metal DC pad.

    Returns:
        function that generates a Cell object
    """
    @pdk.hashme(name, 'length')
    def cell(length=length):
        """Create a PhotoDetector cell.

        Args:
            length (float): length of the diode in um
            pad (bool): flag to add a bond pad

        Returns:
            Cell: photo diode element
        """
        with nd.Cell(hashme=True) as C:
            C.groupname = groupname
            C.default_pins('a0', 'a0')
            C.foundry_spt = []
            bb_length = length+buf
            pdk.addBBmap(name, params=(length))
            nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).put(0, 0, 180)
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0']).put(bb_length)

            pdk.put_stub(['a0', 'c0'])
            pdk.put_boundingbox('org', bb_length, width)
            if icon:
                icon(bb_length, width).put(0)
        return C
    return cell


def Tp_PhotoDetectorRF(
        length=100, width=50,
        name='BBname', groupname='',
        pinwidth=None, xs=None,
        spaceGS=10,
        icon=None):
    """Template for a photodetector with an option to draw a metal RF GSG pad.

    Returns:
        function that generates a Cell object
    """
    @pdk.hashme(name)
    def cell(length=length):
        """Create a photodetector RF cell.

        Returns:
            Cell
        """
        cshift = 0
        with nd.Cell(hashme=True) as C:
            C.default_pins('a0', 'a0')
            C.groupname = groupname
            C.foundry_spt = []
            pdk.addBBmap(name)
            nd.Pin(name='a0', xs=xs['a0'], width=pinwidth['a0']).\
                put(0, 0, 180)
            nd.Pin(name='c0', xs=xs['c0'], width=pinwidth['c0']).\
                put(length, cshift-spaceGS-pinwidth['c0'])
            nd.Pin(name='c1', xs=xs['c1'], width=pinwidth['c1']).\
                put(length, cshift)
            nd.Pin(name='c2', xs=xs['c2'], width=pinwidth['c2']).\
                put(length, cshift+spaceGS+pinwidth['c2'])

            pdk.put_stub(['a0'])
            pdk.put_boundingbox('org', length, width)
            if icon:
                icon(length, width).put(0)
            pdk.put_stub(['c0', 'c1', 'c2'])
        return C
    return cell


