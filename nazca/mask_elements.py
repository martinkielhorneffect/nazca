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

# @author: Ronald Broeke (c) 2016-2017
# @email: ronald.broeke@brightphotonics.eu
#-----------------------------------------------------------------------

"""
Define mask elements such straight and bent waveguide via template functions.

This module is to define PDKs elements via closures on nested functions.
It is not intended for use in mask design directly.
"""

#import sys
#from types import FunctionType
#from collections import defaultdict
from math import sin, cos, acos, sqrt
from itertools import count
import warnings

import numpy as np

from .netlist import Cell, Pin
from . import cfg
from .gds_base import gds_db_unit as gridsize
from .util import md5
from .font import text as textput
import nazca as nd


def layeriter(xs=None, layer=None, dogrowy=False):
    """Generator yielding all layers in a xsection.

    Args:
        xs (str): xsection name
        layer (int): layer number
        growy (bool): temporary switch to also yield growy

    Yields:
        layer, growx, accuracy: iterate over all layers in <xs> and <layer>
    """
    if xs is None and layer is None:
        xs = cfg.default_xs_name
        if xs not in cfg.XSdict.keys():
            _handle_missing_xs(xs)

    if layer is not None:
        layer = nd.get_layer(layer)
        lt = cfg.layer_table
        growx, growy = 0, 0 # a layer not as part of a xs has no grow
        try:
            lineitem = lt[(lt['layer'] == layer[0]) & (lt['datatype'] == layer[1])]
            accuracy = lineitem['accuracy']
            if accuracy.empty:
                raise Exception()
            else:
                accuracy = accuracy.iloc[0]
        except:
            print("Error: accuracy not defined for layer {}.".format(layer))
            accuracy = 0.001
        if dogrowy:
            yield layer, growx, growy, accuracy
        else:
            yield layer, growx, accuracy

    if xs is not None:
        try:
            ML = cfg.XSdict[xs].mask_layers
            if ML.empty:
                warnings.warn("xsection '{0}' has no layers. "\
                    "Continuing export using fallback layer {1}.\n"\
                    "Recommended solutions: Use a different xsection "\
                    "or add layers to this xsection:\n"\
                    "add_layer2xsection('{0}', layer=<num>).".\
                    format(xs, cfg.default_dump_layer), stacklevel=0)
                if dogrowy:
                     yield (cfg.default_dump_layer, 0), 0, 0, 0.1
                else:
                    yield (cfg.default_dump_layer, 0), 0, 0.1
            else:
                layer_info = ML[['layer', 'datatype', 'growx', 'growy', 'accuracy']]
                for i, layer, datatype, growx, growy, accuracy in layer_info.itertuples():
                    if dogrowy:
                        yield (layer, datatype), growx, growy, accuracy
                    else:
                        yield (layer, datatype), growx, accuracy
        except:
            _handle_missing_xs(xs)
            if dogrowy:
                yield (cfg.default_dump_layer, 0), 0, 0, 0.1
            else:
                yield (cfg.default_dump_layer, 0), 0, 0.1


def _handle_missing_xs(xs):
    """Handle missing xsection.

    If xsection <xs> is not defined, create a default xsection with layer.
    A warning is issued, unless the Nazca default xsection is created.

    Returns:
        None
    """
    if xs not in cfg.XSdict.keys():
        if xs != cfg.default_xs_name:
            print("Warning: No xsection named '{0}'. "\
                "Correct the name or add the xsection with:\n"\
                "add_xsection('{0}') to get rid of this warning. "\
                "Already available xsections are {1}.\n".\
                format(xs, list(cfg.XSdict.keys())))
        nd.add_xsection(xs)
        nd.add_layer2xsection(xs,layer=cfg.default_dump_layer)
    return None

#==============================================================================
#   Waveguide element definitions
#==============================================================================
cnt = 0 #ordinal counter for unique naming
def Tp_straight(length=10, width=1.0, xs=None, layer=None, edge1=None,
        edge2=None, name=None):
    """Template for creating parametrized straight waveguide function.

    Args:
        length (float): length of waveguide
        width (float): width of waveguide
        xs (str): xsection of taper
        layer (int | str): layer number or layername
        edge1 (function): optional function F(t) describing edge1 of the waveguide
        edge2 (function): optional function G(t) describing edge2 of the waveguide

    Returns:
        function: Function returning a Cell object with a straight guide
    """

    def cell(length=length, width=width, xs=xs, layer=layer, edge1=edge1,
             edge2=edge2, name=name):
        """Create a straight waveguide element.

        Args:
            length (float): length of waveguide
            width (float): width of waveguide
            xs (str): xsection of waveguide
            layer (int | str): layer number or layername
            edge1 (function): optional function F(t) describing edge1 of the waveguide
            edge2 (function): optional function G(t) describing edge2 of the waveguide

        Returns:
            Cell: straight element
        """

        if name is None:
            name = 'straight'

        with Cell(name=name, cnt=True) as guide:
            guide.instantiate = False
            guide.put_pin(name='a0', width=width, xs=xs, show=True, connect=(0, 0, 180))
            guide.put_pin(name='b0', width=width, xs=xs, show=True, connect=(length, 0, 0))
            guide.length_geo = length

            for lay, growx, growy, acc in layeriter(xs, layer, dogrowy=True):
                if edge1 is None:
                    outline = [(0-growy, 0.5*width+growx), (length+growy, 0.5*width+growx),
                        (length+growy, -0.5*width-growx), (0-growy, -0.5*width-growx)]
                else:
                    if edge2 is None:
                        edge2 = edge1
                    Fp1 = []
                    Fp2 = []
                    for t in np.linspace(0, 1, 20):
                        Fp1.append((length*t, edge1(t)+growx))
                        Fp2.append((length*t, -edge2(t)-growx))
                    outline = Fp1 + list(reversed(Fp2))

                guide.put_polygon(layer=lay, points=outline)
        return guide
    return cell


def Tp_taper(length=100, width1=2, width2=3, xs=None, layer=None, name=None):
    """Template for creating a parametrized linear taper function.

    Args:
        length (float): length of the taper
        width1 (float): width at start
        width2 (float): width at end
        xs (str): xsection of taper
        layer (int | str): layer number or layername

    Returns:
        function: Function returning a Cell object with a linear taper
    """

    def cell(length=length, width1=width1, width2=width2, xs=xs, layer=layer,
             name=name):
        """
        Create a taper element.

        Args:
            length (float): length of taper
            width1 (float): start width of taper
            width2 (float): end width of taper
            xs (str): xsection of taper
            layer (int | str): layer number or layername

        Returns:
            Cell: taper element
        """

        #nonlocal name
        if name is None:
            name = 'taper'
        # Cellname with md5 hash
        #h = md5('{}{}{}{}'.format(length, width1, width2, xs))
        #fullname = name="{}_{}".format(name, h)
        #if fullname in cfg.cellnames.keys():
        #    return cfg.cellnames[fullname]

        with Cell(name=name, cnt=True) as taper:
            taper.instantiate = False
            if abs(length) < gridsize / 2:
                return taper # Empty taper
            if width1 > width2:
                pin = ['b0', 'a0']
                width1, width2 = width2, width1
            else:
                pin = ['a0', 'b0']

            taper.put_pin(name=pin[0], width=width1, xs=xs, show=True, connect=(0, 0, 180))
            taper.put_pin(name=pin[1], width=width2, xs=xs, show=True, connect=(length, 0, 0))
            taper.length_geo = length


            for lay, growx, growy, acc in layeriter(xs, layer, dogrowy=True):
                 # Set widths for this layer
                twidth1 = width1 + 2 * growx
                twidth2 = width2 + 2 * growx
                outline = [(0-growy, twidth1 / 2), (length+growy, twidth2 / 2),
                    (length+growy, -twidth2 / 2), (0-growy, -twidth1 / 2)]
                taper.put_polygon(layer=lay, points=outline)
        return taper
    return cell


def Tp_ptaper(length=100, width1=1.0, width2=3.0, xs=None, layer=None,
              name=None):
    """Template for creating a parametrized parabolic taper  function.

    Args:
        length (float): length of the taper
        width1 (float): width at start
        width2 (float): width at end
        xs (str): xsection of taper
        layer (int | str): layer number or layername

    Returns:
        function: Function returning a Cell object with a ptaper
    """


    def cell(length=length, width1=width1, width2=width2, xs=xs, layer=layer,
             name=name):
        """Create a parabolic taper element.

        Args:
            length (float): length of taper
            width1 (float): start width of taper
            width2 (float): end width of taper
            xs (str): xsection of taper
            layer (int | str): layer number or layername

        Returns:
            Cell: parabolic taper element
        """
        if name is None:
            name = 'ptaper'

        with Cell(name=name, cnt=True) as taper:
            taper.instantiate = False
            if abs(length) < gridsize / 2:
                return taper # Empty taper
            if width1 > width2:
                width1, width2 = width2, width1
                direction = -1
                x0 = length
            else:
                direction = 1
                x0 = 0

            taper.put_pin(name='a0', width=width1, xs=xs, show=True, connect=(0, 0, 180))
            taper.put_pin(name='b0', width=width2, xs=xs, show=True, connect=(length, 0, 0))
            taper.length_geo = length

            for lay, grow, acc in layeriter(xs, layer):
                # Set widths for this layer
                twidth1 = width1 + 2*grow
                twidth2 = width2 + 2*grow
                if abs((twidth1 - twidth2) / (4*length)) < gridsize:
                    points=[
                        (0, 0.5*twidth1), (length, 0.5*twidth2),
                        (length, -0.5*twidth2), (0, -0.5*twidth1)
                    ]
                    taper.put_polygon(layer=lay, points=points)
                    continue
                # y = a * x**2
                a = 4 * length / (twidth2**2 - twidth1**2)
                y2 = y0 = a * (0.5*twidth1)**2
                w_tap1 = twidth1
                ptop = [(x0,  0.5*twidth1)]
                pbot = [(x0, -0.5*twidth1)]
                while y2 - y0 < length:
                    w_tap2 = w_tap1 + 4*acc + 4*sqrt(acc*(w_tap1 + acc))
                    y2 = a * (0.5*w_tap2)**2
                    if y2 - y0 < length:
                        ptop.append((x0+direction*(y2-y0),  0.5*w_tap2))
                        pbot.append((x0+direction*(y2-y0), -0.5*w_tap2))
                        w_tap1 = w_tap2
                    else:
                        break
                ptop.append((x0+direction*length,  0.5*twidth2))
                pbot.append((x0+direction*length, -0.5*twidth2))
                taper.put_polygon(layer=lay, points=ptop + list(reversed(pbot)))
        return taper
    return cell


def __get_offset(xs, width, radius):
    """Get the offset function."""
    offset = 0
    if xs is None:
        return offset
    try:
        F = cfg.XSdict[xs].os
        try:
            offset = float(F)
        except:
            offset = F(width=width, radius=radius)
    except:
        _handle_missing_xs(xs)
        cfg.XSdict[xs].os = 0
    return offset


def Tp_arc(radius=10, width=1.0, angle=90, xs=None, layer=None, offset=None,
           name=None):
    """Template for creating a parametrized circular arc waveguide function.

    Args:
        radius (float): radius at the center line of the arc in um.
        width (float): width of the arc in um.
        angle (float): angle of arc in degree (default = 90).
        xs (str): xsection of taper
        layer (int | str): layer number or layername
        offset (float | function): positive offset reduces radius.
            The offset can be a function F(width, radius) that returns a float

    Returns:
        function: Function returning a Cell object with an arc
    """
    def cell(radius=radius, width=width, angle=angle, xs=xs, layer=layer,
            offset=offset, name=name):
        """Create a circular arc element.

        A straight-bend offset is included when it has been defined in the
        xsection used.

        Args:
            radius (float): radius at the center line of the arc in um.
            width (float): width of the arc in um.
            angle (float): angle of arc in degree (default = 90).
            xs (str): xsection of taper
            layer (int | str): layer number or layername
            offset (float | function): positive offset reduces radius.
                The offset can be a function F(width, radius) that returns a float

        Returns:
            Cell: circular arc element
        """

        if name is None:
            name = 'arc'
        ang = np.radians(angle)
        sign = np.sign(ang)
        radius = abs(radius)

        if offset is None:
            offset = __get_offset(xs, width, radius)

        Nmax = 98 # 2x98 + 4 = 200
        with Cell(name=name, cnt=True) as guide:
            guide.instantiate = False
            guide.put_pin(name='a0', width=width, xs=xs, show=True, connect=(0, 0, 180))
            guide.put_pin(name='b0', width=width, xs=xs, show=True,
                connect=(radius*sin(abs(ang)), sign*radius*(1-cos(ang)), angle))
            guide.length_geo = abs(radius*ang)

            if ang == 0:
                return guide

            for lay, grow, acc in layeriter(xs, layer):
                # Start and end sections are different from other sections
                Ro = radius + width * 0.5 - offset + grow # outer radius
                Ri = radius - width * 0.5 - offset - grow # inner radius
                if Ri < 0:
#                    raise Exception('Inner arc radius too small: inner_radius={:.3f}, offset={:.3f}'.\
#                        format(radius-0.5*width, offset))
                    print("Warning: Inner arc radius too small: inner_radius={:.3f}, offset={:.3f}, xs='{}'".\
                        format(radius-0.5*width, offset, xs))
                    Ri = 0
                da = 2 * acos(Ro / (Ro + acc)) # step angle for accuracy
                N = int((abs(ang) / da) + 1) # +1 for rounding
                da = ang / N # step angle for N
                u = np.linspace(da / 2, ang - da / 2, N)
                Roe = Ro * (da / 2) / sin(da / 2) # effective outer radius
                Rie = Ri * (da / 2) / sin(da / 2) # effective inner radius

                p1 = [(Roe * sin(abs(a)), sign * (radius - Roe * cos(a)))
                    for a in u]
                p2 = [(Rie * sin(np.abs(a)), sign * (radius - Rie * cos(a)))
                    for a in u]

                pstart = [(0, sign * (radius - Ri)), (0, sign * (radius - Ro))]
                pend = [
                    (Ro * sin(abs(ang)),
                    sign * (radius - Ro * cos(ang))),
                    (Ri * sin(abs(ang)),
                    sign * (radius - Ri * cos(ang)))
                ]

                section = list(range(0, len(p1), Nmax))
                for s in section[:-1]:
                    outline = pstart\
                        + p1[s:s + Nmax]\
                        + list(reversed(p2[s:s + Nmax]))
                    guide.put_polygon(layer=lay, points=outline)
                    pstart = [p2[s + Nmax - 1], p1[s + Nmax - 1]]

                outline = pstart\
                     + p1[section[-1]:section[-1] + Nmax]\
                     + pend\
                     + list(reversed(p2[section[-1]:section[-1] + Nmax]))
                guide.put_polygon(layer=lay, points=outline)
        return guide
    return cell


#TODO: not supposed to be defined in this template module.
wherecnt = count()
def whereami(text='here', size=100, pin=None):
    """Show current pointer position as arrow in the layout.

    Args:
        text (str): annotation text
        size (float): size of the annotation

    Returns:
        None
    """
    layer = 500
    nd.cp.push()
    if pin is None:
        pin = nd.cp.here()
    points = [
        (0, 0), (-0.5, 0.5), (-0.4, 0.25), (-1, 0.3), (-1, -0.3),
        (-0.4, -0.25), (-0.5, -0.5)]
    points = [(x*size,y*size) for x, y in points]
    with Cell('I am'.format(wherecnt)) as crisis:
        nd.Polygon(layer=layer, points=points).put(0)
        nd.text(text+' ', height=size/4, align='rc', layer=layer).put(0)

    crisis.put(pin)
    nd.cp.pop()


#==============================================================================
# create elements
#==============================================================================
strt   = Tp_straight()
bend   = Tp_arc()
ptaper = Tp_ptaper()
taper  = Tp_taper()

#These will be removed:
#sw = Tp_straight()
#cw = Tp_arc()
#tw = Tp_ptaper()
#tr = Tp_taper()

#==============================================================================
# #create a default cell and set cp.
#==============================================================================
cfg.defaultcell = Cell(name=cfg.defaultcellname)
cfg.cp = cfg.defaultcell.org


