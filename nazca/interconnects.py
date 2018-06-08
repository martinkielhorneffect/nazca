#!#usr/bin/env python3
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
# 2017 (c) Ronald Broeke
#-----------------------------------------------------------------------
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 11:24:14 2017

@author: rgb
Interconnect module.
"""

import sys
import numpy as np
from numpy import sign
import math as m
import nazca as nd
import nazca.cfg as cfg
import nazca.cp as cp
import nazca.geometries as geom
import nazca.pdk_template_core as tdk



def PosRad(a):
    """Clip angle to [0, 360>"""
    if a < 0:
        pa = a-2*m.pi*(int(a/(2*m.pi)-1))
    else:
        pa = a-2*m.pi*(int(a/(2*m.pi)))
    if pa == 2*m.pi:
        pa = 0
    #printf("pa=", degree(pa), ", a=", a, ", int(a/rad(360))=", int(a/rad(360)), "\n")
    return pa


def NegRad(a):
    """Clip angle to <-360, 0]"""
    pn = PosRad(a)-2*m.pi
    if pn <= -2*m.pi:
        pn = 0
    #printf("pn=", degree(pn), "\n")
    return pn


class Interconnect():
    """Interconnect class for drawing waveguides.

    An Interconnect object can be configured to match a specifc foundry.
    This includes properties like width, xsection, straight-bend offset and others.

    Example:
        Create two Interconnect objects, each for a different kind of waveguide::

            ic1 = Interconnect(width=2.0, radius=20)
            ic2 = Interconnect(width=1.0, radius=50)
    """

    def __init__(self, radius=None, width=None, angle=90, xs='nazca',
            layer=None, adapt_width=False, adapt_xs=False):
        """Contruct an Interconnect object.

        If a xsection is provided in <xs> then values for <radius> and <width>
        will be copied from <xs>.
        If <radius> and/or <width> are explicitly set in __init__ then they will
        take priority over the values in <xs>.
        If <xs> nor <radius> and/or <width> are set they default to
        values in the cfg module

        Args:
            radius (float): default radius in um
            width (float): default waveguide width im um
            angle (float): default angle of a bend (default = 90 degrees)
            xs (str): waveguide xsection (default = 'nazca')
            layer (str): layer to draw interconnect in
            adapt_width (bool): adapt interconnect width to the pin in connects to
                (default = False)
            adapt_xs (bool): adapt interconnect width to the pin in connects to
                 (default = False)

        Returns:
            None
        """

        if radius is None:
            try:
                self.radius = nd.get_xsection(xs).radius
            except:
                self.radius = cfg.default_xs_radius
        else:
            self.radius = radius

        if width is None:
            try:
                self.width = nd.get_xsection(xs).width
            except:
                self.width = cfg.default_xs_width
        else:
            self.width = width

        self.angle = angle #angle
        self.length = 10 #length
        self.xs = xs
        self.layer = layer
        self.instantiate=False
        self.line = nd.Tp_straight(xs=self.xs, layer=self.layer)
        self.Farc = nd.Tp_arc(xs=self.xs, layer=self.layer)
        self.__ptaper = nd.Tp_ptaper(xs=self.xs, layer=self.layer)
        self.__taper = nd.Tp_taper(xs=self.xs, layer=self.layer)
        self.adapt_width = adapt_width
        self.adapt_xs = adapt_xs
        self.arrow = tdk.make_pincell()


    def arc(self, radius=10, width=1.0, angle=90, xs=None,
            layer=None, offset=None, name=None):
        # put arc in a function to make it overloadable in derived classed.
        # (a derived function can't overload a base attribute)
        return self.Farc(radius=radius, width=width, angle=angle, xs=xs,
            layer=layer, offset=offset, name=name)


    def _getpinin(self, pin):
        """Return the pin as is or the default_in pin if an instance is provided."""
        if isinstance(pin, nd.Instance):
            pin = pin.pin[pin.cnode.cell.default_in]
        return pin

    def _getpinout(self, pin):
        """Return the pin as is or the default_out pin if an instance is provided."""
        if isinstance(pin, nd.Instance):
            pin = pin.pin[pin.cnode.cell.default_out]
        return pin

    def _getwidth(self, pin, width, xs, end=False):
        """Return width based on interconnect rules.

        Returns:
            width
        """
        if width is not None:
            return width
        if self.adapt_width and not end:
            if pin is None:
                pin = cp.here()
            if pin.width is not None:
                width = pin.width
            #elif xs is not None and width is None:
            #    width = None
        if width is None:
            try:
                width = nd.get_section(xs).width
            except:
                width = self.width
        return width


    def _getradius(self, pin, radius, xs):
        if self.adapt_xs:
            if pin is None:
                pin = cp.here()
            if xs is None:
                try:
                    xs = pin.xs
                except:
                    xs = self.xs
        if radius is None:
            try:
                radius = nd.get_section(xs).radius
            except:
                radius = self.radius
        return radius


    def _getxs(self, pin, xs):
        """Obtain the xs based on the Interconnect settings."""
        if xs is not None:
            return xs
        if self.adapt_xs:
            if pin is None:
                pin = cp.here()
            try:
                if pin.xs is not None:
                    xs = pin.xs
            except:
                 xs = self.xs
                #raise Exception('No xsection defined in pin. Add a xs attribute.')
        if xs is None:
            xs = self.xs
        return xs


    def strt(self, length=None, width=None, pin=None, xs=None, edge1=None,
             edge2=None, name=None, arrow=True):
        """Create a straight waveguide.

        Args:
            length (float): length of guide in um
            width (float): width of guide in um
            pin (Node): optional Node for modeling info
            xs (str): optionals xsection of guide
            layer (int | str): layer number or layername
            edge1 (function): optional function F(t) describing edge1 of the waveguide
            edge2 (function): optional function G(t) describing edge2 of the waveguide

        Returns:
            Cell: waveguide element

        Example:
            Create and place a straight waveguide::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.strt(length=20)
                guide.put()
                nd.export_plt()
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width = self._getwidth(pin, width, xs)

        if length is None:
            length = self.length

        with nd.Cell('strt_'+xs, instantiate=self.instantiate, cnt=True) as ICcell:
            nd.Pin('a0', width=width, xs=xs).put(0, 0, 180)
            self.line(length=length, width=width, xs=xs, edge1=edge1,
                edge2=edge2, name=name).put(0)
            nd.Pin('b0', width=width, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'])

        if pin is not None:
            cfg.cp = pin
        return ICcell


    def bend(self, radius=None, angle=None, width=None, pin=None, xs=None,
            name=None, arrow=True):
        """Create a bent waveguide (circular arc).

        Args:
            radius (float): radius at the center line of the arc in um
            width (float): width of the arc in um
            angle (float): angle of arc in degree (default = 90)
            pin (Node): optional Node for modeling info
            xs (str): optiinal xsection of bend

        Returns:
            Cell: circularly bent waveguide element

        Example:
            Create and place a bend::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.bend(angle=45)
                guide.put()
                nd.export_plt()
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width = self._getwidth(pin, width, xs)
        radius = self._getradius(pin, radius, xs)

        if angle is None:
            angle = self.angle

        with nd.Cell('bend_'+xs, instantiate=self.instantiate, cnt=True) as ICcell:
            nd.Pin('a0', width=width, xs=xs).put(0, 0, 180)
            self.arc(radius=radius, width=width, angle=angle, xs=xs, name=name).put(0)
            nd.Pin('b0', width=width, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'])

        if pin is not None:
            cfg.cp = pin
        return ICcell


    def ptaper(self, length=None, width1=None, width2=None, pin=None, xs=None,
               name=None, arrow=True):
        """Create a parabolic taper.

        Args:
            length (float): length of taper
            width1 (float): start width of taper
            width2 (float): end width of taper
            pin (Node): optional Node for modeling info
            xs (str): optional xsection of taper

        Returns:
            Cell: parabolic taper element

        Example:
            Create and place a parabolic taper::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.ptaper(length=10, width1=2.0, width2=5.0)
                guide.put()
                nd.export_plt()
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width1 = self._getwidth(pin, width1, xs)
        width2 = self._getwidth(pin, width2, xs, end=True)

        if length is None:
            length = self.length

        with nd.Cell('ptaper_'+xs, instantiate=self.instantiate, cnt=True) as ICcell:
            nd.Pin('a0', width=width1, xs=xs).put(0, 0, 180)
            self.__ptaper(length=length, width1=width1, width2=width2,
                name=name, xs=xs).put(0)
            nd.Pin('b0', width=width2, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'])

        if pin is not None:
            cfg.cp = pin
        return ICcell


    def taper(self, length=None, width1=None, width2=None, xs=None, pin=None,
              name=None, arrow=True):
        """Create a linear taper.

        Args:
            length (float): length of taper
            width1 (float): start width of taper
            width2 (float): end width of taper
            xs (str): optional xsection of taper
            pin (Node): optional Node for modeling info

        Returns:
            Cell: linear taper element

        Example:
            Create and place a linear taper::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.taper(length=10, width1=2.0, width2=5.0)
                guide.put()
                nd.export_plt()
        """
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width1 = self._getwidth(pin, width1, xs)
        width2 = self._getwidth(pin, width2, xs, end=True)

        if length is None:
            length = self.length

        with nd.Cell('taper_'+xs, instantiate=self.instantiate, cnt=True) as ICcell:
            nd.Pin('a0', width=width1, xs=xs).put(0, 0, 180)
            self.__taper(length=length, width1=width1, width2=width2,
                name=name, xs=xs).put(0)
            nd.Pin('b0', width=width2, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'])

        if pin is not None:
            cfg.cp = pin
        return ICcell


    def strt_p2p(self, pin1=None, pin2=None, width=None, xs=None, name=None,
            arrow=True):
        """Create point-to-point straight interconnect.

        Args:
            pin1 (Node | Instance): start of waveguide
            pin2 (Node | Instance): end of waveguide
            width (float): width of waveguide
            xs (str): optional xsection of waveguide

        Returns:
            Cell: straight waveguide element

        Example:
            Create and place a straight guide between points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.strt_p2p(pin1=(0, 0), pin2=(10, 10))
                guide.put()
                nd.export_plt()
        """
        if pin1 is None:
            pin1 = cp.here()
        if pin2 is None:
            pin1, pin2 = cp.here(), pin1

        pin1b, T = nd.parse_pin(pin1)
        pin1 = pin1b.move(*T)
        pin2b, T = nd.parse_pin(pin2, rot=True, default='out')
        pin2 = pin2b.move(*T)

        xs = self._getxs(pin1, xs)
        width = self._getwidth(pin1, width, xs)

        x, y, a = nd.diff(pin1,  pin2)
        length = m.sqrt(x**2 + y**2)
        b = m.degrees(m.atan2(y, x))
        if name is None:
            name = 'sw'
        with nd.Cell(name, instantiate=self.instantiate, cnt=True) as ICcell:
            nd.Pin('a0', width=width, xs=xs).put(0)
            self.line(length=length, width=width, xs=xs).put(ICcell.pin['a0'].rot(180+b))
            nd.Pin('b0', width=width, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'])
        cfg.cp = pin1
        return ICcell


    def rot2ref_solve(self, pin=None, ref=None, angle=0, cw=None):
        """
        Calculate angle to rotate from <pin> to a certain reference direction <ref>.

        Args:
            pin (Node): starting pin (default = cp)
            ref (Node): reference pin (default = org)
            angle (float): rotation with repect to ref in [Degrees] (default = 0)
            cw (bool): angle direction clockwise or counter clockwise (default is shortest bend)

        Returns:
            float: angle to rotate from <pin> to <ref> + <a>
        """

        nd.cp.push()

        if pin is None:
            print('Error: source pin not specified in rot2ref.')
            sys.exit()

        #xs = self._getxs(pin, xs)
        #width  = self._getwidth(pin, width, xs)

        if ref is None:
            ref = cfg.cells[-1].pin['org']
        else:
            ref = self._getpinout(ref)

        x, y, a = nd.diff(ref.rot(angle), pin)
        if a >= 180:
            a -= 360

        if cw:
            if a>0:
                a -= 360
        if not cw and cw is not None:
            if a<0:
                a += 360

        nd.cp.pop()
        return -a


    def rot2ref(self, pin=None, ref=None, angle=0, length1=0, length2=0,
            cw=None, width=None, xs=None, radius=None, name=None, arrow=True):
        """
        Rotate a waveguide from <pin> to a specific <angle> with respect to reference <ref>.

        Args:
            pin (Node): starting pin (default = cp)
            ref (Node): reference pin (default = org)
            angle (float): rotation with repect to ref (default = 0)
            length1 (float): optional straight starting section (default = 0)
            length2 (float): optional straight ending section (default = 0)
            cw (bool): angle direction clockwise or counter clockwise (default is shortest bend)
            width (float): width of waveguide
            xs (str): optional xsection of waveguide
            radius (float): radius at the center line of the bend in um

        Returns:
            Cell: waveguide element rotating from <pin> into to the desired direction

        Example:
            Rotate to angle 125 degree w.r.t. 'org' after a straight guide of 100 um::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.rot2ref(angle=125, length1=10)
                guide.put()
                nd.export_plt()
        """

        if pin is None:
            pin = nd.cp.here()
        else:
            pin = self._getpinout(pin)

        xs = self._getxs(pin, xs)
        width  = self._getwidth(pin, width, xs)
        radius  = self._getradius(pin, radius, xs)
        a = self.rot2ref_solve(pin, ref, angle, cw)

        if name is None:
            name = 'rot2ref'
        with nd.Cell(name, instantiate=self.instantiate, cnt=True) as ICcell:
            nd.Pin('a0', width=width, xs=xs).put(0, 0, 180)
            self.line(length=length1, width=width, xs=xs).put(0)
            self.arc(angle=a, radius=radius, width=width, xs=xs).put()
            self.line(length=length2, width=width, xs=xs).put()
            nd.Pin('b0', width=width, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'])

        if pin is not None:
            cfg.cp = pin
        return ICcell


    def sbend(self, radius=None, width=None, pin=None, xs=None, offset=20,
            Ltot=0, name=None, arrow=True):
        """Create an s-bend interconnect.

        Args:
            radius (float): bend radius at the center line of the arc in um
            width (float): width of the interconnect in um
            pin (Node): optional Node for modeling info
            xs (str): xsection of taper
            offset (float): lateral offset of the sbend in um
            Ltot (float): optional total forward length of the sbend in um

        Returns:
            Cell: sbend element

        Example:
            Create and place a sbend waveguide::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.sbend(offset=20)
                guide.put()
                nd.export_plt()
        """

        #TODO: check Ltot and add 'length' var
        xs = self._getxs(pin, xs)
        width  = self._getwidth(pin, width, xs)
        radius  = self._getradius(pin, radius, xs)

        L, La, Lb = 0, 0, 0
        Amax = m.radians(90)

        if 2*radius*(1-m.cos(Amax)) > abs(offset):
            A = m.acos(1-abs(offset)/(2*radius))
            if offset > 0:
                A = -A
        else:
            L = (abs(offset) - abs(2*radius*(1-m.cos(Amax))))
            A = sign(-offset)*Amax

        Lx = (2*radius+L) * m.sin(abs(A))
        dLx = abs(Ltot)-Lx

        if Ltot < 0 and dLx > 0:
            La = dLx
        elif Ltot > 0 and dLx > 0:
            Lb = dLx
        else:
            La=0
            Lb=0

        if name is None:
            name = 'cw'
        with nd.Cell('name', instantiate=self.instantiate, cnt=True) as ICcell:
            nd.Pin('a0', width=width, xs=xs).put(0, 0, 180)
            s1 = self.line(La, width, xs=xs).put(0)
            b1 = self.arc(radius, angle=-m.degrees(A), width=width, xs=xs).put()
            s2 = self.line(L, width).put()
            b2 = self.arc(radius, angle=m.degrees(A), width=width, xs=xs).put()
            s3 = self.line(Lb, width, xs=xs).put()
            nd.Pin('b0', width=width, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'])

        if pin is not None:
            cfg.cp = pin
        return ICcell


    def sbend_p2p(self, pin1=None, pin2=None, width=None, radius=None, Amax=90,
            xs=None, doStrFirst=1, Lstart=0, BendEndFlag=1, name=None,
            arrow=True):
        """Create point-to-point s-bend interconnect.

        The direction of the end pin is ignored.

        Args:
            pin1 (Node | Instance): start pin (default = cp)
            pin2 (Node | Instance): end pin
            width (float): width of the interconnect in um
            radius (float): bend radius of the interconnect in um
            xs (str): optional xsection of sbend
            Lstart (float): straight waveguide length at beginning (positive value)
                or end (negative value) of sbend

        Returns:
            Cell: sbend element

        Example:
            Create and place a sbend to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.sbend_p2p(pin1=(0), pin2=(40, 20))
                guide.put()
                nd.export_plt()
        """
        if pin1 is None:
            pin1 = cp.here()
        if pin2 is None:
            pin1, pin2 = cp.here(), pin1

        pin1b, T = nd.parse_pin(pin1)
        pin1 = pin1b.move(*T)
        pin2b, T = nd.parse_pin(pin2, rot=True, default='out')
        pin2 = pin2b.move(*T)

        xs = self._getxs(pin1, xs)
        width  = self._getwidth(pin1, width, xs)
        radius  = self._getradius(pin1, radius, xs)

        A = 0
        Ls =0
        Ltap = 0
        Amax = m.radians(Amax)

        dx, dy, da = nd.diff(pin1,  pin2)

        if Lstart < 0:
            Lstart = abs(Lstart)
            doStrFirst = 0

        # get dy from start and end position
        if Lstart < Ltap:
            Lstart = Ltap

        if 2*radius*(1-m.cos(Amax))+2*Ltap*m.sin(Amax) > abs(dy):
            #not enough offset for Amax--> no vertical straight guide
            if Ltap == 0:
                A = m.acos(1-abs(dy)/(2*radius))
            else:
                tel=0
                A = 0
                da=0.2
                damin=1e-8
                while abs(da) > damin and tel < 100:
                    if Ltap*m.sin(A)-radius*m.cos(A) < (abs(dy)-2*radius)/2.0:
                        A += da
                    else:
                        A -= da
                        da /= 2.0
                        A += da
                    tel += 1
                if A < 10*damin:
                    A=0

            A = sign(dy)*A
            Lstr = 0
        else: # use Amax angle
            A = sign(dy)*Amax
            Lstr = (abs(dy) - abs(2*radius*(1-m.cos(Amax))) -2*Ltap) # abs(m.sin(Amax))

        Lfit = dx -2*radius*abs(m.sin(A)) - (Lstr+2*Ltap)*abs(m.cos(A)) - Lstart
        La, Lb = 0, 0
        if doStrFirst == 1:
            La = Lstart-Ltap
            Lb = Lfit-Ltap
        else:
            Lb = Lstart-Ltap
            La = Lfit-Ltap

        if La < 0 or Lb < 0:
            print("(interconnect): No solution for sbend_p2p.")
            return self.strt_p2p(pin1, pin2)

        else:
            if A == 0:
                Lstr += 4*Ltap

            if Lfit >= 0:
                if name is None:
                    name = 'sbend_p2p'
                with nd.Cell(name, instantiate=self.instantiate, cnt=True) as ICcell:
                    nd.Pin('a0', width=width, xs=xs).put((0, 0, 180))
                    s1 = self.line(La, width, xs=xs).put()
                    b1 = self.arc(radius=radius, angle=m.degrees(A), width=width, xs=xs).put()
                    s2 = self.line(Lstr, width, xs=xs).put()
                    b2 = self.arc(radius=radius, angle=-m.degrees(A)*sign(dx), width=width, xs=xs).put()
                    s3 = self.line(Lb, width, xs=xs).put()
                    nd.Pin('b0', width=width, xs=xs).put()
                    if arrow:
                        self.arrow.put(ICcell.pin['a0'])
                        self.arrow.put(ICcell.pin['b0'])
                #nd.connect(pin2, ICcell)
                cfg.cp = pin1
                return ICcell

            else: # solve start- or end-point connection with an extra s-bend
                if name is None:
                    name = 'sbend'
                if doStrFirst == 0:
                    with nd.Cell("{}_{}".format(name, cnt), instantiate=self.instantiate) as ICcell:
                        nd.Pin('a0', width=width, xs=xs).put(0, 0, 180)
                        s1 = self.line(La, width, xs=xs).put()
                        b1 = self.arc(radius=radius, angle=m.degrees(A), width=width, xs=xs).put()
                        out = nd.Pin('b0', width=width, xs=xs).put()
                        if arrow:
                            self.arrow.put(ICcell.pin['a0'])
                            self.arrow.put(ICcell.pin['b0'])
             ####s2 = ml::bp_sbend_bend_p2p(in0->B1@out0, out0->end@cin : w,R, -1e-6)
                    #nd.connect(pin2, ICcell)
                    cfg.cp = pin1
                    return ICcell

                else:
                    with nd.Cell("{}_{}".format(name, cnt), instantiate=self.instantiate) as ICcell:
                        nd.Pin('a0', width=width, xs=xs).put(0, 0, 180)
                        s1 = self.line(Lb, width, xs=xs).put()
                        b1 = self.arc(radius=radius, angle=m.degrees(-A), width=width, xs=xs).put()
                        nd.Pin('b0', width=width, xs=xs).put()
                        if arrow:
                            self.arrow.put(ICcell.pin['a0'])
                            self.arrow.put(ICcell.pin['b0'])
                    #nd.connect(pin2, ICcell)
                    cfg.cp = pin1
                    return ICcell


    def __bend_strt_bend_solve(self, pin1=None, pin2=None, radius1=None,
             radius2=None, ictype='shortest'):
        """Calculate geometry for a bend_strt_bend interconnect.

        Returns:
            dict: A dictionary with bend_strt_bend solutions
        """

        pin1 = self._getpinout(pin1)
        pin2 = self._getpinin(pin2)

        Ltap = 0
        A =  pin1.move(Ltap, 0, 0) # to calculate the geometry with Ltap
        B =  pin2.move(-Ltap, 0, 180)
        dx, dy, da = nd.diff(A, B)

        # Calculate circle centers. Note that pin1 is conceptually put at (0,0,0)
        c1Lx, c1Ly = 0, radius1
        c1Rx, c1Ry = 0, -radius1
        c2Lx, c2Ly = dx-radius2*m.sin(m.radians(da)), dy+radius2*m.cos(m.radians(da))
        c2Rx, c2Ry = dx+radius2*m.sin(m.radians(da)), dy-radius2*m.cos(m.radians(da))

        solutions = {}
        if ictype in ['rr', 'rl', 'lr' 'll']:
            shapes = [ictype]
        else:
            shapes = ['rr', 'rl', 'lr', 'll']

        options = ['rr', 'rl', 'lr', 'll', 'shortest', 'all']
        if ictype not in options:
            print("Warning: ictype '{}' not defined, switching to 'shortest'.".\
                format(ictype))
            print("Valid options are {}.".format(options))
            ictype = 'shortest'

        for shape in shapes:
            if shape is 'rr':
                sx, sy = c2Rx-c1Rx, c2Ry-c1Ry
                d1, d2 = 1, -1
            if shape is 'rl':
                sx, sy = c2Lx-c1Rx, c2Ly-c1Ry
                d1, d2 = 1, 1
            if shape is 'lr':
                sx, sy = c2Rx-c1Rx, c2Ry-c1Ly
                d1, d2 = -1, -1
            if shape is 'll':
                sx, sy = c2Lx-c1Lx, c2Ly-c1Ly
                d1, d2 = -1, 1

            rr = d1*radius1 + d2*radius2
            s = m.sqrt(sx**2+sy**2)  # (sx,sy): (x,y)-distannce between circle centers
            if rr > s:
                print("WARNING (module interconnects): Radii to large.")
                found = False
                return {} #found, 0, 0, 0, 0
            else:
                gs = m.atan2(sy, sx) # angle through the circle centres at the start
                if abs(rr/s) <= 1:
                    found = True
                    gb = m.asin(rr/s) # angle through the circle centers after placing connecting straight horizontal
                else:
                    found = False
                    return {} #found, 0, 0, 0, 0

                gt = -gs+gb # angle of rotation of axis through circle centers from to put connection between circles horizontal.
                t1 = m.radians(0)+gt # angle of spoke that points to start-point bsb on the circle
                t2 = m.radians(da)+gt # angle of spoke that points to end-point bsb on the circle
                if d1 == 1:
                    b = NegRad(-t1) # angle of drawn self.arc on start self.arc
                else:
                    b = PosRad(-t1) # angle of drawn self.arc on start self.arc

                if d2 == -1:
                    e = NegRad(t2) # angle of drawn self.arc on end self.arc
                else:
                    e = PosRad(t2) # angle of drawn self.arc on start self.arc


                L = s*m.cos(gb) # length of straight self.line
                Ltot = L + radius1*abs(b)+radius2*abs(e) # total connection length

                if b == 0:
                    L += Ltap
                else:
                    L -= Ltap
                if e == 0:
                    L += Ltap
                else:
                    L -= Ltap
                if L < 0:
                    found = True

            solutions[shape] = (found, Ltot, L, b, e)

        if ictype is 'shortest':
            shortest = None
            LMin = 1e10
            for shape, geo in solutions.items() :
                found, Ltot, L, b, e = geo
                if found:
                    if Ltot < LMin:
                        LMin = Ltot
                        shortest = shape
            return {shortest: solutions[shortest]}

        elif ictype is 'all':
            return solutions

        elif ictype in ['rr', 'rl', 'lr', 'll']:
            return {ictype: solutions[ictype]}

        else:
            exit('Something wrong, no solution found.')


    def bend_strt_bend_p2p(self, pin1=None, pin2=None, radius=None, radius1=None, radius2=None,
            width=None, xs=None, ictype='shortest', name=None, arrow=True):
        """Generate a point-to-point bend-straight-bend interconnect.

        Args:
            pin1 (Node | Instance): start pin (default = cp)
            pin2 (Node | Instance): end pin
            radius1 (float): optional first bend radius in um
            radius2 (float): optional second bend radius im um
            width (float): optional waveguide width in um
            xs (str): optional xsection
            ictype (str): interconnection type (default = 'shortest')
                options: 'shortest', 'll', 'lr', 'rl', rr', 'all'

        Returns:
            Cell: bend_strt_bend element

        Example:
            Create and place a bend-straight-bend guide to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.bend_strt_bend_p2p(pin1=(0, 0, 0), pin2=(40, 20, 90))
                guide.put()
                nd.export_plt()
        """
        if pin1 is None:
            pin1 = cp.here()
        if pin2 is None:
            pin1, pin2 = cp.here(), pin1

        pin1b, T = nd.parse_pin(pin1)
        pin1 = pin1b.move(*T)
        pin2b, T = nd.parse_pin(pin2, rot=True, default='out')
        pin2 = pin2b.move(*T)

        xs = self._getxs(pin1, xs)
        width  = self._getwidth(pin1, width, xs)
        if radius is not None:
            if radius1 is None:
                radius1 = radius
            if radius2 is None:
                radius2 = radius
        radius1 = self._getradius(pin1, radius1, xs)
        radius2 = self._getradius(pin1, radius2, xs)

        curves = self.__bend_strt_bend_solve(pin1, pin2, radius1, radius2, ictype=ictype)

        instantiate = self.instantiate
        if ictype is 'all':
            instantiate = True

        cells = []
        for shape, geo in curves.items():
            noSolution, Ltot, L, b, e = geo
            if name is None:
                name = 'bend_strt_bend'
            with nd.Cell("{}_{}".format(name, shape),
                    instantiate=instantiate, cnt=True) as ICcell:
                nd.Pin('a0', width=width, xs=xs).put(0, 0, 180)
                b1 = self.arc(radius1, angle=m.degrees(b), width=width, xs=xs).put()
                s1 = self.line(L, width, xs=xs).put()
                b2 = self.arc(radius2, angle=m.degrees(e), width=width, xs=xs).put()
                nd.Pin('b0', width=width, xs=xs).put()
                if arrow:
                    self.arrow.put(ICcell.pin['a0'])
                    self.arrow.put(ICcell.pin['b0'])
            #nd.connect(pin2, ICcell)
            cfg.cp = pin1
            cells.append(ICcell)

        if not cells:
            print("(interconnect): No solution for bend_strt_bend_p2p.")
            return self.strt_p2p(pin1, pin2)
        elif len(cells) == 1:
            return cells[0]
        else:
            with nd.Cell(name, cnt=True) as ICgroup:
                p = nd.Pin('a0', width=width, xs=xs).put(0)
                for cell in cells:
                    cell.put(p.rot(180))
            #nd.connect(pin2, ICcell)
            cfg.cp = pin1
            return ICgroup


    def bend_strt_bend(self, pin=None, radius=None, radius1=None, radius2=None,
        width=None, xs=None, ictype='shortest', name=None, arrow=True):
        """Generate a bend-straight-bend connection starting at the current pointer.

        This is the same connection as 'bend_strt_bend_p2p' with pin1 = cp.
        """
        return self.bend_strt_bend_p2p(
            pin1=cp.here(), pin2=pin,
            radius=radius, radius1=radius1, radius2=radius2,
            width=width, xs=xs, ictype=ictype, name=name, arrow=arrow)


    def __strt_bend_strt_p2p_solve(self, pin1, pin2, radius):
        """Solve geometry for a strt_bend_strt interconnection.

        Returns:
            dict: A dictionary with solutions
        """
        pin1 = self._getpinout(pin1)
        pin2 = self._getpinin(pin2)

        dx, dy, da = nd.diff(pin1, pin2.rot(180))
        if da >= 180:
            da-= 360
        g = (180-da)/2.0

        message = ''
        L1 = 0
        L2 = 0

        if abs(180-da) > 1e-5 and abs(180+da) > 1e-5:
            x0 = dy/m.tan(m.radians(180-da)) + dx
            solution = True

            dx1 = abs(radius/m.tan(m.radians(g)))

            if sign(dy)*da < 0 or sign(dy)*da >= 180:
                solution = False
                message = "strt_bend_strt: Wrong direction end-point. Switching to bend_strt_bend."
            elif  abs(g-90) < 1e-4:
                solution = False
                message = "strt_bend_strt: pointers in-line. Switching to bend_strt_bend."
            else:
                Ltap = 0
                s = m.sqrt((x0-dx)**2 + dy**2)
                L1 = x0 - dx1
                L2 = s - dx1
                L1 -= Ltap
                L2 -= Ltap
                if L1>0 and L2>0:
                    solution = True
                else: solution = False

        else:
            solution = False
            message = "strt_bent_strt: angle not possible. Switching to bend_strt_bend."

        return (solution, L1, L2, da, message)


    def strt_bend_strt_p2p(self, pin1=None, pin2=None, radius=None, width=None,
            xs=None, name=None, arrow=True):
        """Create point-to-point straight-bend-straight interconnect.

        Args:
            pin1 (Node | Instance): start pin (default = cp)
            pin2 (Node | Instance): end pin
            radius (float): optional bend radius in um
            width (float): optional waveguide width in um
            xs (str): optional xsection

        Returns:
            Cell: strt_bend_strt element

        Example:
            Create and place a straight-bend-straight guide to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.strt_bend_strt_p2p(pin1=(0, 0, 0), pin2=(40, 20, 90))
                guide.put()
                nd.export_plt()
        """
        # For the calculation always rotate coordinate system to point start
        # coordinat in positive x-axis.
        if pin1 is None:
            pin1 = cp.here()
        if pin2 is None:
            pin1, pin2 = cp.here(), pin1

        pin1b, T = nd.parse_pin(pin1)
        pin1 = pin1b.move(*T)
        pin2b, T = nd.parse_pin(pin2, rot=True, default='out')
        pin2 = pin2b.move(*T)

        xs = self._getxs(pin1, xs)
        width  = self._getwidth(pin1, width, xs)
        radius  = self._getradius(pin1, radius, xs)

        shape = self.__strt_bend_strt_p2p_solve(pin1, pin2, radius)
        solution, L1, L2, da, message = shape

        if L1 > 0 and L2 > 0:
            if name is None:
                name = 'scs'
            with nd.Cell(name, instantiate=self.instantiate, cnt=True) as ICcell:
                nd.Pin('a0', width=width, xs=xs).put(0, 0, 180)
                s1 = self.line(L1, width, xs=xs).put()
                b1 = self.arc(radius, angle=da, width=width, xs=xs).put()
                s2 = self.line(L2, width, xs=xs).put()
                nd.Pin('b0', width=width, xs=xs).put()
                if arrow:
                    self.arrow.put(ICcell.pin['a0'])
                    self.arrow.put(ICcell.pin['b0'])
            #nd.connect(pin2, ICcell)
            cfg.cp = pin1
            return ICcell

        else:
            message = "strt_bend_strt: Negative straight guide. Switching to bend_strt_bend."

        # goto bend_strt_bend:
        print(message)
        return self.bend_strt_bend_p2p(pin1, pin2, radius1=radius,
            radius2=radius, width=width, xs=xs, name=name)


    def ubend_p2p(self, pin1=None, pin2=None, radius=None, width=None, xs=None,
                  length=0, name=None, arrow=True):
        """Create point-to-point u-bend interconnect.

        Args:
            pin1 (Node | Instance): start pin (default = cp)
            pin2 (Node | Instance): end pin
            radius (float): optional bend radius in um
            width (float): optional waveguide width in um
            xs (str): optional xsection of ubend
            length (float): extra straight section for longer ubend (default = 0)

        Returns:
            Cell: ubend element

        Example:
            Create and place a ubend to connect two specific points::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.ubend_p2p(pin1=(0, 0, 0), pin2=(10, 20, 90), length=10)
                guide.put()
                nd.export_plt()
        """
        # draw a 180 degree u-shape connection
        # orientation of the output pin does not matter.
        if pin1 is None:
            pin1 = cp.here()
        if pin2 is None:
            pin1, pin2 = cp.here(), pin1

        pin1b, T = nd.parse_pin(pin1)
        pin1 = pin1b.move(*T)
        pin2b, T = nd.parse_pin(pin2, rot=True, default='out')
        pin2 = pin2b.move(*T)

        xs = self._getxs(pin1, xs)
        width  = self._getwidth(pin1, width, xs)
        radius  = self._getradius(pin1, radius, xs)

        dx, dy, da = nd.diff(pin1, pin2)
        if dx < 0:
            L2 = length-dx
            L1 = length
        else:
            L1 = length+dx
            L2 = length
        if name is None:
            name = 'scs'
        with nd.Cell(name, instantiate=self.instantiate, cnt=True) as ICcell:
            nd.Pin('a0', width=width, xs=xs).put((0, 0, 180))
            s1 = self.line(L1, width, xs=xs).put(ICcell.pin['a0'].rot(180))
            s2 = self.line(L2, width, xs=xs).put(ICcell.pin['a0'].move(-dx, -dy, 180))
            b1 = self.bend_strt_bend_p2p(s1.pin['b0'].rot(0), s2.pin['b0'].rot(0),
                radius1=radius, radius2=radius, width=width, xs=xs, arrow=False).put()
            nd.Pin('b0', width=width, xs=xs).put()
            if arrow:
                self.arrow.put(ICcell.pin['a0'])
                self.arrow.put(ICcell.pin['b0'])
        cfg.cp = pin1
        return ICcell


    def mamba(self, points, radius=None, width=None, pin=None, xs=None,
              N=1, pitch=10, offset=0, polyline=True, showpins=False,
              name=None, arrow=True):
        """Creates a snake-like interconnect along a list of points (x, y).

        Start and end of a mamba are pins 'a0' and 'b0', respectively.
        To connect a mamba in absolute cell coordinates use: mamba(...).put('org', 0).

        Args:
            points: list of (x, y) positions to guide the mamba
            radius (float): optional waveguide radius (default self.radius)
            width (float): optional waveguide width (default self.width)
            pin (Node): optional Node for modeling info
            xs (str): optional xsection of mamba
            N (int): number of parallel guides in the mamba
            pitch (float): pitch of the guide if N>1
            offset (float): lateral offset in the position of all guides
            polyline (bool): boolean determining if the mamba is also drawn as polyline
                (default = True)
            showpins (bool): show the points as dots in the layout (default = False)

        Returns:
            Cell: mamba element through <points>

        Example:
            Create a mamba and attach the first point to the current pin::

                import nazca as nd
                from nazca.interconnects import Interconnect
                ic = Interconnect(width=2.0, radius=10.0)

                guide = ic.mamba(points=[(10, 10), (20, 20), (20, 50), (10, 50)])
                guide.put(0) # put first mamba point 'a0' on a pin
                guide.put('org', 0) # put mamba 'org' in 0 for absolute coordinates
                nd.export_plt()

            Hence to put a mamba at absolute coordinates of <points> in the cell::

                guide.put('org', 0)
        """
        nd.cp.push()
        pin = self._getpinout(pin)
        xs = self._getxs(pin, xs)
        width = self._getwidth(pin, width, xs)
        radius = self._getradius(pin, radius, xs)

        ring = nd.Polygon(points=geom.circle(radius=width/2), layer=self.layer)

        #create points along the mamba for interconnects:
        p1, p2 = [], []
        size = len(points)
        for i in range(size-1):
            dx = points[i+1][0] - points[i][0]
            dy = points[i+1][1] - points[i][1]
            a = np.degrees(m.atan2(dy, dx))
            p1.append((points[i][0], points[i][1], a))
            p2.append((points[i+1][0], points[i+1][1], a+180))

        #print('p1[0]:', p1[0])
        if name is None:
            name = 'mamba'
        with nd.Cell(name, instantiate=False, cnt=True) as ICcell:
            ICcell.default_pins('pla0', 'plb0')
            start = nd.Pin(width=width, xs=xs).put(p1[0])
            nd.Pin('pla0', width=width, xs=xs).put(start.rot(180))

            #loop over guides
            for num, dis in enumerate([pitch*(n-0.5*(N-1))+offset for n in range(N)]):
                last = start.move(0, -dis)
                nd.Pin('a'+str(num), xs=xs).put(last.rot(180))
                if showpins:
                    ring.put(last)
                for i in range(size-2): #loop over points
                    if showpins:
                        pin0 = nd.Pin().put(p1[i+1])
                        ring.put(pin0.move(0, -dis))

                    pin1 = last
                    pin2 = nd.Pin(width=width, xs=xs).put(p2[i+1]).move(0, dis)

                    cp.push()
                    shape = self.__strt_bend_strt_p2p_solve(pin1, pin2, radius)
                    cp.pop()
                    solution, L1, L2, da, message = shape
                    if solution is True:
                        self.strt(length=L1, pin=last).put(last)
                        self.bend(angle=da, radius=radius).put()
                        last = cp.here()
                        if i is size-3:
                            self.strt(length=L2, pin=last).put()
                    else:
                        if i < size-3:
                            self.bend(radius=radius, pin=last, angle=da).put()
                            last = cp.here()
                        else:
                            self.bend_strt_bend_p2p(pin1, pin2, radius=radius).put(last)

                    if showpins:
                        plast = nd.Pin().put(p2[-1])
                        ring.put(plast.move(0, -dis))

                nd.Pin('b'+str(num)).put(cp.here())

            end = nd.Pin(width=width, xs=xs).put(p2[-1])
            nd.Pin('plb0').put(end.rot(180))
            if arrow:
                self.arrow.put(ICcell.pin['pla0'])
                self.arrow.put(ICcell.pin['plb0'])

            if polyline is True:
                if pin is None:
                    nd.Polyline(points=points, width=2, layer=1111).put(0)
                else:
                    nd.Polyline(points=points, width=2, layer=1111).put(p1[0][0], p1[0][1])

        nd.cp.pop()
        return ICcell
