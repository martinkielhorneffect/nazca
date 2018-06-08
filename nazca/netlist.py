#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------
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
# @author: Ronald Broeke (c) 2016-2017, 2011 (c) Xaveer Leijtens,
# @email: ronald.broeke@brightphotonics.eu
#
"""
Nazca classes to construct Nodes, Pointers and Cells,
and from those the graphs/netlists of the design.
"""

import sys
from itertools import count
from collections import defaultdict
import copy as COPY
from math import sin, cos, acos, pi

from numpy.linalg import inv
from numpy import dot
import numpy as np

from . import cfg
from .mask_layers import get_layer
from .geometries import ring


Et = count(0) #counter for default cell names:
elmlist = []
Pt = count(0)

def to_mat(direction=1, chain=0, x=0, y=0, a=0):
    """Transform vector into matrix."""
    a = (a+chain*180.0) / 180. * pi
    if direction == 1: #identity
        return np.array([[cos(a), sin(a), 0.0],
                         [-sin(a), cos(a), 0.0],
                         [x, y, 1.0]])
    else: #inverse
        return np.array([[cos(a), -sin(a), 0.0],
                         [sin(a), cos(a), 0.0],
                         [-x*cos(a)-y*sin(a), x*sin(a)-y*cos(a), 1.0]])

def inverse(mat):
    """Return inverse matrix"""
    return inv(mat)


class Pointer():
    """A pointer with state information.

    The pointer has the positional information of a node.
    The positional information is kept in a matrix that holds the homogeneous
    coordinates. Structures are drawn in local coordinates and are easily
    converted to the global coordinate system by a matrix multiplication. By
    using homogeneous coordinates, translations are also described by a
    matrix multiplication.

    .. figure:: _static/Pointer.png
        :height: 100px
        :alt: alternate text

       *fig: pointer*

    :argument x: x-coordinate of the pointer position [µm]
    :argument y: y-coordinate of the pointer position [µm]
    :argument a: angle of the pointer [degrees]
    """

    def __init__(self, x=0.0, y=0.0, a=0.0):
        """Constructor.

        Args:
            x (float): x coordinate
            y (float): y coordinate
            a (float): angle coordinate
        """
        self.flipstate = False
        self.mat = None
        self.goto(x, y, a)
        self.chain = 0


    def __str__(self):
        xy = self.get_xy()
        return ("%s:\t(x = %.2f, y = %.2f, a = %.2f°)"
                % (self.__class__.__name__, xy[0], xy[1], self.a))


    def goto(self, x, y, a=0.0):  # Absolute
        """
        Move pointer to a coordinate relative to the org.

        Args:
            x (float): x-coordinate of the pointer position [µm]
            y (float): y-coordinate of the pointer position [µm]
            a (float): angle of the pointer [degrees]

        Returns:
            None
        """
        a = a / 180 * pi
        self.mat = np.array([
            [cos(a), sin(a), 0],
            [-sin(a), cos(a), 0],
            [x, y, 1]
            ])
        self.flipstate = False

#        self.goto(ptr.get_xya())
#

    def move(self, x=0, y=0, a=0.0):  # Relative
        """
        Move pointer relative to current location.

        Args:
            x (float): x-coordinate of the pointer position [µm]
            y (float): y-coordinate of the pointer position [µm]
            a (float): angle of the pointer [degrees]

        Returns:
            Pointer: self
        """
        a = a / 180 * pi
        mat = np.array([[cos(a), sin(a), 0.0],
                        [-sin(a), cos(a), 0.0],
                        [x, y, 1.0]])
        self.mat = dot(mat, self.mat)
        return self


    def set_a(self, a=0.0):  # Relative
        """
        Move angle absolute.

        Args:
            a (float): angle of the pointer [degrees]

        Returns:
            Pointer: self
        """
        a = a / 180 * pi
        x = self.x
        y = self.y
        mat = np.array([[cos(a), sin(a), 0.0],
                        [-sin(a), cos(a), 0.0],
                        [x, y, 1.0]])
        self.mat = dot(mat, self.mat)
        return self


    def move_ptr(self, ptr):  # Relative
        """
        Move pointer relative by a pointer 'ptr' to current location.

        Args:
            ptr (Pointer): move by value in pointer

        Returns:
            Pointer: self
        """
        x, y, a = ptr.get_xya()
        a = a / 180 * pi
        mat = np.array([[cos(a), sin(a), 0.0],
                        [-sin(a), cos(a), 0.0],
                        [x, y, 1.0]])
        self.mat = dot(mat, self.mat)
        return self


    def multiply_ptr(self, ptr):  # Relative
        """
        Multiply the pointer by the matrix in <ptr>.

        Args:
            ptr (Pointer): multiply by <ptr>

        Returns:
            Pointer: self
        """
        self.mat = dot(self.mat, ptr.mat)
        return self


    def inv(self):  # Relative
        """Inverse the matrix in the pointer. Returns the pointer.

        Returns:
            Pointer: self
        """
        self.mat = inv(self.mat)
        return self

    def trans(self, mat):  # Relative
        """Translate pointer by matrix <mat>.

        Return:
            Pointer: translated by matrix <mat>
        """
        return dot(mat, self.mat)


    def chain180(self, chain):
        """Set the chain property of the pointer.

        Returns:
            None
        """
        self.chain = chain


    def set_mat(self, t):  # Relative
        """Fill pointer matrix based on vector <t> = (x, y, a).

        Returns:
             None
        """
        self.mat = t


    def skip(self, s):  # move forward (negative: backward)
        """Translate a pointer in the direction it is pointing in.

        Returns:
            Pointer: self
        """
        self.move(s, 0, 0)
        return self


    def offset(self, s):  # offset positive/negative to left/right.
        """Translate a pointer perpendicular to the direction it is pointing in.

        Returns:
            None
        """
        self.move(0, s, 0)


    def flip(self):
        """Flip (mirror) the state of the pointer.

        Returns:
            None
        """
        self.mat = dot(self.mat, [[1, 0, 0], [0, -1, 0], [0, 0, 1]])
        self.flipstate = not self.flip


    def rotate(self, a):
        """
        Rotate pointer by angle a.

        Args:
            a (float): angle of the pointer [degrees]

        Returns:
            Pointer: self
        """
        self.move(0, 0, a)
        return self


    @property
    def x(self):
        """Return position x of the pointer."""
        return self.mat[2, 0]


    @property
    def y(self):
        """Return the position y of the pointer."""
        return self.mat[2, 1]


    @property
    def a(self):
        """Return angle a of the pointer."""
        a = min(1,max(self.mat[0, 0], -1)) #clip [-1, 1]
        if self.mat[0, 1] >= 0:
            return acos(a) * 180 / pi
        else:
            return 360 - acos(a) * 180 / pi


    def get_xy(self, x=0, y=0):
        """Return the (x, y) position of the pointer."""
        return tuple(dot([x, y, 1], self.mat))[0:2]


    def get_xya(self):
        return (self.x, self.y, self.a)


    def get_xy_a(self):
        """Return the pointer position as (x, y, a)."""
        return ((self.x, self.y), self.a)


    def xya(self):
        """Return the pointer position as (x, y, a)."""
        return (self.x, self.y, self.a)


    def copy(self, flip=False):
        """Copy pointer.

        Returns:
            Pointer: copy of the Pointer
        """
        # One level deep copy. The state information is still copied by
        # reference only. This is probably the desired way of working.
        copy = COPY.copy(self)
        if flip:
            copy.y = -copy.y
            copy.a = -copy.a
        return copy


class smatrix:
    """S-matrix tbd."""
    #TE-TE
    #TM-TM
    #TE-TM
    #TM-TE
    #phase
    #amplitude
    def __init__(self, tete=(0, 0), tmtm=(0, 0), tetm=(0, 0), tmte=(0, 0)):
        self.tete = tete
        self.tmtm = tmtm
        self.tetm = tetm
        self.tmte = tmte



class Node():
    """
    Node class for creating node objects.

    The netlist of Nodes and connections between the nodes constructs the photonic
    circuit and or layout.f
    """
    def __init__(self, name=None):
        """Contruct a Node.

        Args:
            name (str): optional Node name.
        """
        #print('______________________new node', name)
        self.nb_geo = [] #nearest neighbors geometry
        self.nb_opt = [] #nearest neighbors, optical connection (in Smatrix)
        self.nb_ele = [] #nearest neighbors, electronic connection
        self.nb_cnode = [] #cell tree, only for used for cnodes
        self.cnode = None # cnode of the node. A cnode is its own cnode
        self.parent_cnode = None # TODO: in use?
        self.cell = None # TODO: in use? Cell object of the cnode
        self.pointer = None # Pointer of the Node. Contains the node position after solving.
        self.xs = None # String name reference to the xs of the node.
        self.width = 0 # width of the node
        self.instance = False # store if node is in a cell (False) or instance (True)
        self.type = None
        self.show = False # show pin in layout.

        if name is None:
            self.name = str(next(Pt))
        else:
            self.name = str(name)


    def __repr__(self):
        coor = "({t[0]:0.3f}, {t[1]:0.3f}, {t[2]:0.3f})".format(t=self.xya())
        try:
            cell = "'{}'".format(self.cnode.cell.cell_name)
        except:
            cell = self.cell
        return "<Node(name='{}') object in cell {}, xs='{}', width={}, xya={}>".format(
            self.name, cell, self.xs, self.width, coor)


    def copy(self):
        """Copy a node.

        Returns:
            Node: copy of the node."""
        # do NOT use deepcopy: extremely slow.
        node = Node()
        node.pointer = self.pointer.copy()
        node.xs = self.xs
        node.width = self.width
        node.instance = self.instance
        node.type = self.type

        cnode = cfg.cells[-1].cnode
        if self.cnode is not cnode and self.cnode.parent_cnode is not cnode:
            raise Exception('Not allowed to create a node outside the scope of max one level deep.')

        node.cnode = cfg.cells[-1].cnode #place node in active cell.
        node.cnode.parent_cnode = cfg.cells[-1].parent_cnode
        return node


    @property
    def chain(self):
        """Attribute to indicate if a Node is a chain connector."""
        return self.pointer.chain


    def goto(self, x=0, y=0, a=0):
        """
        Goto position (<x>, <y>, <a>) with repect to celll origin 'org'.

        Returns:
            Node: a new Node after translating by vector (<x>, <y>, <a>)
        """
        node = Node()
        node.pointer = COPY.copy(self.pointer)
        node.xs = self.xs
        node.width = self.width
        node.type = self.type

        #node.pointer.chain = 0
        node.cnode = cfg.cells[-1].cnode
        node.parent_cnode = cfg.cells[-1].parent_cnode
        connect_geo(cfg.cells[-1].cnode, node, (x, y, a), 1, 0)
        #print('cfg.cells[-1].cnode.cell.cell_name:', cfg.cells[-1].cnode.cell.cell_name)
        cfg.cp = node
        return node


#==============================================================================
#     Create pointer operations on Node level for syntax convenience
#==============================================================================
    def move(self, x=0, y=0, a=0):
        """
        Move pointer in Node by (<x>, <y>, <a>).

        Returns:
            Node: a new Node after translating by vector (<x>, <y>, <a>)
        """
        node = self.copy()
        connect_geo(self, node, (x, y, a), 0, 0)
        #cfg.cp = node
        return node


    def rotate(self, a=0):
        """
        Rotate pointer by <a>.

        Args:
            a (float): angle of rotation in degrees

        Returns:
            Node: a new Node after translating by vector (0, 0, <a>)
        """
        node = self.copy()
        connect_geo(self, node, (0, 0, a), 0, 0, DRC=False)
        #cfg.cp = node
        return node
    rot = rotate


    def skip(self, x=0):
        """
        Skip pointer in Node a distance <x> in direction of the pointer.

        Returns:
            Node: a new Node after translating by vector (<x>, 0, 0)
        """
        node = self.copy()
        connect_geo(self, node, (x, 0, 0), 0, 0 )
        #cfg.cp = node
        return node


    def shift(self, x=0, y=0):
        """
        Shift pointer (<x>, <y>). keep orietation.

        Returns:
            Node: a new Node after translating by vector (<x>, <y,> 0)
        """
        node = self.copy()
        connect_geo(self, node, (x, 0, 0), 0, 0 )
        #cfg.cp = node
        return node


    def offset(self, y=0):
        """
        Offset the pointer of the Node by <y>
        .
        Returns:
            Node: new Node after translating by vector (0, <y>, 0)
        """
        node = self.copy()
        connect_geo(self, node, (0, -y, 0), 0, 0 )
        #cfg.cp = node
        return node
    os = offset


   # def roll(self):
   #     """Roll around the x-axis. (Mirror in the xaxis)"""
   #     node = self.copy()
   #     connect_geo(self, node, (0, 0, 0), 0, 0 )
   #     #cfg.cp = node
   #     return node


    #def mirror(self):
    #    node = self.copy()
    #    connect_geo(self, node, (0, 0, 0), 0, 0 )
    #    #cfg.cp = node
    #    return node


    def xya(self):
        """
        Get pointer position of the Node.

        Returns:
            (x, y, a): pointer position
        """
        return self.pointer.get_xya()


#==============================================================================
# netlist generators
#==============================================================================
    def opt_nb_iter(self):
        """Get optical neighbours.

        Yields:
            Node: next neighbor
        """
        for nn in self.nb_opt:
            if nn[2] == 1:
                yield nn


    def ele_nb_iter(self):
        """Get electrical neighbours.

        Yields:
            Node: next neighbor
        """
        for nn in self.nb_ele:
            if nn[2] == 1:
                yield nn

    def geo_nb_iter(self):
        """Get geometrical neighbours.

        Yields:
            Node: next neighbor
        """
        for nn in self.nb_geo:
            #if nn[2]==1:
            yield nn


    def cnode_nb_iter(self):
        """Get cnodes at the cell level below the cell level of the node.

        Yields:
            Node: next cell

        """
        for nn in self.nb_cnode:
            if nn[2] == 1:
                yield nn[0]


    def print_neighbours(self):
        """Print neighbours of the node in the netlists.

        Returns:
            str: information on neighbors of this node
        """
        out = '\n'.join([
            '---neighbours geo:', '\n'.join(map(str, self.nb_geo)),
            '---neighbours opt:', '\n'.join(map(str, self.nb_opt)),
            '---neighbours ele:', '\n'.join(map(str, self.nb_ele))]
                       )
        return out


    def get_info(self):
        """Get information string.

        Returns:
            str: information on the node.
        """
        s = 'Node---------------------:\n  name: {}\n  cell: {}\n  node: {}'.format(
            self.name,
            self.cnode.cell.cell_name,
            self)
        s += '\n--cnode:\n    name: {}\n    cnode: {}'.format(
            self.cnode.name,
            self.cnode)
        #s += '\n--parent_cnode:\n    cell: {1}\n    cnode: {0}'.format(
        #     self.cnode.parent_cnode,
        #     self.cnode.parent_cnode.cell.cell_name)
        return s



#==============================================================================
# Functions to defining an edge (geometrical, optical, electrical) between nodes.
#==============================================================================

def connect_geo(pin1, pin2, t=(0, 0, 0), chain_node=0, chain_connect=None, DRC=True):
    """Generate a geometrical edge between two Nodes.

    Connects pin n1 and n2 with transformation t.

    Args:
        pin1 (Node): 1st pin
        pin2 (Node): 2nd pin
        chain_node: chain states for n1 and/or n2 if None.
        chain_connect: overrules chain connect rules if None.
            chain_connect = 0 -> chain pin
            chain_connect = 1 -> no-chain pin
            (default = 0)

    chain_pin1 & chain_pin2 = chain_connect

    Returns:
        None
    """

    #print('connect_geo |  n1, n2 =', n1, n2)
    #print('chain_node=',chain_node ,' chain_connect=',chain_connect)

    #create pointer in nodes if none exists.
    #TODO: is this needed? Pointer added on each new pin?
    if pin1.pointer is None:
        pin1.pointer = Pointer()
        pin1.pointer.chain = chain_node

    if pin2.pointer is None:
        pin2.pointer = Pointer()
        pin2.pointer.chain = chain_node

    if chain_connect is None:
        if pin1.cnode is pin2.cnode: # pins in same cell
            chain_connect = 0
        else:
            chain_connect = pin1.pointer.chain or pin2.pointer.chain

    mat = to_mat(1, chain_connect, *t)

    #pin1.nb_geo.append((pin2, mat, 1, chain_connect))
    #pin2.nb_geo.append((pin1, inverse(mat), -1, chain_connect))
    pin1.nb_geo.append((pin2, mat))
    pin2.nb_geo.append((pin1, inverse(mat)))

    DRC = False
    if DRC:
        if abs(t[0]) < 1e-4 and abs(t[1]) < 1e-4: #points at same location
            #print (t, pin1.xs, pin2.xs )
            if pin1.xs is not None and pin2.xs is not None: #points both have a xs
                if pin1 is not pin1.cnode and pin2 is not pin2.cnode:
                    if pin1.xs != pin2.xs: # xs DRC
                        print('DRC: xs mismatch: {} != {}'.format(pin1.xs, pin2.xs))
                        put_polygon(connect=pin1, layer=cfg.drc_xs, points=ring(*cfg.drc_ring_xs))
                    elif abs(t[2]) > 1e-6 and abs(abs(t[2])-180) > 1e-6: #angle DRC
                        print('DRC: non-zero angle connection: {} degrees.'.format(t[2]))
                        put_polygon(connect=pin1, layer=cfg.drc_xs, points=ring(6, 3))

def connect_opt(n1, n2, s):
    """Generate optical connection.

    Returns:
        None
    """
    n1.nb_opt.append((n2, s, 1))
    n2.nb_opt.append((n1, s, -1))


def connect_ele(n1, n2, val):
    """Generate electrical connection.

    Returns:
        None
    """
    n1.nb_ele.append((n2, val, 1))
    n2.nb_ele.append((n1, val, -1))


def connect_cnode(n1, n2, s):
    """Generate cell connection.

    Returns:
        None
    """
    n1.nb_cnode.append((n2, s, 1, 0))
    n2.nb_cnode.append((n1, s, -1, 0))


def parse_pin(C1=None, C2=None, C3=None, C4=None, instance=None, rot=False, default='out'):
    """
    Parse the part of the put statement after a potential first string arg.

    Args:
        C1: tuple of floats, len<=3; C2 = pinname
        C1 to C[i]: 1<=i<=3 are floats; C[i+1] = pinname
        instance (Instance object): Instance object to connect to
        rot (bool): Add an extra rotation for end pins in p2p connections if True

    Returns:
        pin, (x, y, a): describes a point at 'translation' from 'pin'.
            default = 'org', (0, 0, 0).

    Examples:
        Allowed input formats:
        ()
        (Node)
        ((0))
        ((0,0))
        ((0,0,0))
        (0)
        (0,0)
        (0,0,0)

        # 'pin' can be referred to by a Node, Instance or string name.
        (pin)
        (None, pin)
        ((0), pin)
        ((0,0), pin)
        ((0,0,0), pin)
        (0, pin)
        (0,0, pin)
        (0,0,0, pin)
    """

    if C1 is None:
        return cfg.cp, (0, 0, 0)

    P = None
    if isinstance(C1, tuple):
       P = C2
       try:
           x = float(C1[0])
       except:
           x = 0
       try:
           y = float(C1[1])
       except:
           y = 0
       try:
           a = float(C1[2])
       except:
           a = 0

    else:
        try:
            x = float(C1)
            try:
                y = float(C2)
                try:
                    a = float(C3)
                    P = C4
                except:
                    a = 0
                    P = C3
            except:
                y, a = 0, 0
                P = C2
        except:
            x, y, a = 0, 0, 0
            P = C1

    if P is None: #position w.r.t. 'org'
        if instance is None:
            if rot:
                a += 180
            return cfg.cells[-1].cnode, (x, y, a)
        else:
            return instance.cnode.parent_cnode, (x, y, a)

    elif isinstance(P, str):
        return cfg.cells[-1].pin[P], (x, y, a)
    elif isinstance(P, Node):
        return P, (x, y, a)
    elif isinstance(P, Instance):
        if default == 'out':
            return P.pin[P.cnode.cell.default_out], (x, y, a)
        else:
            return P.pin[P.cnode.cell.default_in], (x, y, a)
    else:
        raise Exception('pin could not be parsed {}, {}, {}, {}'.format(
            C1, C2, C3, C4))


class Instance():
    """Class to store pins of cell instantiations.
    """

    def __init__(self):
        """Construct an Instance."""
        self.name = None
        self.cnode = None
        self.parent_cnode = None
        self.pin = dict()
        self.default_in = None
        self.default_out = None
        self.org = None
        self.flip = None
        self.array = None
        #TODO: add instantiate?

    def __repr__(self):
        return "<Instance() object of cell '{}' in cell '{}'>".format(
            self.cnode.cell.cell_name,  self.cnode.parent_cnode.cell.cell_name)


    def raise_pins(self, namesin=None, namesout=None, show=True):
        """Copy pins of a cell instance to its parent cell.

        If a cell A is "put" in cell P, then an instance of cell A is created
        inside "parent" cell P. Each instance of cell A, say A1, A2, ...,
        contains a copy of all pins of cell A to represent cell A in P at a
        particular location, rotation, flip state and/or scaling. raise_pins
        automatically copies all pin attributes.

        The instance its pins are available inside parent P, but the pins are
        themselves not pins of cell P. Pins have to be explicitly set in P
        with nd.Pin.put(). Alternatively, all or part of the pins in an instance
        can be set at once in P by "raise_pins", avoiding many nd.Pin().put()
        statements.

        Note:
            raise_pins is a method of Class Instance, not of Cell. Hence,
            it operates on a cell that has been put.

        Args:
            namesin (list of str): list of pin names to raise
                (Default = all pins of the instance)
            namesout (list of str): list of new names of the raised pins
                (Default = <namesin>)
            show (bool): show the pins in the layout with an arrow (default = True)

        Returns:
            None

        Example:
            Raise (copy) all pins of an instance to its parent cell.

            The version without raise_pins looks like this::

                # put Pin
                import nazca as nd

                with nd.Cell('new_cell') as new_cell:
                    instance1 = nd.example.put(100, 200)
                    nd.Pin('a0').put(instance1.pin['a0'])
                    nd.Pin('b0').put(instance1.pin['b0'])

                instance2 = new_cell.put(0)
                nd.bend(angle=90).put(instance2.pin['a0'])
                nd.export_plt()

            Now with raise_pins on <instance1>,
            which also copies all pin attributes::

                # raise pins
                import nazca as nd

                with nd.Cell('new_cell') as new_cell:
                    instance1 = nd.example.put(100, 200)
                    instance1.raise_pins()

                instance2 = new_cell.put(0)
                nd.bend(angle=90).put(instance2.pin['a0'])
                nd.export_plt()

            Raise and rename pins 'a0' and 'b0' of <instance1>. Note that the
            cell default pins, 'a0' and 'b0', are automatically added to
            <new_cell>::

                # raise and rename pins
                import nazca as nd

                with nd.Cell('new_cell') as new_cell:
                    instance1 = nd.example.put(0)
                    instance1.raise_pins(['a0', 'b0'], ['new_a0', 'new_b0'])

                instance2 = new_cell.put(0)
                nd.bend(angle=90).put(instance2.pin['new_a0'])
                nd.bend(angle=-90).put(instance2.pin['a0'])
                nd.export_plt()
        """

        if namesin is None:
            for name, node in self.pin.items():
                Pin(name, width=node.width, xs=node.xs, type=node.type,
                    chain=node.chain, show=show).put(node)
        else:
            if namesout is None:
                namesout = namesin
            for namein, nameout in zip(namesin, namesout):
                if namein in self.pin.keys():
                    node = self.pin[namein]
                    Pin(nameout, width=node.width, xs=node.xs, type=node.type,
                        chain=node.chain, show=show).put(node)
                else:
                    print("Warning 'raise_pins': pin '{}' not defined.".\
                        format(namein))

num = 0
class Cell():
    """Cell object, corresponding to a GDS cell output (if instantiated).

    A cell has a default pins:
        'org': cell origin (non-chain)
        'a0' : input pin   (chain)
        'b0' : output pin  (chain)

    It is up to the user to place 'a0' and 'b0' at sensible positions in the cell.
    If not provided, they will be place automatically at the origin.

    When starting a new cell: cp = C.pin['org']

    A cell typically contain objects:
        1. pins: A position w.r.t. the cell origin that can be used to connect to.
            method: Pin().put()

        2. cell instances: A reference to another cell
            method: Cell().put()

        3. polygons: A polygon assigned to a specific GDS layer.
            method: Polygon().put()

        4. polylines: A polyline/path assigned to a specific GDS layer.
            method: Polyline().put()

        5. annotations:
            method: Annotation().put()
    """

    def __init__(self, name='cell', celltype='element', instantiate=True,
            hashme=False, cnt=False):
        """Construct a Cell.

        Args:
            name (str): cell name (default = 'cell')
            celltype (str): type of cell, option 'element', 'mask' (default = 'element')
            instantiate (bool): flag if the cell is instantiated (default = True)
        Returns:
            None

        Example:
            Create a cell::

                import nazca as nd

                with nd.Cell('mycell') as C:
                    nd.strt(length=40).put(0)
                    nd.bend(angle=15).put()

                C.put(0)

                nd.export_plt()
        """

        global num
        num += 1

        self.hashme = hashme

        if self.hashme:
            name = cfg.hash_name
            basename = cfg.hash_basename
            paramsname = cfg.hash_paramsname
            self.parameters = cfg.hash_params
        else:
            basename = cfg.gds_cellname_cleanup(name)
            paramsname = cfg.gds_cellname_cleanup(name)
            if cnt:
                name = '{}_{}'.format(name, num)
                name = cfg.gds_cellname_cleanup(name)
            self.parameters = None

        if name in cfg.cellnames.keys():
            if instantiate is True:
                name0 = cfg.gds_cellname_cleanup(name)
                name = '{}_{}'.format(name0, num)
                cfg.cellnames[name] = self
                print("Warning: duplicate cell name '{}' renamed to '{}'.".format(
                    name0, name))
            #else: #Treat as new element
        else:
            name = cfg.gds_cellname_cleanup(name)
            cfg.cellnames[name] = self

        cfg.cells.append(self)
        cfg.self = cfg.cells[-1]
        self.bbinfo = dict() # store metadata on the cell
        self.closed = False
        self.id = next(Et)
        self.cnode = Node(next(Pt))
        self.cnode.cnode = self.cnode
        self.cnode.cell = self
        self.cnode.pointer = Pointer()
        self.cnode.flip = False
        self.cnode.ininstance = False
        self.org = self.cnode
        self.parent_cnode = None

        self.polygons = []
        self.polylines = []
        self.gdsfiles = []
        self.annotations = []
        self.name = self.id
        self.cell_basename = basename     # name: base
        self.cell_paramsname = paramsname # name: base + params
        self.cell_name = name             # name: base + params + hash
        self.instantiate = instantiate
        self.foundry_spt = None # list of BB parameters for use in spt export only.

        #pin handling:
        self.oldcp = cfg.cp
        self.pin = {'org': self.cnode}
        if celltype is 'mask':
            self.default_in = 'org'
            self.default_out = 'org'
        else: #'element':
            self.default_in = 'a0'
            self.default_out = 'b0'

        cfg.cp = self.org


    def __repr__(self):
        return "<Cell(name='{}', instantiate={})>".\
            format(self.cell_name, self.instantiate)


    def default_pins(self, pinin=None, pinout=None):
        """Set default input and output pin.

        Args:
            pinin (Node): pin to set as default cell input pin
            pinout (Node): pin to set as default cell output pin

        Returns:
            None
        """
        if pinin is None:
            pinin = self.default_in
        else:
            self.default_in = pinin

        if pinout is None:
            pinout = self.default_out
        else:
            self.default_out = pinout


    def __enter__(self):
        return self


    def __exit__(self, type, value, traceback):
        self.close()
        return False


    def parse_instance_connection(self, C1, C2, C3, C4, C5, instance=None):
        """Parse up to five connection arguments into a pin1, pin2, (x,y,a) format.

        Args:
            C1 (float | Node | instance): Parsing information
            C2 (float | Node | instance): Parsing information
            C3 (float | Node | instance): Parsing information
            C4 (float | Node | instance): Parsing information
            C5 (float | Node | instance): Parsing information
            instance (Instance): instance object to connect to

        Returns:
            Node, Node, (x, y, a): pin, pin, translation

        Notes:
            The connection has the following syntax:

            [pin_inst [, translation]]

            or

            [pin_inst, ] translation  [, pin_cell]

            * translation = x | x,y | x, y, a | (x) | (x, y) | (x, y, z), where
              translation is w.r.t. <pin_cell>
            * pin_inst: pinname of instance (default = 'a0')
            * pin_cell: pinname of cell (default = 'org')

            Options:

            * ()
            * (node)
            * (instance)
            * (0,)
            * (0, 0)
            * (0, 0, 0)
            * ((0))
            * ((0, 0))
            * ((0, 0, 0))
            * (0, pin)
            * (0, 0, pin_cell)
            * (0, 0, 0, pin_cell)
            * ((0), pin_cell)
            * ((0, 0), pin_cell)
            * ((0, 0, 0), pin_cell)

            * (pin_inst)
            * (pin_inst, 0)
            * (pin_inst, 0, 0)
            * (pin_inst, 0, 0, 0)
            * (pin_inst, (0))
            * (pin_inst, (0, 0))
            * (pin_inst, (0, 0, 0))

            * (pin_inst, pin_cell)
            * (pin_inst, 0, pin)
            * (pin_inst, 0, 0, pin_cell)
            * (pin_inst, 0, 0, 0, pin_cell)
            * (pin_inst, (0), pin_cell)
            * (pin_inst, (0, 0), pin_cell)
            * (pin_inst, (0, 0, 0), pin_cell)
        """

        if C1 is None:
            return instance.pin[self.default_in], cfg.cp, (0, 0, 0)
        elif isinstance(C1, Node):
            if not C1.cnode.ininstance and C1.cnode.cell.closed:
                raise Exception("You are trying to connect to a close Cell object:\n"\
                    "  $ foo.put(cell.pin['a0']\n"\
                    "  Connect to an instance of the cell instead, hence\n"\
                    "  $ instance = cell.put()\n"\
                    "  $ foo.put(instance.pin['a0])\n")
            return instance.pin[self.default_in], C1, (0, 0, 0)
        elif isinstance(C1, Instance):
            if not C1.cnode.ininstance and C1.cnode.cell.closed:
                raise Exception("Trying to connect to a close cell. Try using an instance instead.")
            return instance.pin[self.default_in], C1.pin[C1.cnode.cell.default_out], (0, 0, 0)

        elif isinstance(C1, str):
            pin2, T = parse_pin(C2, C3, C4, C5, instance)
            return instance.pin[C1], pin2, T
        else:
            pin2, T = parse_pin(C1, C2, C3, C4, instance)
            return instance.pin[self.default_in], pin2, T


    def parse_pin_connection(self, C1, C2, C3, C4):
        """
        Parse pin connection polygon, polyline, annotation and gds methods.

        Args:
             C1, C2, C3, C4 (float | Node | instance): Parsing information

        Returns:
            Node, (x, y, a): pin, translation

        Example:
            In case of provinding floats:
            * ()
            * (node)
            * (0,)
            * (0, 0)
            * (0, 0 ,0)
            * ((0))
            * ((0, 0))
            * ((0, 0 ,0))
            * (0, 'b0')
            * (0, 0, 'b0')
            * (0, 0 ,0, 'b0')
            * ((0), 'b0')
            * ((0, 0), 'b0')
            * ((0, 0 ,0), 'b0')
        """

        #TODO: check if pin connection is in scope (parent or sibling)
        if C1 is None or C1 == ():
            return cfg.cp, (0, 0, 0)
        elif isinstance(C1, Node):  #Node
            if not C1.cnode.ininstance and C1.cnode.cell.closed:
                raise Exception("Trying to connect to a close cell. Instantiate the cell first and connect to the instance.")
            return C1, (0, 0, 0)
        elif isinstance(C1, Instance):  #Node
            if not C1.cnode.ininstance and C1.cnode.cell.closed:
                raise Exception("Trying to connect to a close cell. Instantiate the cell first and connect to the instance.")
            return C1.pin[C1.cnode.cell.default_out], (0, 0, 0)

        if isinstance(C1, tuple):
            try:
                x = float(C1[0])
            except:
                x = 0
            try:
                y = float(C1[1])
            except:
                y = 0
            try:
                a = float(C1[2])
            except:
                a = 0
            if isinstance(C2, str):
                return self.pin[C2], (x, y, a)
            else:
                return self.org, (x, y, a)
        else:
            P = None
            try:
                x = float(C1)
                try:
                    y = float(C2)
                    try:
                        a = float(C3)
                        P = C4
                    except:
                        a = 0
                        P = C3
                except:
                    y, a = 0, 0
                    P = C2
            except:
                x, y, a = 0, 0, 0
                P = C1

            if isinstance(P, str):
                return self.pin[P], (x, y, a)
            else:
                return self.org, (x, y, a)


    def put_pin(self, name=None, connect=None, xs=None, width=0, type=None,
            chain=1, show=False):
        """Add new pin and connect it to an existing pin.

        Returns:
            Node: pin being put
        """

        node_new = Node(next(Pt))
        node_new.pointer = Pointer()
        node_new.pointer.chain = chain
        node_new.cnode = self.cnode
        node_new.xs = xs
        node_new.width = width
        node_new.type = type
        node_new.show = show

        if name is not None: # add a node
            self.pin[name] = node_new

        C1, C2, C3, C4 = None, None, None, None
        try: C1 = connect[0]
        except: C1 = connect
        try: C2 = connect[1]
        except: pass
        try: C3 = connect[2]
        except: pass
        try: C4 = connect[3]
        except: pass

        node_exist, T = self.parse_pin_connection(C1, C2, C3, C4)
        connect_geo(node_exist, node_new, t=T, chain_node=1, chain_connect=0)
        return node_new


    def get_pin(self, pinname=None, connect=None):
        """"Parse pinname and connect for put_polygon and put_gds methods.

        Returns:
            Node"""

        if pinname is None:
            pinname = 'org'
        if pinname not in cfg.cells[-1].pin.keys(): # connect to org and add new pin
            if connect is None: # org
                node = self.put_pin(pinname, cfg.cells[-1].org)
            elif isinstance(connect, tuple): #(x,y,a)
                node = self.put_pin(pinname, connect)
                #node = cfg.cells[-1].org.move(*connect)
                #node.name = pinname
            elif isinstance(connect, Node):  #Node
                self.put_pin(pinname, connect)
                node = connect
            return node
        else: #existing pin
            if connect is None: # org
                node = cfg.cells[-1].org
            elif isinstance(connect, tuple): #(x,y,a)
                node = cfg.cells[-1].pin[pinname].move(*connect)
            elif isinstance(connect, Node):  #Node
                node = connect #maybe set parent node !!!!
            return node


    def put_polygon(self, pinname='', connect=(0,0,0),\
            layer=cfg.default_xs['layer'], points=None):
        """Add a polygon to the cell. Anchor it to a pin by the pinname (string).

        Returns:
            None
        """
        #print('put_polygon, pinname = ', pinname)

        layer = get_layer(layer)

        node = self.put_pin(None, connect)
        if points is not None:
            self.polygons.append((node, layer, points))
            #print('polygons: layer:{}:{}'.format(layer, points))
        return None


    def put_polyline(self, pinname=None, connect=None, layer=cfg.default_xs['layer'],
            width=None, pathtype=None, points=None):
        """Add a polyline to the cell. Anchor it to a pin by the pinname (string).

        Returns:
            None"""
        #print('put_polyline, pinname = ', pinname)

        #if isinstance(layer, str):
        layer = get_layer(layer)

        node = self.put_pin(None, connect)
        if points is not None:
            self.polylines.append((node, layer, width, pathtype, points))
            #print('polylines: layer:{}, width:{}, {}'.format(layer, width, points))
        return None


    def put_annotation(self, pinname=None, connect=None, layer=None, text=None):
        """Add an annotation to the cell.

        Args:
            pinname (str): pinname to attach the annotation to.

        Returns:
            None
        """
        #print('put_annotation, pinname = ', pinname)

        #if isinstance(layer, str):
        layer = get_layer(layer)

        node = self.put_pin(pinname, connect)
        if text is not None:
            self.annotations.append((node, layer, text))
            #print('annotation: layer:{}, text:{}'.format(layer, text))
        return None


    def put_gds(self, pinname=None, connect=None, filename=None, cellname=None,
            newcellname=None, layermap=None, cellmap=None, scale=1.0):
        """Add a GDS to the cell. Anchor it to a pin by the pinname (string).

        Example::

            cell_1.put_gds(pinname='a0', filename='path_to_file.gds',
                cellname='name', newcellname='new')
        """
        node = self.put_pin(pinname, connect)
        self.gdsfiles.append((node, filename, cellname, newcellname, layermap,
            cellmap, scale))
        #print('put_gds node {}, pinname {} = '.format(node, pinname))
        return None


    I = 0
    def copy_cell2instance(self, parent_cnode, flip, array):
        """Copy all nodes from a Cell into an Instance.

        Args:
            parent_node (Node): main cell node of the cell to copy.
            flip (bool): flip state of the instance.
            array (list): creates an array instance upon GDS output if not None.
                Format the array as [col#, [dx1, dy1], row#, [dx2, dy2]]
                (default = None)

        Returns:
            Instance: instance created from <parent_node>
        """
        Cell.I += 1 #keep instance counter unique
        instance = Instance()

        try:
            if array is not None:
                A = array
                array = [ A[0], [ A[1][0]*A[0], A[1][1]*A[0] ], \
                          A[2], [ A[3][0]*A[2], A[3][1]*A[2] ] ]
                if not self.instantiate:
                    print(
"Warning: Trying to set array {} on cell '{}' with instantiate = False."\
" Only a single instance will be visible."\
" Set instantiate=True to obtain the full array GDS.".\
                    format(array, self.cell_name))
        except:
            print("Error: incompatible array format {}.".format(array))
            array = None

        # cnode:
        node = Node()
        node.name = '{}_cnode_{}'.format(self.cell_name, Cell.I)
        node.pointer = self.cnode.pointer.copy()
        node.cnode = node
        node.xs = self.cnode.xs
        node.width = self.cnode.width
        node.type = self.cnode.type
        node.cell = self
        node.ininstance = True
        node.flip = flip
        node.array = array
        instance.cnode = node
        instance.cnode.parent_cnode = parent_cnode
        instance.pin['org'] = node
        instance.org = node

        connect_cnode(parent_cnode, instance.cnode, None)

        # all other nodes:
        for key, pin in self.pin.items():
            if key is not 'org':
                node = Node()
                node.xs = pin.xs
                node.width = pin.width
                node.type = pin.type
                node.name = '{}_{}_{}'.format(self.cell_name, key, Cell.I)
                node.cnode = instance.cnode
                instance.pin[key] = node
                # reconstruct geometry and connectivity
                # 1. org-to-node tree
                x, y, a = pin.pointer.get_xya()
                if flip:
                    y, a = -y, -a
                connect_geo(instance.cnode, node, (x, y, a))
                # 2. solved pointer position
                node.pointer = pin.pointer.copy()

        instance.pinin = instance.pin[self.default_in]
        instance.pinout = instance.pin[self.default_out]
        return instance


    def put_instance(self, C1, C2, C3, C4, C5, flip=False, array=None, **kwargs):
        """Instantiate a cell and put it in the active cell.

        Returns:
            Instance: Instance object being put
        """
        parent_cnode = cfg.cells[-1].cnode
        if self.cnode is parent_cnode:
            print('ERROR: can not connect cell \'{}\' to itself.'.format(cfg.cells[-1].cell_name))

        #print('instance flip = {}, {}'.format(flip, self.cell_name))
        instance = self.copy_cell2instance(parent_cnode, flip, array)

        node1, node2, T = self.parse_instance_connection(C1, C2, C3, C4, C5, instance=instance)
        connect_geo(node2, node1, T)

        #update cp:
        p = kwargs.get('newcp', None)
        if p is None:
            p = kwargs.get('cp', None) #backward compatible
        if p is None:
            cfg.cp = instance.pin[self.default_out]
        else:
            cfg.cp = instance.pin[p]

        for ID in cfg._trace_ids:
            cfg.traces[ID].append(instance)

        return instance


    def put(self, *args, flip=False, array=None, **kwargs):
        """Instantiate a Cell object and put it in the active cell.

        Args:
            *args (float | Node | Instance): a set of max. 5 unnamed parameters
                interpreted as a connection between two pins.
                (default = (), which connects to the current pin)
            flip (bool): mirror the instance (default = False)
            cp (str): Set current pin cp to a specific pin after putting (default = 'b0')
                Example: cp='b1'
            array (list): creates an array instance upon GDS output if not None.
                Format the array as [col#, [dx1, dy1], row#, [dx2, dy2]]
                (default = None)

        Returns:
            Instance: reference to the instance of the cell being put

        Examples:
            Below are different syntax examples for put.

            In the examples we use the following objects:
                * nd.example: a pre-defined cell in Nazca
                * 'a0': default input pin name of "example"
                * 'b0': default output pin name of "example"


            1. connect the default pin of cell "example" to current pin (cp)::

                import nazca as nd
                nd.example.put() # -> put('a0', cp)
                nd.export_plt()

            2. connect pin 'b0' of cell "example" to current pin (cp)::

                import nazca as nd
                nd.example.put('b0') # -> put('b0', cp)
                nd.export_plt()

            3. connec default pin of cell "example" to its parent cell at org + (10, 20, 30)::

                import nazca as nd
                nd.example.put(10, 20, 30) # -> put('a0').move(10, 20, 30)
                nd.export_plt()

            4. connect default pin of cell "example" to instance C at pin 'b0'::

                import nazca as nd
                C = nd.example.put(0)
                nd.example.put(C.pin['b0']) # -> put('a0', C.pin['b0'])
                nd.export_plt()

            5. connect pin 'b0' of cell "example" to instance C at pin 'b0'::

                import nazca as nd
                C = nd.example.put(0)
                nd.example.put('b0', C.pin['b0']) # -> put('b0', C.pin['b0'])
                nd.export_plt()

            6. connect pin 'b0' of cell "example" to instance C at its default out pin, i.e.
               refering to instance C only without pin attribute is interpreted as
               connecting to the default output pin of C.::

                import nazca as nd
                C = nd.example.put(0)
                nd.example.put('b0', C) # -> put('b0', C.pin['b0'])
                nd.export_plt()
        """

        C1, C2, C3, C4, C5 = None, None, None, None, None
        try: C1 = args[0]
        except: pass
        try: C2 = args[1]
        except: pass
        try: C3 = args[2]
        except: pass
        try: C4 = args[3]
        except: pass
        try: C5 = args[4]
        except: pass
        return self.put_instance(C1, C2, C3, C4, C5, flip, array, **kwargs)


    def solve(self):
        """Solve the graph of the cell."""
        Netlist().solvecell(self)
        #for k, v in self.pin.items():
        #    print('--', k, v.pointer)


    def __str__(self):
        """Print which structures are contained in the cell."""
        s = 'cell info of \'{}\''.format(self.cell_name)
        s += '\n--pin     {:2}x              = {}'.format(len(self.pin.keys()), self.pin.keys())
        s += '\n--polygon {:2}x (layer, N)   = '.format(len(self.polygons))

        L = []
        for g in self.polygons:
            L.append('({g[1]}, {N})'.format(g=g, N=len(g[2])))
        s += ', '.join(L)
        s += '\n--gds     {:2}x (file, cell) = '.format(len(self.gdsfiles))

        L = []
        for g in self.gdsfiles:
            L.append('({g[1]}, {g[2]})'.format(g=g))
        s += ', '.join(L)

        L = []
        #s += '\n--inst    {:2}x (x,y, cell)  = '.format(len(self.instances))
        #for i in self.instances:
        #    L.append('({}, {})'.format(i[0], i[1].cell_name))
        #s += ', '.join(L)
        return s + '\n'


    def close(self):
        """Close the cell.

        Solve the geometry of the nodes in the Cell.
        Set the cp back to the position before opening the cell.
        If default input and/or output pins have not been set yet they will be
        set on 'org', pointing in opposite directions.

        Returns:
            Cell: self
        """
        #add default ports if missing:
        if self.closed:
            print('Cell \'{}\' already closed.'.format(self.cell_name))

        if self.default_in not in self.pin:
            self.put_pin(self.default_in, (0, 0, 180))
            #print("Warning: no default pin '{}' set in '{}', setting it to {}.".\
            #    format(self.default_in, self.cell_name, (0, 0, 180)))
        if self.default_out not in self.pin:
            self.put_pin(self.default_out, (0))
            #print("Warning: no default pin '{}' set in '{}', setting it to {}.".\
            #    format(self.default_out, self.cell_name, (0, 0, 0)))
            #self.put_pin(self.default_out, cfg.cp)

        self.pinin = self.pin[self.default_in]
        self.pinout = self.pin[self.default_out]
        self.closed = True
        self.solve()
        cfg.cp = self.oldcp
        cfg.cells.pop()
        if cfg.cells:
            cfg.self = cfg.cells[-1]
        return self


def diff(node1, node2):
    """
    Calculate the difference between two nodes.

    Args:
        node1 (Node): start node
        node2 (Node): end node

    Returns:
        (dx, dy, da): difference vector <node2> - <node2> between their pointers
    """
    cfg.cells[-1].solve()
    ptr1 = node1.pointer.copy()
    ptr2 = node2.pointer.copy()
    ptr2.multiply_ptr(ptr1.inv())
    cfg.cells[-1].closeflag = False
    #print(ptr2.get_xya())
    return ptr2.get_xya()


def bbinfo(item=None, data=None):
    """Attach metadata to the Active cell in dictionary bbinfo.'

    Args:
        item (hashable object): dictionary key
        data: The value of the metadata
    """
    cfg.cells[-1].bbinfo[item] = data


def get_transformation(start=None, end=None):
    """Return the transformation matrix that connects two pointers."""
    translate = dot(inverse(start.pointer.mat), end.pointer.rotate(180).mat)
    return translate


#==============================================================================
# Class structures to represent mask structures that can be put.
#==============================================================================

class Pin():
    """Pin class"""
    def __init__(self, name=None, xs=None, width=None, type=None, pin=None, show=False, chain=1):
        """Construct a Pin object.

        Args:
            name (str):
            xs (str): xsection name assigned to the pin.
            width (float): width in um of the connection the pin represents
            type (str): extra proporty for the pin
            pin (Node):
            chain (int):indicate if pin is a chain or connector (1) or not (0)
                (default = 1)
        """
        self.name = name
        self.xs = None
        self.width = None
        self.type = type
        self.pin = pin
        self.chain = chain
        self.show = show

        if pin is not None:
            self.xs = pin.xs
            self.width = pin.width
        if self.xs is None:
            self.xs = xs
        if self.width is None:
            self.width = width

    def __repr__(self):
        return "<Pin() object, name='{}', xs='{}', width={}, type={}, chain={}, show={}>".\
            format(self.name, self.xs, self.width, self.type, self.chain, self.show)

    def put(self, *args):
        if self.pin is not None:
            cfg.cp = self.pin
        return cfg.cells[-1].put_pin(self.name, connect=args, xs=self.xs,
            width=self.width, type=self.type, chain=self.chain, show=self.show)


class Polygon():
    """Polygon class."""
    def __init__(self, pinname=None, layer=None, points=None):
        """Construct a Polygon object.

        Args:
            pinname (str): name of the pin (obsolete)
            layer (int): layer number to put the Polygon in
            points (list): list of points [(x1, y1), (x2, y2), ...]

        Returns:
            None
        """
        self.pinname = pinname
        if layer is None:
            self.layer = cfg.default_xs['layer'][0]
        else:
            self.layer = layer
        self.points = points

    def __repr__(self):
        try:
            size = len(self.points)
        except:
            size = 0;
        return "<Polygon() object, layer={}, points:{}>".format(self.layer, size)

    def put(self, *args):
        return cfg.cells[-1].put_polygon(self.pinname, connect=args,
            layer=self.layer, points=self.points)


class Polyline():
    """Polyline (path) class."""

    def __init__(self, pinname=None, layer=None, width=None, pathtype=0,
                points=None):
        """Construct a Polyline object.

        Args:
            pinname (str): name of the pin (obsolete)
            layer (int): layer number to put the Polyline in
            width (float): width of the polyline
            pathtype (int): gds type of path: 1, 2, 3 or 4
            points (list): list of points [(x1, y1), (x2, y2), ...]

        Returns:
            None
        """
        self.pinname = pinname
        self.layer = layer
        if width is not None:
            self.width = width
        else:
            self.width = 0.20
        self.pathtype = pathtype
        self.points = points

    def __repr__(self):
        try:
            size = len(self.points)
        except:
            size = 0
        return "<Polyline() object, layer={}, width={}, pathtype={} points:{}>".\
            format(self.layer, self.width, self.pathtype, size)

    def put(self, *args):
        return cfg.cells[-1].put_polyline(self.pinname, connect=args,
            layer=self.layer, width=self.width, pathtype=self.pathtype,
            points=self.points)


class Annotation():
    """Annotation class."""
    def __init__(self, pinname=None, layer=None, text=''):
        """Construct an annotation object.

        Args:
            pinname (str): name of the pin (obsolete)
            layer (int): layer number to put the Polyline in
            text (str): annotation text

        Returns:
            None
        """
        self.pinname = pinname
        self.layer = layer
        self.text = text

    def __repr__(self):
        if len(self.text) >= 10:
            text = self.text[:7]+'...'
        else:
            text = self.text
        return "<Annotation() object, layer={}, text='{}'>".format(
            self.layer, text)

    def put(self, *args):
        return cfg.cells[-1].put_annotation(self.pinname, connect=args,
            layer=self.layer, text=self.text)



#==============================================================================
# Functions to place mask structures in the active cell.
# These become obsolete
# The class alternatives are used for put-syntax consistency with Cell
#==============================================================================

def put_pin(name=None, connect=None, xs=None, width=None, type=None, chain=1):
    """put_pin in the active cell."""
    return cfg.cells[-1].put_pin(name, connect, xs, width, type, chain)

def put_polygon(pinname=None, connect=None, layer=cfg.default_xs['layer'][0], points=None):
    """put_polygon in the active cell."""
    return cfg.cells[-1].put_polygon(pinname, connect, layer, points)

def put_polyline(pinname=None, connect=None, layer=None, width=None,
        pathtype=0, points=None):
    """put_polyline in the active cell."""
    return cfg.cells[-1].put_polyline(pinname, connect, layer, width, pathtype, points)

def put_annotation(pinname=None, connect=None, layer=None, text=None):
    """put_annotation in the active cell."""
    return cfg.cells[-1].put_annotation(pinname, connect, layer, text)

def put_gds(pinname=None, connect=None, filename=None, cellname=None,
        newcellname=None, layermap=None, cellmap=None, scale=1.0):
    """put_polygon inthe active cell."""
    return cfg.cells[-1].put_gds(pinname, connect, filename, cellname,
        newcellname, layermap, cellmap, scale)




import nazca.gds_import as gdsimp

def __get_cellnames(filename, cellname):
    """Get a list of all filenames.

    Args:
        filename: gds filename

    Returns:
        list of str: list of cellnames in <filename> under cell name <cellname>
    """
    strm = gdsimp.GDSII_stream(filename)
    return strm.cell_branch(cellname)
    #return strm.cells[cellname].snames

def gds_filter(filename, cellmap=None, layermap=None):
    """Filter layers and cells form gds and export a copy to file."""
    #TODO filter out cellnames by pattern.
    print('Filtering \'{}\'...'.format(filename))
    strm = gdsimp.GDSII_stream(filename=filename, cellmap=cellmap, layermap=layermap)
    newfilename = '{}_filtered.gds'.format(filename[:-4])
    strm.GDSII_write(newfilename)
    print('...Wrote \'{}\''.format(newfilename))

def print_structure(filename, cellname, level=0):
    """Print gds structure in ascii art.

    Rerturns:
        None
    """
    strm = gdsimp.GDSII_stream(filename)
    strm.print_structure(cellname, level)


def load_gds(filename, cellname, newcellname=None, layermap=None,
        cellmap=None, scale=1.0, prefix='', instantiate=False):
    """Load a GDS cell <celllname> (and its instances) from <filename> into a Nazca cell.

    This function can prefix cell names in <filename> to avoid name classes between
    present cells and cells in <filename>.

    Args:
        filename (str): filename of gds to load
        cellname (str): cellname in filename to load
        newcellname (str): new name for <cellname>
        layermap (dict): mapping {old_layernumber: new_layernumber}
        cellmap (dict): mapping {old_cellname: new_cellname}
        scale (float): scaling factor if the cell (not yet implemented)
        prefix (str): optional string to avoid name clashes
        instantiate (bool): instantiate the GDS (default = False)

    Returns:
        Cell: Cell containing the loaded gds cell(tree)

    Example:
        gds_load(filename='path_to_file.gds', cellname='name')
    """

    #print('load_gds: {} -- {}'.format(filename, cellname))
    global num
    num += 1

    if cellmap is None:
        cellmap = dict()
        if newcellname is not None:
            cellmap[cellname] = newcellname

    #check if cellnames in loaded GDS cell tree clash with existing cell names.
    names = __get_cellnames(filename, cellname)
    for name in sorted(names):
        if name not in cellmap:
            cellmap[name] = '{}{}'.format(prefix, name)
        if cellmap[name] in cfg.cellnames.keys():
            raise Exception(
"""Error in load_gds: cell name '{0}' in file '{1}' is already in use in the design.
  You have five options to solve this issue:
    1. give a unique cell name, i.e. load_gds(..., newcellname=<uniquename>
    2. rename cell '{0}' of the existing cell.
    3. rename cell '{0}' in file '{1}'.
  and for gds trees:
    4. apply a cellmap dict. i.e. load_gds(..., cellmap=<cellmap_dict>).
    5. apply a prefix, i.e. load_gds(..., prefix='new').
""".format(cellmap[name], filename, prefix))
        cfg.cellnames[cellmap[name]] = 'loaded_gds_cell'

    with Cell('load_gds', celltype='mask', instantiate=instantiate) as maskcell:
        put_gds(pinname='org',
            connect=(0, 0, 0),
            filename=filename,
            cellname=cellname,
            newcellname=cellmap[cellname],
            layermap=layermap,
            cellmap=cellmap,
            scale=scale)
    return maskcell


def show_cp(r=10, w=1):
    Polygon(layer=1, points=ring(radius=r, width=w)).put()


class Netlist:
    """Class for building the netlist."""

    def __init__(self):
        """ Construct a Netlist object. Singleton.
        """
        self.nodes_visited = set()
        self.cnodes_visited = set()


    def solvecell_core(self, node, cnode0=None):
        """
        Calculate the geometrical positions of nodes in a cell via the netlist.

        This includes all nodes in all instantiated cells one level deep.
        A starting node 'node' has to be provided.
        cnode0 acts as the reference to the cell being solved.
        A good starting node is a cnode with its pointer set at (0, 0, 0).

        Nodes that are being processed, in scope, adhere to the following condition:
        - the node resides in same cell as the starting node (parent cnode0).
        - the node resides in an instance of the cell of the starting node.

        Args:
            node (Node): start node for solving the cell.
            cnode0 (Node): cnode of the cell to solve (automatically obtained from node).

        Returns:
            None
        """

        if cnode0 is None:
            cnode0 = node.cnode
        self.nodes_visited.add(node)

        neighbours = node.geo_nb_iter()
        for nextnode, mat in neighbours:
            if nextnode not in self.nodes_visited:
                #print(node)
                #node in main cell or in instance
                if (nextnode.cnode is cnode0) or \
                        (nextnode.cnode.parent_cnode is cnode0):
                    nextnode.pointer.set_mat(node.pointer.trans(mat))
                    self.solvecell_core(nextnode, cnode0)
                else:
                    print("...Warning: skipping node in cell '{}': it is not in scope of cell '{}'.".format(
                        nextnode.cnode.cell.cell_name, cnode0.cell.cell_name))

            #TODO: else: check geo consistency


    def solvecell(self, cell):
        """Set cell origin at (0, 0, 0) and solve the cell."""
        #print('SOLVING cell: {}'.format(cell.cell_name))
        cell.cnode.pointer.goto(0, 0, 0)
        self.solvecell_core(cell.cnode)
        return None


    def polygon_iter(self, cnode, infolevel=0):
        """Iterater for all polygons in a cell (cnode).

        Yields:
            pointer, layer, points:
                Next Polygon and info to generate it
        """
        polygons = cnode.cell.polygons
        for polygon in polygons:
            node, layer, points = polygon
            if infolevel > 1:
                print('{}polygon: ML={}, xya={}'.format('  '*infolevel,
                    layer, node.pointer.get_xya()))
            yield node.pointer.copy(), layer, points


    def polyline_iter(self, cnode, infolevel=0):
        """Iterator for all polylines (paths) in a cell (cnode).

        Yields:
            pointer, width, layer, path type points:
                Next Polyline and info to generate it

        """
        polylines = cnode.cell.polylines
        for polyline in polylines:
            node, layer, width, pathtype, points = polyline
            if infolevel > 1:
                print('{}polyline: ML={}, xya={}'.format('  '*infolevel,
                    layer, node.pointer.get_xya()))
            yield node.pointer.copy(), width, layer, pathtype, points


    def annotation_iter(self, cnode, infolevel=0):
        """Iterator for all annotations in a cell (cnode).

        Yields:
            pointer, layer, text:
                Next Annotation and info to generate it
        """
        annotations = cnode.cell.annotations
        for annotation in annotations:
            node, layer, text = annotation
            if infolevel > 1:
                print('{}annotation: ML={}, xya={}'.format('  '*infolevel,
                    layer, node.pointer.get_xya()))
            yield node.pointer.copy(), layer, text


    def gdsfile_iter(self, cnode, level=0, infolevel=0):
        """Iterator for all GDS files instantiated in a cell (cnode).

        Yields:
            pointer, file, cell, newcell, layermap, cellmap, scale:
                Next gds file and info to generate it
        """
        gdsfiles = cnode.cell.gdsfiles
        for gdsinfo in gdsfiles:
            pin, file, cell, newcell, layermap, cellmap, scale = gdsinfo
            out = pin.pointer.copy(), file, cell, newcell, layermap, cellmap, scale
            if infolevel > 1:
                print('{}gds-cell: cell={}, xya={}'.format('  '*infolevel,
                    newcell, cnode.pointer.get_xya()))
            yield out


    def celltree_iter(self, cnode, level=0, position=None, infolevel=0):
        """
        Iterate over the cell hierarchy from the topcell at 'cnode'.
        Deep-first search.

        The tree decents from a 'cell cnode' into an 'instance cnode',
        To subsequently decent into an instance we
        - look up the cell object the instance is derived from,
        - loop up the 'cell cnode' of that cell object,
        - use that cnode to yield and seed the next cell level.

        Yields:
            cnode, level, position, flip:
                Yields next cell, its level in the hierarchy and position.
        """

        flip = cnode.flip
        cnode = cnode.cell.cnode #get cell cnode for cell or instance
        self.cnodes_visited.add(cnode)

        if position is None: #coordinate
            position = Pointer()
        if infolevel > 1:
            print('{}CT: cnode={}, name={}, level={}'.format(
                '  '*infolevel, cnode, cnode.cell.cell_name, level))

        yield cnode, level, position, flip

        cnode_iter = cnode.cnode_nb_iter()
        for nextcnode in cnode_iter:
            if cnode.cell is not None:
                elm1 = cnode.cell.name
            else:
                elm1 = 'None'
            if nextcnode.cell is not None:
                elm2 = nextcnode.cell.name
            else:
                elm2 = 'None'
            if elm1 is not elm2:
                yield from self.celltree_iter(
                    nextcnode,
                    level+1,
                    nextcnode.pointer,
                    infolevel)

#-----------------------------------------------------------------------------
if __name__ == '__main__':
    import nazca.layout as layout
# Don't want to make each pin into a node.
    pin_ph = Cell(name='pin_ph')
    pin_ph.put_gds(filename='pin_symbol.gds', cellname='pin_ph')
    p = pin_ph.close()
    # The cell
    C = Cell(name='one')
    C.put_polygon(layer=5, points=[[10, 10], [20, 10], [20, 30], [10, 30]])
    C.put_pin(name='a0', connect=(10, 20, 180))
    C.put_pin(name='b0', connect=(12, 10, -90))
    C.put_pin(name='b1', connect=(15, 10, -90))
    C.put_pin(name='b2', connect=(18, 10, -90))
    C.put_pin(name='c0', connect=(20, 20, 0))
    C.put_gds(pinname='a0', filename='pin_symbol.gds', cellname='pin_ph')
    # Adding just a reference to another cell is missing?
    # Repeat for now
#    p.put()
    C.put_gds(pinname='b0', filename='pin_symbol.gds', cellname='pin_ph')
    C.put_gds(pinname='b1', filename='pin_symbol.gds', cellname='pin_ph')
    C.put_gds(pinname='b2', filename='pin_symbol.gds', cellname='pin_ph')
    C.put_gds(pinname='c0', filename='pin_symbol.gds', cellname='pin_ph')
    C.close()

    layout.export_gds(C)
