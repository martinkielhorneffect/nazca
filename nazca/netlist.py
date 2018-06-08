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

from itertools import count
from collections import defaultdict
import copy as COPY
from math import sin, cos, acos, pi

from numpy.linalg import inv
from numpy import dot
import numpy as np
from pprint import pprint

from . import cfg
from .mask_layers import get_layer
from .geometries import ring
from . import gds_import as gdsimp
from . import gds_base as gb
import nazca.trace as trace
from nazca.pdk_template_core import put_boundingbox, bbox_pinnames

Et = count(0) #counter for default cell names:
elmlist = []
Pt = count(0)

def to_mat(chain=0, x=0, y=0, a=0, inverse=False):
    """Transform vector into matrix."""
    a = (a+chain*180.0) / 180. * pi
    if not inverse:
        return np.array([[cos(a), sin(a), 0.0],
                         [-sin(a), cos(a), 0.0],
                         [x, y, 1.0]])
    else:
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
        """Move pointer to a coordinate relative to the org.

        Args:
            x (float): x-coordinate of the pointer position [µm]
            y (float): y-coordinate of the pointer position [µm]
            a (float): angle of the pointer [degrees]

        Returns:
            None
        """
        a = pi*a/180
        self.mat = np.array([
            [cos(a), sin(a), 0],
            [-sin(a), cos(a), 0],
            [x, y, 1]
            ])
        self.flipstate = False

#        self.goto(ptr.get_xya())
#

    def move(self, x=0, y=0, a=0.0):  # Relative
        """Move pointer relative to current location.

        Args:
            x (float): x-coordinate of the pointer position [µm]
            y (float): y-coordinate of the pointer position [µm]
            a (float): angle of the pointer [degrees]

        Returns:
            Pointer: self
        """
        a = pi*a/180
        mat = np.array([[cos(a), sin(a), 0.0],
                        [-sin(a), cos(a), 0.0],
                        [x, y, 1.0]])
        self.mat = dot(mat, self.mat)
        return self


    def set_a(self, a=0.0):
        """Set angle absolute. Keep x, y position

        Args:
            a (float): angle of the pointer [degrees]

        Returns:
            Pointer: self
        """
        a = pi*a/180
        x = self.x
        y = self.y
        self.mat = np.array([[cos(a), sin(a), 0.0],
                             [-sin(a), cos(a), 0.0],
                             [x, y, 1.0]])
        #self.mat = dot(mat, self.mat)
        return self


    def move_ptr(self, ptr):  # Relative
        """Move pointer relative by a pointer 'ptr' to current location.

        Args:
            ptr (Pointer): move by value in pointer

        Returns:
            Pointer: self
        """
        x, y, a = ptr.get_xya()
        a = pi*a/180
        mat = np.array([[cos(a), sin(a), 0.0],
                        [-sin(a), cos(a), 0.0],
                        [x, y, 1.0]])
        self.mat = dot(mat, self.mat)
        return self


    def multiply_ptr(self, ptr):  # Relative
        """Multiply the pointer by the matrix in <ptr>.

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
        self.flipstate = not self.flipstate


    def rotate(self, a):
        """Rotate pointer by angle a.

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
            return acos(a)*180/pi
        else:
            return 360 - acos(a)*180/pi


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
    """Node class for creating node objects.

    The netlist of Nodes and connections between the nodes (edges) construct
    the photonic circuit and or layout.
    """
    def __init__(self, name=None):
        """Contruct a Node.

        Args:
            name (str): optional Node name (default is an interger counter)

        Returns:
            None
        """
        #print('______________________new node', name)
        self.id = next(Pt)
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
        self.remark = None

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
        return "<Node(id={}, name='{}') object in cell {}, xs='{}', width={}, xya={}, remark='{}'>".\
            format(self.id, self.name, cell, self.xs, self.width, coor, self.remark)


    def copy(self):
        """Copy a Node object.

        Returns:
            Node: copy of the node with the same attribute settings as the original"""
        # do NOT use deepcopy: extremely slow.
        node = Node()
        node.pointer = self.pointer.copy()
        node.xs = self.xs
        node.width = self.width
        node.instance = self.instance
        node.type = self.type
        node.remark = self.remark

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
        """Goto position (<x>, <y>, <a>) with repect to celll origin 'org'.

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
        """Move pointer in Node by (<x>, <y>, <a>).

        Returns:
            Node: a new Node after translating by vector (<x>, <y>, <a>)
        """
        node = self.copy()
        connect_geo(self, node, (x, y, a), 0, 0)
        #cfg.cp = node
        return node


    def rotate(self, a=0):
        """Rotate pointer by <a>.

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
        """Skip pointer in Node a distance <x> in direction of the pointer.

        Returns:
            Node: a new Node after translating by vector (<x>, 0, 0)
        """
        node = self.copy()
        connect_geo(self, node, (x, 0, 0), 0, 0 )
        #cfg.cp = node
        return node


    def shift(self, x=0, y=0):
        """Shift pointer (<x>, <y>) and keep orientation.

        Returns:
            Node: a new Node after translating by vector (<x>, <y,> 0)
        """
        node = self.copy()
        connect_geo(self, node, (x, 0, 0), 0, 0 )
        #cfg.cp = node
        return node


    def offset(self, y=0):
        """Offset the pointer of the Node by <y>.

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
        """Get pointer position of the Node.

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
def connect_geo(pin1, pin2, t=(0, 0, 0), chain_node=0, chain_connect=None,
        DRC=False, solve=True):
    """Generate a geometrical edge between two Nodes.

    Connects pin1 and pin2 with transformation geometrical t.

    Args:
        pin1 (Node): 1st pin
        pin2 (Node): 2nd pin
        t (tuple): edge value: transformation (x, y, a) to get from pin1 to pin2
        chain_node (int): 1 if a chain node, 0 if not (default = 0)
        chain_connect: overrules chain connect rules if None.
            chain_connect = 0 -> chain pin,
            chain_connect = 1 -> no-chain pin,
        (default = None)
        DRC (bool): apply DRC on the edge (default = True)


    Returns:
        None
    """
    #create a pointer in the nodes if none exists.
    assert pin1.pointer is not None
    #print("1: {}, 2: {}".format(pin1.pointer, pin2.pointer))
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

    mat = to_mat(chain_connect, *t)
    pin1.nb_geo.append((pin2, mat))
    pin2.nb_geo.append((pin1, inverse(mat)))

    if cfg.solve_direct and solve:
        if pin2.cnode is pin1.cnode:
            pin2.pointer.set_mat(pin1.pointer.trans(mat))
        elif pin2.cnode.parent_cnode is pin1.cnode:
            # case where instance pin is created from cell pin in instance creation
            pin2.pointer.set_mat(pin1.pointer.trans(mat))
        elif pin1.cnode.parent_cnode is pin2.cnode:
            # case for raise pin where cell pin is created from instance pin
            pin2.pointer.set_mat(pin1.pointer.trans(mat))
        elif pin1.cnode.parent_cnode is pin2.cnode.parent_cnode:
            # case for instance pin connected together in same parent
            pin2.pointer.set_mat(pin1.pointer.trans(mat))
        else:
            raise Exception("Error: pin connection not in scope.")

    if DRC:
        if abs(t[0]) < 1e-4 and abs(t[1]) < 1e-4: #points at same location
            #print (t, pin1.xs, pin2.xs )
            if pin1.xs is not None and pin2.xs is not None: #points both have a xs
                if pin1 is not pin1.cnode and pin2 is not pin2.cnode:
                    if pin1.xs != pin2.xs: # xs DRC
                        print('DRC: xs mismatch: {} != {}'.format(pin1.xs, pin2.xs))
                        Polygon(points=ring(*cfg.drc_ring_xs), layer=cfg.drc_xs, ).put(pin1)
                    elif abs(t[2]) > 1e-6 and abs(abs(t[2])-180) > 1e-6: #angle DRC
                        print('DRC: non-zero angle connection: {} degrees.'.format(t[2]))
                        Polygon(points=ring(6, 3), layer=cfg.drc_xs).put(pin1)


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
        raise Exception('pin could not be parsed {}, {}, {}, {}'.\
            format(C1, C2, C3, C4))


def validate_basename(basename):
    """Check if basename is unique on substring level.

    This is only relevant for pcell replacement where no mistake between
    different basenames is allowed.

    Args:
        basename (str): basename of a cell's name

    Returns:
        bool: True if basename is valid
    """

    size1 = len(basename)
    for bname in cfg.basenames:
        size2 = len(bname)
        size = min(size1, size2)
        #print(bname[0:size], basename[0:size])
        if bname[0:size] == basename[0:size]:
            print("Error: basename overlap {} - {}".format(bname, basename))
            return False
    return True


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

    @property
    def length_geo(self):
        try:
            lgeo = self.cnode.cell.length_geo
        except:
            lgeo =0

        return lgeo

    def __repr__(self):
        return "<Instance() object of cell '{}' in cell '{}'>".format(
            self.cnode.cell.cell_name,  self.cnode.parent_cnode.cell.cell_name)

    def ic_pins(self):
        """Generator over interconnect pins, filter out org and bounding box.

        Yields:
            str, Pin: iterator over pins in cell
        """
        for name, p in self.pin.items():
            if (name != 'org') and (name not in bbox_pinnames):
                yield name, p


    def raise_pins(self, namesin=None, namesout=None, show=True):
        """Copy pins of a cell instance to its parent cell.

        If a cell A is "put" in cell P, then an instance of cell A is created
        inside "parent" cell P. Each instance of cell A, say A1, A2, ...,
        contains a copy of all pins of cell A to represent cell A in P at a
        particular location, rotation, flip state and/or scaling. raise_pins
        automatically copies all pin attributes.

        The instance pins are themselves not pins of cell P.
        Pins have to be explicitly set in P with nd.Pin.put().
        Alternatively, all or part of the pins in an instance
        can be set at once in P by "raise_pins", avoiding many nd.Pin().put()
        statements. This can very useful when an instance has a large number
        of pins that have to be "raised".

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
            In the example a 90 degree bend is connected to port 'a0'
            of the new_cell.

            Create pins without the `raise_pins` method looks like this::

                # create and put Pins explicitly one by one:
                import nazca as nd

                with nd.Cell('new_cell') as new_cell:
                    instance1 = nd.example_cell().put(100, 200)
                    nd.Pin('a0').put(instance1.pin['a0'])
                    nd.Pin('b0').put(instance1.pin['b0'])

                instance2 = new_cell.put(0)
                nd.bend(angle=90).put(instance2.pin['a0'])
                nd.export_plt()

            Now **raise** pins with `raise_pins` on <instance1> in the new_cell.
            This also copies all pin attributes::

                # raise pins
                import nazca as nd

                with nd.Cell('new_cell') as new_cell:
                    instance1 = nd.example_cell().put(100, 200)
                    instance1.raise_pins()

                instance2 = new_cell.put(0)
                nd.bend(angle=90).put(instance2.pin['a0'])
                nd.export_plt()

            Now **raise and rename** pins 'a0' and 'b0' of <instance1>
            with method `raise_pins`. Note that the
            cell default pins, 'a0' and 'b0', are automatically added to
            <new_cell>::

                # raise and rename pins
                import nazca as nd

                with nd.Cell('new_cell') as new_cell:
                    instance1 = nd.example_cell().put(0)
                    instance1.raise_pins(['a0', 'b0'], ['new_a0', 'new_b0'])

                instance2 = new_cell.put(0)
                nd.bend(angle=90).put(instance2.pin['new_a0'])
                nd.bend(angle=-90).put(instance2.pin['a0'])
                nd.export_plt()
        """

        if namesin is None:
            for name, node in self.pin.items():
                Pin(name, width=node.width, xs=node.xs, type=node.type,
                    chain=node.chain, show=show, remark=node.remark).put(node)
        else:
            if namesout is None:
                namesout = namesin
            for namein, nameout in zip(namesin, namesout):
                if namein in self.pin.keys():
                    node = self.pin[namein]
                    Pin(nameout, width=node.width, xs=node.xs, type=node.type,
                        chain=node.chain, show=show, remark=node.remark).put(node)
                else:
                    print("Warning 'raise_pins': pin '{}' not defined.".\
                        format(namein))

num = 0
class Cell():
    """Cell object, corresponding to a GDS cell when eported to GDS.

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
            autobbox=False, hashme=False, cnt=False):
        """Construct a Cell.

        Args:
            name (str): cell name (default = 'cell')
            celltype (str): type of cell, option 'element', 'mask' (default = 'element')
            instantiate (bool): flag if the cell is instantiated (default = True)
            hashme (bool): set True if cell obtains information from hasme decorator
            cnt (bool): Append ordinal counter to cellname (default = False)

        Returns:
            None

        Example:
            Create a cell with cell name 'mycell' and variable name C.
            Put the cell in the mask layout and export the layout::

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
            self.func_id = cfg.hash_func_id
        else:
            basename = cfg.gds_cellname_cleanup(name)
            if cfg.validate_basename:
                validate_basename(basename)
            paramsname = cfg.gds_cellname_cleanup(name)
            if cnt:
                name = '{}_{}'.format(name, num)
                name = cfg.gds_cellname_cleanup(name)
            self.parameters = None
            self.func_id = None

        name = cfg.gds_cellname_cleanup(name)
        if name in cfg.cellnames.keys():
            #TODO: what if new name also exists?
            name2 = '{}_{}'.format(name, num)
            cfg.cellnames[name2] = self
            if instantiate is True:
                print("Warning: duplicate cell name '{}' renamed to '{}'.".\
                    format(name, name2))
            name = name2
        else:
            cfg.cellnames[name] = self

        cfg.cells.append(self)
        cfg.self = cfg.cells[-1]
        self.bbinfo = dict() # store metadata on the cell
        self.closed = False
        self.id = next(Et)
        self.cnode = Node(next(Pt))
        self.cnode.cnode = self.cnode
        self.cnode.cell = self
        self.cnode.pointer = Pointer(0, 0, 0)
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
        self.bbox = None
        self.autobbox = autobbox

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


    def filter(layers):
        """Create a new cell after applying a layer filter.

        Args:
            layers (list of layer): layers to keep

        Returns:
            Cell: new cell with subset of original layers

        To be implemented
        """


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

    def ic_pins(self):
        """Generator over interconnect pins, filter out org and bounding box.

        Yields:
            str, Pin: iterator over pins in cell
        """
        for name, p in self.pin.items():
            if (name != 'org') and (name not in bbox_pinnames):
                yield name, p

    def __enter__(self):
        return self


    def __exit__(self, type, value, traceback):
        self.close()
        return False


    def parse_instance_connection(self, C1, C2, C3, C4, C5, instance=None):
        """Connect a (closed) cell as an Instance object into the active cell.

        Parse up to five connection arguments to extract: pin1, pin2 and (x,y,a).

        Args:
            C1 (float | Node | Instance): connection information
            C2 (float | Node | Instance): connection information
            C3 (float | Node | Instance): connection information
            C4 (float | Node | Instance): connection information
            C5 (float | Node | Instance): connection information
            instance (Instance): instance object to connect to

        Returns:
            Node, Node, (x, y, a): instance pin, cell pin, connection

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
        """Parse pin connection polygon, polyline, annotation and gds methods.

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


    def _put_pin(self, name=None, connect=None, xs=None, width=0, type=None,
            chain=1, show=False, remark=None):
        """Add new pin and connect it to an existing pin.

        Returns:
            Node: pin being put
        """
        node_new = Node(name) #next(Pt)
        node_new.pointer = Pointer()
        node_new.pointer.chain = chain
        node_new.cnode = self.cnode
        node_new.xs = xs
        node_new.width = width
        node_new.type = type
        node_new.show = show
        node_new.remark = remark

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
                node = self._put_pin(pinname, cfg.cells[-1].org)
            elif isinstance(connect, tuple): #(x,y,a)
                node = self._put_pin(pinname, connect)
                #node = cfg.cells[-1].org.move(*connect)
                #node.name = pinname
            elif isinstance(connect, Node):  #Node
                self._put_pin(pinname, connect)
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


    def _put_polygon(self, connect=None, polygon=None):
        """Add a polygon to the cell.

        Args:
            connect (float, float, float): position (x, y, a) of polygon in cell.
            polygon (Polygon)

        Returns:
            None
        """
        node = self._put_pin(None, connect)
        if polygon is not None:
            self.polygons.append((node, polygon))
        return None


    def _put_polyline(self, connect=None, polyline=None):
        """Add a polyline to the cell.

        Args:
            connect (float, float, float): position (x, y, a) of polygon in cell.
            polyline (Polyline)

        Returns:
            None"""
        node = self._put_pin(None, connect)
        if polyline is not None:
            self.polylines.append((node, polyline))
        return None


    def _put_annotation(self, connect=None, annotation=None):
        """Add an annotation to the cell.

        Args:
            pinname (str): pinname to attach the annotation to.

        Returns:
            None
        """
        node = self._put_pin(None, connect)
        if annotation is not None:
            self.annotations.append((node, annotation))
        return None


    def _put_gds(self, connect=None, filename=None, cellname=None,
            newcellname=None, layermap=None, cellmap=None, scale=1.0, strm=None):
        """Add a GDS to the cell. Anchor it to a pin by the pinname (string).

        Example::

            cell_1.put_gds(filename='path_to_file.gds',
                cellname='name', newcellname='new')
        """
        node = self._put_pin(connect=connect)
        self.gdsfiles.append((node, filename, cellname, newcellname, layermap,
            cellmap, scale, strm))
        #print('put_gds node {}, pinname {} = '.format(node, pinname))
        return None


    I = 0
    def _copy_cell2instance(self, parent_cnode, flip, flop, array):
        """Copy all nodes from the Cell into an Instance.

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

        flipflop = flip != flop
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
            print("Error: incompatible array format {} when placing cell '{}' in '{}'."\
                " Array will be ignored.".\
                format(array, self.cell_name, parent_cnode.cell.cell_name))
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
        node.flip = flipflop
        node.array = array
        instance.cnode = node
        instance.cnode.parent_cnode = parent_cnode
        instance.pin['org'] = node
        instance.org = node

        connect_cnode(parent_cnode, instance.cnode, None)

        # all other nodes:
        for key, pin in self.pin.items():
            if True: #the name 'org' maybe outside the cnode. #key is not 'org':
                node = Node()
                node.xs = pin.xs
                node.width = pin.width
                node.type = pin.type
                node.show = pin.show
                node.remark = pin.remark
                node.name = '{}_{}_{}'.format(self.cell_name, key, Cell.I)
                node.cnode = instance.cnode
                instance.pin[key] = node
                # reconstruct connectivity
                # 1. org-to-node connectivity
                x, y, a = pin.pointer.get_xya()
                if flipflop:
                    y, a = -y, -a
                if flop:
                    a += 180
                connect_geo(instance.cnode, node, (x, y, a), solve=False)
                # 2. copy pointer properties (stored position xya is irrelavant)
                node.pointer = pin.pointer.copy()

        instance.pinin = instance.pin[self.default_in]
        instance.pinout = instance.pin[self.default_out]
        return instance


    def put(self, *args, flip=False, flop=False, array=None, **kwargs):
        """Instantiate a Cell object and put it in the active cell.

        Args:
            *args (float | Node | Instance): a set of max. 5 unnamed parameters
                interpreted as a connection between two pins.
                (default = (), which connects to the current pin)
            flip (bool): mirror the instance in the vector line (default = False)
                flip is preformed after any translation and/or rotation.
                Hint: for connecting waveguides you only need flip (not flop).
            flop (bool): mirror the instance in the line perpendicular to the vector
                (default = False)
                flop is performed after any translation and/or rotation.
                A flop is the same as rot(180) + a flip
                Hint: for mask assembly or connecting bbox pins you may find both
                flip and flop useful.
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
                nd.example_cell().put() # -> put('a0', cp)
                nd.export_plt()

            2. connect pin 'b0' of cell "example_cell()" to current pin (cp)::

                import nazca as nd
                nd.example_cell().put('b0') # -> put('b0', cp)
                nd.export_plt()

            3. connec default pin of cell "example_cell()" to its parent cell at org + (10, 20, 30)::

                import nazca as nd
                nd.example_cell().put(10, 20, 30) # -> put('a0').move(10, 20, 30)
                nd.export_plt()

            4. connect default pin of cell "example_cell()" to instance C at pin 'b0'::

                import nazca as nd
                C = nd.example_cell().put(0)
                nd.example.put(C.pin['b0']) # -> put('a0', C.pin['b0'])
                nd.export_plt()

            5. connect pin 'b0' of cell "example_cell()" to instance C at pin 'b0'::

                import nazca as nd
                C = nd.example_cell().put(0)
                nd.example_cell().put('b0', C.pin['b0']) # -> put('b0', C.pin['b0'])
                nd.export_plt()

            6. connect pin 'b0' of cell "example_cell()" to instance C at its default out pin, i.e.
               refering to instance C only without pin attribute is interpreted as
               connecting to the default output pin of C.::

                import nazca as nd
                C = nd.example_cell().put(0)
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

        active_cell = cfg.cells[-1]
        parent_cnode = active_cell.cnode
        if self.cnode is parent_cnode:
            print('ERROR: can not connect cell \'{}\' to itself.'.\
                format(active_cell.cell_name))

        instance = self._copy_cell2instance(parent_cnode, flip, flop, array)
        nodeI, nodeC, T = self.parse_instance_connection(C1, C2, C3, C4, C5,
            instance=instance)
        connect_geo(nodeC, nodeI, T)

        if cfg.solve_direct:
    #        print("\nI: '{}'".format(instance.cnode.cell.cell_name))
    #        print("nodeI.pointer:", nodeI.pointer)
    #        print("cnodeI:", nodeI.nb_geo[0][0])
    #        print("cnodeI:", nodeI.nb_geo[0][0])
            mat = nodeI.nb_geo[0][1]
            #postion cnodeI wrt parent cnode over the connected instance node nodeI:
            instance.cnode.pointer.set_mat(nodeI.pointer.trans(mat))
            for pp in instance.cnode.nb_geo:
                pp[0].pointer.set_mat(instance.cnode.pointer.trans(pp[1]))
    #            print(pp[0].name, pp[0].pointer)

            if False:
                for name, P in instance.pin.items():
                    ni = ''
                    if P == nodeI:
                        ni = "*"
                    print(" {}{} nb:{}".format(ni, name, len(P.nb_geo)))
                    for pp in P.nb_geo:
                        if pp[0] == instance.cnode:
                            print('    cnodeI', pp[0].name)
                        elif pp[0] == nodeI:
                            print('    nodeI', pp[0].name)
                        else:
                            print('    other', pp[0].name)

        #update cp:
        p = kwargs.get('newcp', None)
        if p is None:
            p = kwargs.get('cp', None) #backward compatible
        if p is None:
            cfg.cp = instance.pin[self.default_out]
        else:
            cfg.cp = instance.pin[p]

        #for ID in trace.trace_id_list:
        if trace.trace_id_list:
            trace.trace_append(instance)

        return instance


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


    def _poly_iter(self):
        """Generator to iterate over all polylines and polygons.

        Yields:
            tuple: (org, points)
        """
        for org, polygon in self.polygons:
            yield org, polygon.points
        for org, polyline in self.polylines:
            yield org, polyline.points


    def _solve(self, bbox=True):
        """Solve the graph of the cell and calculate the bounding box.

        Args:
            bbox (bool): Calculate the bounding box if True.


        Retuurns:
            None
        """
        if not cfg.solve_direct:
            Netlist().solvecell(self)

        if bbox:
            limit = 1e8
            self.bbox = (limit, limit, -limit, -limit)
            flip = False
            s = 1
            if flip:
                s = -1

            elements = 0
            for org, points in self._poly_iter():
                elements += 1
                if flip:
                    org = org.copy()
                    org.flip()
                [x0, y0, a0] = org.xya()
                a = s*np.radians(a0)
                for u, v in points:
                    x, y = x0+cos(a)*u-sin(a)*v, y0+s*(sin(a)*u+cos(a)*v)
                    self.bbox = (
                         min(self.bbox[0], x), min(self.bbox[1], y),
                         max(self.bbox[2], x), max(self.bbox[3], y))

            for cnode in self.cnode.cnode_nb_iter():
                elements += 1
                org = cnode.pointer.copy()
                if cnode.flip:
                    org.flip()
                [x0, y0, a0] = org.xya()
                a = s*np.radians(a0)
                if cnode.cell.bbox is None:
                    self.bbox_complete = False
                    continue
                x1, y1, x2, y2 = cnode.cell.bbox
                xbl, xbr, xtl, xtr = cos(a)*x1-sin(a)*y1, cos(a)*x2-sin(a)*y2, cos(a)*x2-sin(a)*y1, cos(a)*x1-sin(a)*y2
                ybl, ybr, ytl, ytr = sin(a)*x1+cos(a)*y1, sin(a)*x2+cos(a)*y2, sin(a)*x2+cos(a)*y1, sin(a)*x1+cos(a)*y2

                xmin = min(xbl, xbr, xtl, xtr)
                ymin = min(ybl, ybr, ytl, ytr)
                xmax = max(xbl, xbr, xtl, xtr)
                ymax = max(ybl, ybr, ytl, ytr)
                if cnode.array is not None:
                    Nx, (x1, y1), Ny, (x2, y2) = cnode.array
                    Xmin = min(0, (Nx-1)*x1/Nx, (Ny-1)*x2/Ny, (Nx-1)*x1/Nx+(Ny-1)*x2/Ny)
                    Ymin = min(0, (Nx-1)*y1/Nx, (Ny-1)*y2/Ny, (Nx-1)*y1/Nx+(Ny-1)*y2/Ny)
                    Xmax = max(0, (Nx-1)*x1/Nx, (Ny-1)*x2/Ny, (Nx-1)*x1/Nx+(Ny-1)*x2/Ny)
                    Ymax = max(0, (Nx-1)*y1/Nx, (Ny-1)*y2/Ny, (Nx-1)*y1/Nx+(Ny-1)*y2/Ny)
                else:
                    Xmin, Ymin, Xmax, Ymax = 0, 0, 0, 0

                self.bbox = (
                    min(self.bbox[0], x0 + xmin + Xmin),
                    min(self.bbox[1], y0 + ymin + Ymin),
                    max(self.bbox[2], x0 + xmax + Xmax),
                    max(self.bbox[3], y0 + ymax + Ymax))

            if elements == 0:
                self.bbox_complete = False
                self.bbox = None

            #print("-- bbox: {}: ({:.3f}, {:.3f}, {:.3f}, {:.3f})".\
            #     format(self.cell_name, *self.bbox ))
        return None


    def close(self):
        """Close the cell.

        Solve the geometry of the nodes in the Cell.
        Set the cp back to the position before opening the cell.
        If default input and/or output pins have not been set yet they will be
        set on 'org', pointing in opposite directions.

        Returns:
            Cell: self
        """
        if self.closed:
            print('Cell \'{}\' already closed.'.format(self.cell_name))

        #add default ports if missing:
        if self.default_in not in self.pin:
            self._put_pin(self.default_in, (0, 0, 180))
            #print("Warning: no default pin '{}' set in '{}', setting it to {}.".\
            #    format(self.default_in, self.cell_name, (0, 0, 180)))
        if self.default_out not in self.pin:
            self._put_pin(self.default_out, (0))
            #print("Warning: no default pin '{}' set in '{}', setting it to {}.".\
            #    format(self.default_out, self.cell_name, (0, 0, 0)))
            #self.put_pin(self.default_out, cfg.cp)

        self.pinin = self.pin[self.default_in]
        self.pinout = self.pin[self.default_out]
        self.bbox_complete = True
        self._solve(bbox=True)
        if self.autobbox:
            if self.bbox is None:
                print("Warning: Can't determine a bounding box for cell '{}',"\
                    " e.g. due to only non-native GDS subcell(s) and/or empty subcell(s)."\
                    " Try option native=True in gds loading to get the complete bbox"\
                    " or set bbox=False to get rid of this message.".\
                    format(self.cell_name))
            else:
                if not self.bbox_complete:
                    print("Warning: Not all instances in cell '{}' have a known bbox,"\
                        " e.g. due to a non-native GDS subcell or an empty subcell."
                        " Therefore, the bbox maybe smaller than the actual structures."\
                        " Try option native=True in gds loading to get the complete bbox"\
                        " or set bbox=False to get rid of this message.".\
                        format(self.cell_name))
                length = self.bbox[2] - self.bbox[0]
                width = self.bbox[3] - self.bbox[1]
                put_boundingbox('org', length, width,
                    move=(self.bbox[0], self.bbox[1]+0.5*width, 0))
        self._solve(bbox=False)
        self.closed = True

        cfg.cp = self.oldcp
        cfg.cells.pop()
        if cfg.cells:
            cfg.self = cfg.cells[-1]
        return self


def diff(node1, node2):
    """Calculate the geometrical difference between two nodes.

    Args:
        node1 (Node): start node
        node2 (Node): end node

    Returns:
        (dx, dy, da): difference vector <node2> - <node2> between their pointers
    """
    if not cfg.solve_direct:
        cfg.cells[-1]._solve()
    ptr1 = node1.pointer.copy()
    ptr2 = node2.pointer.copy()
    ptr2.multiply_ptr(ptr1.inv())
    if not cfg.solve_direct:
        cfg.cells[-1].closeflag = False
    #print('nd.diff:', ptr2.get_xya())
    dx, dy, da = ptr2.get_xya()
    if da == 360:
        da = 0
    return dx, dy, da


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
    def __init__(self, name=None, xs=None, width=None, type=None, pin=None,
            show=False, remark=None, chain=1):
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
        self.remark = remark

        if pin is not None:
            self.xs = pin.xs
            self.width = pin.width
        if self.xs is None:
            self.xs = xs
        if self.width is None:
            self.width = width

    def __repr__(self):
        return "<Pin() object, name='{}', xs='{}', width={}, type={}, chain={}, show={}, remark={}>".\
            format(self.name, self.xs, self.width, self.type, self.chain,
                self.show, self.remark)

    def put(self, *args):
        """Put a Pin object in the layout."""
        if self.pin is not None:
            cfg.cp = self.pin
        return cfg.cells[-1]._put_pin(self.name, connect=args, xs=self.xs,
            width=self.width, type=self.type, chain=self.chain, show=self.show,
            remark=self.remark)


class Polygon():
    """Polygon class."""
    def __init__(self, points=None, layer=None):
        """Construct a Polygon object.

        Args:
            pinname (str): name of the pin (obsolete)
            layer (int): layer number to put the Polygon in
            points (list): list of points [(x1, y1), (x2, y2), ...]

        Returns:
            None
        """
        self.layer = get_layer(layer)
        self.points = points
        xy = list(zip(*points))
        self.bbox = (min(xy[0]), min(xy[1]), max(xy[0]), max(xy[1]))

    def __repr__(self):
        try:
            size = len(self.points)
        except:
            size = 0;
        return "<Polygon() object, layer={}, points={}, bbox={}>".\
            format(self.layer, size, self.bbox)

    def put(self, *args):
        """Put a Polygon object in the layout."""
        return cfg.cells[-1]._put_polygon(connect=args, polygon=self)


class Polyline():
    """Polyline (path) class."""

    def __init__(self, points=None, layer=None, width=None, pathtype=0):
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
        self.layer = get_layer(layer)
        if width is not None:
            self.width = width
        else:
            self.width = 0.20
        self.pathtype = pathtype
        self.points = points
        xy = list(zip(*points))
        self.bbox = (min(xy[0]), min(xy[1]), max(xy[0]), max(xy[1]))

    def __repr__(self):
        try:
            size = len(self.points)
        except:
            size = 0
        return "<Polyline() object, layer={}, width={}, pathtype={}, points={}, bbox={}>".\
            format(self.layer, self.width, self.pathtype, size, self.bbox)

    def put(self, *args):
        """Put a Polyline object in the layout."""
        return cfg.cells[-1]._put_polyline(connect=args, polyline=self)


class Annotation():
    """Annotation class."""
    def __init__(self, layer=None, text=''):
        """Construct an annotation object.

        Args:
            pinname (str): name of the pin (obsolete)
            layer (int): layer number to put the Polyline in
            text (str): annotation text

        Returns:
            None
        """
        self.layer = get_layer(layer)
        self.text = text

    def __repr__(self):
        if len(self.text) >= 10:
            text = self.text[:7]+'...'
        else:
            text = self.text
        return "<Annotation() object, layer={}, text='{}'>".format(
            self.layer, text)

    def put(self, *args):
        """Put an Annotation object in the layout."""
        return cfg.cells[-1]._put_annotation(connect=args, annotation=self)


def _scan_branch(strm, cellname, celllist=None, level=0):
    """Generator for cell names in a branch. Deepest first."""

    cell_rec = strm.cells[cellname]
    snames = cell_rec.snames
    #print("cell:{}, snames:{}".format(cellname, snames))
    if snames is not None:
        level += 1
        for sname in snames:
            yield from _scan_branch(strm, sname, level=level)
            #print("{}{}".format("  "*level, sname))
        level -= 1
    yield(cellname)


def _gds2native(strm, topcellname=None, cellmap=None, bbox=False):
    """Get a list of all filenames.

    If no cellname is provided it will select the topcell if there is only
    one topcell.

    Args:
        filename (str): gds filename
        cellname (str): cellname (and under) to process

        native: (bool):  create native Nazca cells from GDS

    Returns:
        list of str: list of cellnames in <filename> under cell name <cellname>
    """
    NC = {}
    elements_unknown = set()
    cells_visited = set()
    branch_iter = _scan_branch(strm, topcellname)
    for cellname in branch_iter:
        if cellname in cells_visited:
            continue
        cells_visited.add(cellname)
        cellrecord = strm.cells[cellname]
        with Cell(cellname) as C:
           for elem in cellrecord.elements:
                if elem.etype == gb.GDS_record.BOUNDARY:
                    LD, XY = elem.polygon
                    XY = [(x/1000, y/1000) for x, y, in zip(*[iter(XY)]*2)]
                    Polygon(points=XY, layer=LD).put(0)
                elif elem.etype == gb.GDS_record.PATH:
                    LD, XY = elem.polyline
                    XY = [(x/1000, y/1000) for x, y, in zip(*[iter(XY)]*2)]
                    Polyline(points=XY, layer=LD).put(0)
                elif elem.etype == gb.GDS_record.TEXT:
                    LD, XY, TEXT = elem.annotation
                    XY = [(x/1000, y/1000) for x, y, in zip(*[iter(XY)]*2)]
                    Annotation(layer=LD, text=TEXT).put(XY)
                elif elem.etype == gb.GDS_record.BOX:
                    elements_unknown.add('BOX')
                elif elem.etype == gb.GDS_record.SREF:
                    name, trans, mag, angle, XY = elem.instance
                    if trans == 1:
                        flip = True
                    else:
                        flip = False
                    NC[name].put(XY[0]/1000, XY[1]/1000, angle, flip=flip, scale=mag)
                elif elem.etype == gb.GDS_record.AREF:
                    name, trans, mag, angle, col, row, XY = elem.array
                    if trans == 1:
                        flip = True
                    else:
                        flip = False
                    x, y = XY[0]/1000, XY[1]/1000
                    x1, y1 = XY[2]/1000, XY[3]/1000
                    x2, y2 = XY[4]/1000, XY[5]/1000
                    dx1, dy1 = (x1-x)/col, (y1-y)/col
                    dx2, dy2 = (x2-x)/row, (y2-y)/row
                    NC[name].put(x, y, angle,
                        array=(col, (dx1, dy1), row, (dx2, dy2)),
                        flip=flip, scale=mag)
                else:
                    elements_unknown.add(elem.etype)
           if (cellname == topcellname) and bbox:
               C.autobbox = True
        NC[cellname] = C

    if elements_unknown:
        print("Warning: gds2nazca. Unknown elements in cell (branch) '{}': {}".\
            format(cellname, elements_unknown))
    return C


def gds_filter(filename, cellmap=None, layermap=None):
    """Filter layers and cells from gds and export a copy to file."""
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


def load_gds(filename, cellname=None, newcellname=None, layermap=None,
        cellmap=None, scale=1.0, prefix='', instantiate=False, native=True,
        bbox=False, connect=None):
    """Load a GDS cell <celllname> (and its instances) from <filename> into a Nazca cell.

    This method checks for cellname clashes, i.e. a cell(tree)
    loaded into the layout can not contain a cell name that already exists
    in the layout. There are three ways to avoid name clashes in case they occur:

    1. Set a <newcellname>: changes only the name of the top cell;
    Note <newcellname> is ignored if a <cellmap> is provided,
    i.e. add any a topcell renaming to the cellmap in this case.

    2. Provide a <cellmap>: maps original cell names one-by-one to new cell names.
    Cells names omitted will not be renamed.

    3. Provide a cell name <prefix>: applies to all cells in a branch,
    except those in <newcellname> or <cellmap>.

    Note, if you need to create a building block from GDS with pins and stubs,
    and you have the pin and xsection info available already (for example in a file)
    you may prefer to use method 'load_gdsBB' instead.

    Args:
        filename (str): filename of gds to load
        cellname (str): cellname in filename to load
        newcellname (str): new name for <cellname>
        layermap (dict): mapping {old_layernumber: new_layernumber}
        cellmap (dict): mapping {old_cellname: new_cellname}
        scale (float): scaling factor if the cell (not yet implemented)
        prefix (str): optional string to avoid name clashes (default = '')
        instantiate (bool): instantiate the GDS (default = False)
        native (bool): create native Nazca cells from GDS
        connect (float, float, float): move origin (x, y, a) when placing (when not Native)

    Returns:
        Cell: Nazca Cell containing the loaded gds cell(tree)

    Example::

        load_gds(filename='path_to_file.gds', cellname='name')
    """

    #TODO: The load_GDS method loads the gds in a stream object to check the cellnames
    #    After that, the stream is disregarded. Only the filename, cell reference
    #    are stored. At mask export time the GDS is loaded again. For large GDS
    #    it may be wiser to store the stream in memory for reuse at export time.
    #print('load_gds: {} -- {}'.format(filename, cellname))
    global num
    num += 1

    strm = gdsimp.GDSII_stream(filename)

    #find cellname of branch to mount in design tree
    if cellname is None:
        topcell = strm.topcell()
        if len(topcell) == 1 and topcell is not None:
            topcellname = topcell.pop()
        else:
            raise Exception("load_gds: No cellname provided to load from file '{}'."\
                " Trying the topcell, but multiple topcells exist."\
                " Please, specify the cellname to load."\
                " Available are:\n{}".\
                format(filename, '\n'.join(list(strm.cells))))
    else:
        topcellname = cellname

    #Build cellmap:
    if cellmap is None:
        cellmap = dict()
        if newcellname is not None:
            cellmap[topcellname] = newcellname

    #check if cellnames in loaded GDS cell tree clash with existing cell names.
    names = strm.cell_branch(topcellname)
    for name in sorted(names):
        if name not in cellmap:
            cellmap[name] = '{}{}'.format(prefix, name)
        if (cellmap[name] in cfg.cellnames.keys()) and not native: #native solves theh conflicts itself
            if cellmap[name] is not name:
                rename = "is renamed to '{0}' which ".format(cellmap[name])
            else:
                rename = ''
            raise Exception(
"""
Error in load_gds: cell name overlap: '{0}'.

Cell name '{3}' in file '{1}' {4}is already in use in the design.

  Five suggestions to solve this issue by giving a <uniquename>:
    In 'load_gds':
    1. Rename the cell (topcell only): load_gds(..., newcellname=<uniquename>)
    2. Apply a cellmap (0 or more cells in a branch): load_gds(..., cellmap=<cellmap_dict>)
    3. Apply a prefix (to all cells in a branch): load_gds(..., prefix='new'). Present prefix is '{2}'
    Other:
    4. Rename cell '{3}' in gds file '{1}'
    5. Rename existing cell name '{0}' in the design
""".format(cellmap[name], filename, prefix, name, rename))

        if not native: # or Cell will see it as existing.
            cfg.cellnames[cellmap[name]] = 'loaded_gds_cell'


    strm = gdsimp.GDSII_stream(filename, layermap=layermap, cellmap=cellmap)

    #create cell:
    if native:
        if scale != 1:
            print("Warning: Scaling not applied for load_gds with native=True."\
                " scale={} will be ignored in cell '{}'".format(scale, topcellname))
        maskcell = _gds2native(strm=strm, topcellname=cellmap[topcellname],
            cellmap=cellmap, bbox=bbox)
    else:
        with Cell('load_gds', celltype='mask', instantiate=instantiate) as maskcell:
            maskcell._put_gds(
                connect=connect,
                filename=filename,
                cellname=topcellname,
                newcellname=cellmap[topcellname],
                layermap=layermap,
                cellmap=cellmap,
                scale=scale,
                strm=strm)
            maskcell.autobbox = bbox

    return maskcell



def load_gds_raw(filename, cellname=None, newcellname=None, layermap=None,
        cellmap=None, scale=1.0, prefix='', instantiate=False,
        bbox=False, connect=None, **kwargs):
    """Load GDS and always force native=False

    load_gds_raw(...) is the same as load_gds(..., native=False).
    See load_gds() for a detailed method description.

    Returns:
        Cell: Nazca Cell containing the loaded gds cell(tree)
    """
    return load_gds_raw(filename=filename,
        cellname=cellname, newcellname=newcellname,
        layermap=layermap, cellmap=cellmap,
        scale=scale, prefix=prefix,
        instantiate=instantiate,
        native=False,
        bbox=bbox, connect=connect)


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
        """Calculate the geometrical positions of nodes in a cell via the netlist.

        This includes all nodes in the cell and all node in instantiated cells
        one level deep. Hence a node in scope adheres to one of the following
        conditions:

        - the node resides in same cell as the starting node (parent cnode0).
        - the node resides in an instance of the cell of the starting node.

        A starting <node> has to be provided from where to solve.
        A good starting point in each cell is the cnode with its pointer
        position set at (x, y, a) = (0, 0, 0).

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
                if (nextnode.cnode is cnode0):
                    nextnode.pointer.set_mat(node.pointer.trans(mat))
                    self.solvecell_core(nextnode, cnode0)
                elif (nextnode.cnode.parent_cnode is cnode0):
                    nextnode.pointer.set_mat(node.pointer.trans(mat))
                    self.solvecell_core(nextnode, cnode0)
                else:
                    print("...Warning: skipping node in cell '{}': it is not in scope of cell '{}'.".format(
                        nextnode.cnode.cell.cell_name, cnode0.cell.cell_name))
            #TODO: else: check geo consistency


    def solvecell(self, cell):
        """Solve the cell with. cnode position will be at (x, y, a) = (0, 0, 0).

        Args:
            cell (Cell): cell to solve

        Returns:
            None
        """
        if cfg.solve_direct:
            #print("SOLVE-LIVE_skipping solvecell")
            return None
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
        for node, pgon in cnode.cell.polygons:
            if infolevel > 1:
                print('{}polygon: ML={}, xya={}, bbox={}'.format('  '*infolevel,
                    layer, node.pointer.get_xya(), bbox))
            yield node.pointer.copy(), pgon


    def polyline_iter(self, cnode, infolevel=0):
        """Iterator for all polylines (paths) in a cell (cnode).

        Yields:
            pointer, width, layer, path type points:
                Next Polyline and info to generate it

        """
        for node, pline in cnode.cell.polylines:
            if infolevel > 1:
                print('{}polyline: ML={}, xya={}, bbox={}'.format('  '*infolevel,
                    layer, node.pointer.get_xya(), bbox))
            yield node.pointer.copy(), pline


    def annotation_iter(self, cnode, infolevel=0):
        """Iterator for all annotations in a cell (cnode).

        Yields:
            pointer, layer, text:
                Next Annotation and info to generate it
        """
        for node, anno in cnode.cell.annotations:
            if infolevel > 1:
                print('{}annotation: ML={}, xya={}'.format('  '*infolevel,
                    layer, node.pointer.get_xya()))
            yield node.pointer.copy(), anno


    def gdsfile_iter(self, cnode, level=0, infolevel=0):
        """Iterator for all GDS files instantiated in a cell (cnode).

        Yields:
            pointer, file, cell, newcell, layermap, cellmap, scale:
                Next gds file and info to generate it
        """
        gdsfiles = cnode.cell.gdsfiles
        for gdsinfo in gdsfiles:
            pin, file, cell, newcell, layermap, cellmap, scale, strm = gdsinfo
            out = pin.pointer.copy(), file, cell, newcell, layermap, cellmap,\
                scale, strm
            if infolevel > 1:
                print('{}gds-cell: cell={}, xya={}'.format('  '*infolevel,
                    newcell, cnode.pointer.get_xya()))
            yield out


    def celltree_iter(self, cnode, level=0, position=None, flat=False, infolevel=0):
        """Iterate over the cell hierarchy from the topcell at 'cnode', deep-first.

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
        cnodeC = cnode.cell.cnode #get cell cnode for cell or instance
        self.cnodes_visited.add(cnodeC)

        if position is None: #coordinate
            position = Pointer()
        if infolevel > 1:
            print('{}CT: cnode={}, name={}, level={}'.format(
                '  '*infolevel, cnodeC, cnodeC.cell.cell_name, level))
        try:
            array = cnode.array
        except:
            array = None
        if (array is not None) and \
                ((cnodeC.cell.instantiate is False) or (flat is True)):
            nx, (dx1, dy1), ny, (dx2, dy2) = array
            dx1 = dx1/nx
            dy1 = dy1/nx
            dx2 = dx2/ny
            dy2 = dy2/ny
            for xi in range(nx):
                for yi in range(ny):
                    org2 = position.copy()
                    if array is not None:
                        a = org2.a
                        org2.set_a(0)
                        org2.move(xi*dx1+yi*dx2, xi*dy1+yi*dy2)
                        org2.set_a(a)
                    yield cnodeC, level, org2, flip
                    cnode_iter = cnodeC.cnode_nb_iter()
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
                                flat,
                                infolevel)

        else:
            yield cnodeC, level, position, flip
            cnode_iter = cnodeC.cnode_nb_iter()
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
                        flat,
                        infolevel)


#-----------------------------------------------------------------------------
if __name__ == '__main__':
    pass