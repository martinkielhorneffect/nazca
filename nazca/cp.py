#!/usr/bin/env python3
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
# -*- coding: utf-8 -*-
# 2017 (c)  Ronald Broeke


"""
The 'cp' module provides short syntax methods to operate on current pointer cp.
"""

from . import cfg
from .netlist import show_cp

stack = []

#def self(p):
#    return cfg.cp

#here = cfg.cp

def show(d=10, w=1):
    show_cp(d/2.0, w)

def here():
    return cfg.cp

def push():
    stack.append(cfg.cp)
    return cfg.cp

def pop():
    cfg.cp = stack.pop()
    return cfg.cp

def goto(x=0, y=0, a=0):
    cfg.cp.goto(x, y, a)
    return cfg.cp

def move(x=0, y=0, a=0):
    cfg.cp = cfg.cp.move(x, y, a)
    return cfg.cp

def shift(x=0, y=0):
    cfg.cp = cfg.cp.move(x, y, 0)
    return cfg.cp

def move2(x=0, y=0, a=0):
    cfg.cp = cfg.cp.move2(x, y, a)
    return cfg.cp

def rotate(a=0):
    cfg.cp = cfg.cp.rotate(a)
    return cfg.cp
rot = rotate

def skip(x=0):
    cfg.cp = cfg.cp.skip(x)
    return cfg.cp

def offset(x=0):
    cfg.cp = cfg.cp.offset(x)
    return cfg.cp
os = offset

def roll():
    print('cp.roll not implemented.')
    #cfg.cp.roll()
    return cfg.cp

def mirror():
    print('cp.mirror not implemented.')
    #cfg.cp.reflect()
    return cfg.cp

#TODO: this give the coordinates local to the defined cell.
#  Could be better to make is with respect to the active cell scope.
def get_xya():
    return cfg.cp.pointer.get_xya()

def get_xy():
    return cfg.cp.pointer.get_xy()

def get_a():
    return cfg.cp.pointer.get_a()
