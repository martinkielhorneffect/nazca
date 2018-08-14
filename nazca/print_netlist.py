#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 20:46:24 2018

@author: rgb
"""

def foo(topcell):
    """Testing to create and print Nazca tree/netlist with pins to stdout."""
    NL = Netlist()
    cells = NL.celltree_iter2(topcell.cnode, flat=False, infolevel=0)
    for action, params in cells:
        if action == 'new':
            cnode, create, level, trans, flip, levup = params
            if create: # take care of cell hierarchy opening
                print(cnode.cell.cell_name)
                names = [k for k in cnode.cell.pin.keys()]
                print("  ", names)
                for name, node in cnode.cell.pin.items():
                    n = [(E[0].cnode.cell.cell_name, E[0].name) for E in node.nb_geo]
                    print("{} -> {}".format(name, n))