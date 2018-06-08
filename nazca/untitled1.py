#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 17:30:52 2017

@author: rgb
"""



import nazca as nd
from nazca.interconnects import Interconnect
from math import degrees

R = 30
ic = Interconnect(width=2.0, radius=R)

solve = ic.sbend_p2p_solve(pin1=(20, 0, 0), pin2=(60, 100, 0), radius=R)
ic.strt(length=solve['length1']).put(40, 0, 0)
ic.bend(angle=degrees(solve['angle'])).put()
ic.strt(length=solve['length2']).put()
ic.bend(angle=-degrees(solve['angle'])).put()
ic.strt(length=solve['length3']).put()


ic.sbend_p2p((0, 0, 0), (40, 100, 0)).put()

nd.export_gds()