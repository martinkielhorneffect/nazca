#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 13:09:43 2018

@author: rgb
"""

import nazca as nd


nd.load_gds('foo1.gds').put(0, 0)
nd.load_gds('foo3.gds').put(5, 0)

nd.export_gds()