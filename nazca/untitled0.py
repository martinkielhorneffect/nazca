#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 22:49:28 2018

@author: rgb
"""


# import the nazca module and call it 'nd'
import nazca as nd
# import the SMART Photonics PDK and call it 'smart'
import smart as smart
# import the interconnects module and call it 'IC'
import nazca.interconnects as IC
# Define a class of interconnects for the deep cross-section:
icd = IC.Interconnect(xs='Deep', width=1.5, radius=100)
# Define a class of interconnects for the metal cross-section:
icm = IC.Interconnect(xs='Metal', width=20, radius=20)

with nd.Cell('MZI') as mzi:
    left_mmi = smart.mmi1x2_dp().put()
    top_eopm = smart.eopm_dc_dp(length=1000, pads=True).put(left_mmi.pin['b0'].move(135,50))
    bot_eopm = smart.eopm_dc_dp(length=1000, pads=True).put('b0',left_mmi.pin['b1'].move(135,-50))
    right_mmi = smart.mmi1x2_dp().put('b1', top_eopm.pin['b0'].move(135,-50))
    icd.sbend_p2p(left_mmi.pin['b0'],
                  top_eopm.pin['a0']).put()
    icd.sbend_p2p(left_mmi.pin['b1'], bot_eopm.pin['b0']).put()
    icd.sbend_p2p(right_mmi.pin['b1'], top_eopm.pin['b0']).put()
    icd.sbend_p2p(right_mmi.pin['b0'], bot_eopm.pin['a0']).put()
    nd.Pin('a0').put(left_mmi.pin['a0'])
    nd.Pin('b0').put(right_mmi.pin['a0'])
    nd.Pin('c0').put(top_eopm.pin['c0'])
    nd.Pin('c1').put(bot_eopm.pin['c0'])

nd.export_plt(mzi)