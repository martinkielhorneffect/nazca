#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 09:45:51 2017

@author: rgb
"""

import matplotlib.pyplot as plt
import numpy as np

columns = ('Last', 'High', 'Low', 'Chg.', 'Chg. %', 'Time', 'T?')
rows = ['Gold', 'Silver', 'Copper', 'Aluminum']

data_list = np.random.randint(10,90, size=(len(rows), len(columns)))


fig = plt.figure(1)
fig.subplots_adjust(left=0.2,top=0.8, wspace=1)


cellcolours = np.empty_like(data_list, dtype='object')
cellcolours[1,3] = 'r'
cellcolours[2,4] = '#A2A3F4'

#Table - Main table
ax = plt.subplot2grid((1,3), (0,0), colspan=3, rowspan=1)
ax.table(cellText=data_list,
          rowLabels=rows,
          colLabels=columns, loc="upper center",
          cellColours=cellcolours)

ax.axis("off")

fig.set_size_inches(w=6, h=5)
plt.show()