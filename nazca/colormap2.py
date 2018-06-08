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
Parsing Klayout lyp files into csv via DataFrames: lyp2csv.
Generate a matplotlib colormap from the csv table: csv2colormap
"""

"""
Structure in .lyp XML file:


layer-properties-tabs -> present in case of multiple tabs only
    layer-properties -> represent a tab view in Klayout
        properties
            GROUP
                LAYER2
                LAYER2
            GROUP
                GROUP
                    LAYER2
                    LAYER2
                GROUP
                    LAYER2
                    LAYER2
            LAYER1
            lAYER1

LAYER1:
    <properties>
    <...lyp_attr>
    <name>   -> layer name
    <source> -> layer/datatype

LAYER2:
    <group-members>
        <...lyp_attr>
        <name>   -> layer name
        <source> -> layer/datatype

GROUP:
    <group-members>
        <...lyp_attr>
        <name>   -> group-name
        <source> -> */*
"""

#TODO: add csv2lyp export.
from collections import defaultdict
import pandas as pd
import xml.etree.ElementTree as Tree


lyp_attr   = ['layer', 'datatype', 'source', 'fill-color', 'frame-color',
    'frame-brightness', 'fill-brightness', 'dither-pattern', 'valid',
    'visible', 'transparent', 'width', 'marked', 'animation']
nazca_attr = ['layer', 'datatype', 'name', 'fill_color', 'frame_color',
    'frame_brightness', 'fill_brightness', 'dither_pattern', 'valid',
    'visible', 'transparent', 'width', 'marked', 'animation']


doPrint = False
#==============================================================================
# lyp2csv
#==============================================================================
d = defaultdict(list)
depth=0
def __parse_properties(lev1):
    """Parse lyp tags <properties> and <group-member> levels."""
    global d
    global depth

    depth += 1
    d['depth'].append(depth)
    for lev2 in lev1:
        tag = lev2.tag
        value = lev1.find(tag).text
        if doPrint:
            print('  '*depth + tag + ': ' + str(value))
        if tag == 'group-members':
            __parse_properties(lev2)
        else:
            if tag == 'source':
                value = value[:-2].split('/')
                d['layer'].append(value[0])
                d['datatype'].append(value[1])
            else:
                d[tag].append(value)
    depth -= 1
    return None


def __parse_tab(lev1):
    """Parse a lyp tab section."""
    global d, path

    d = defaultdict(list) #clear global dict for new tab
    tabname = 'nazca'
    for lev2 in lev1: # 'layer-properties'
        tag = lev2.tag
        if doPrint:
            print('  tag1:', tag)
        if tag == 'properties':
            __parse_properties(lev2)
        elif tag == 'name': # group name
            if lev1.find(tag).text is not None:
                tabname = lev1.find(tag).text

    df = pd.DataFrame.from_dict(d)
    df.rename(columns=dict(zip(lyp_attr, nazca_attr)), inplace=True)

    csv_filename = "{}{}{}{}".format(path, 'table_colors_', tabname, '.csv')
    df[['depth']+nazca_attr].to_csv(csv_filename, index=False)
    #if doPrint:
    print("Exported tab '{}' to file '{}'.".format(tabname, csv_filename))
    return tabname, csv_filename


def lyp2csv(lypfile):
    """Convert a Klayout .lyp file into a .csv Nazca color table.

    This function makes it possible to use a Klayout color definition (in xml)
    inside Nazca for Matplotlib output. A .lyp file may contain multiple tabs
    and each tab is exported to a separate .csv file.

    Note that display logic in Klayout and Matplotlib
    are not the same.  Hence, a layout will have a different "style"
    for things like:

        * tranparency and opaqueness
        * visualization of edges
        * stipple

    Args:
        lypfile (str): name of Klayout .lyp file

    Returns:
        dict: dictionary {<tabname>: <csv_filename>}
    """

    #tree = Tree.parse(__clean_lyp(filename))

    tree = Tree.parse(lypfile)
    root = tree.getroot()

    tabs = {}
    # Parse tags <layer-properties-tabs> and <layer-properties>:
    if root.tag == 'layer-properties-tabs':
        for lev1 in root: # <layer-properties>
            if doPrint:
                print('lev1 in root:', lev1.tag)
            tab, filename = __parse_tab(lev1)
            tabs[tab] = filename
    else: # <layer-properties>
        tab, filename = __parse_tab(root)
        tabs[tab] = filename

    return tabs


if __name__ == '__main__':
    from pprint import pprint
    fab = 'smart'
    if fab == 'demo':
        path = 'demofab/'
        file = path + 'demofab_klayout_colors.lyp'
    elif fab =='smart':
        path = '../../pdk/smart/smart/'
        file = path + 'smart-nazca.lyp'
    elif fab =='hhi':
        path = '../../pdk/hhi/hhi/'
        file = path + 'hhi-3.lyp'

    tabs = lyp2csv(file, infolevel=3)
    pprint(tabs)

    pprint(lyp2csv(file))
