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
# Test replacement of gds cells from one file with those from another.
# This will be used for replacing black-box building blocks with their real
# implementation.
#
# (c) 2016-2017 Katarzyna Lawniczuk, Ronald Broeke
#==============================================================================

"""
Module for replacing cells, static and parametric cells.
"""


import time
from collections import Counter
import pprint

import nazca as nd



def __findBlackCells(blackBasename, cells=[]):
    """
    Find all "black" cells having prefix <blackBasename> in given list of <cells>.

    Args:
        blackBasename (str): string match first part (basename) of black cells
        cells (list): list of cell names

    Returns:
        list of cell names in <cells> starting with <blackprefix>
    """

    matchingcells = []
    for name in cells:
        if name.startswith(blackBasename):
            matchingcells.append(name)
    return matchingcells


def __readPcellParams(allCellsOf1type, gdsin, infolevel=0):
    """
    Return a dictionary with all Pcell parameter values per cell.

    Args:
        allCellsOf1type (str): pcell base name (without params string suffix)
        gdsin (str): gds input file name
        infolevel (int): amount of debug info that is printed, 0 is minimum.

    Returns:
        dict: cells and their parameters dict {cellname: {var: value, var: value})

    """

    if infolevel > 1:
        print("Read params from '{}'".format(allCellsOf1type))
    allPcellParams = {} # dictionary of {cell: PcellStrParams}
    for cellname in allCellsOf1type:
        PcellStrParams = {}
        PcellValueParams = {}

        for lay, pos, text in nd.get_cell_annotation(gdsin.cells[cellname]):
            if 'Parameter' in text:
                PcellStrParams = nd.string_to_parameters(text)

        for varname, valuestr in PcellStrParams.items():
            try:
                value = int(valuestr.split(' ')[0])
            except:
                value = float(valuestr.split(' ')[0])
            PcellValueParams[varname] = value

        allPcellParams[cellname] = PcellValueParams
        #if infolevel > 1:
        #    print('{}\n{}\n{}', cellname, PcellStrParams, PcellValueParams)

    return allPcellParams


def __createWhitePcellLib(gdsin, whitelibrary=None, black2whiteMap=None,
        prefixb='black_', prefixw='wb_', infolevel=0):
    """
    Create a Pcell library based in cells <gdsin>.

    Args:
        gdsin (str): file name of input GDS.
        gdsout (str): optional file name of output GDS.
        whitelibrary (str): optional file name of generated white cell GDS library.
        black2whitemap (dict): mapping black cell names to white cell functions {name: function}.
        prefixb: blackbox prefix (default = 'black_')
        prefixw (str): whitebox prefix  (default = 'wb_')
        infolevel (int): amount of debug info printed

    Returns:
        str: filename of white gds Pcell library.
    """

    if whitelibrary is None:
        timestr = time.strftime("%Y-%m-%d")
        whitelibrary = "{}_white_lib_{}.gds".format(gdsin[:-4], timestr)

    # Using the black2whiteMap dictionary:
    # - create a list of all black cellnames in the map,
    # - extract their Pcell parameters,
    # - create a matching list of newly generated white cells using the params.
    # Note the white cell is (should be) generated with a black cell inside.
    # Note that the replaced white cells have the same name for the gds ref
    gdsinstream = nd.GDSII_stream(gdsin)
    allInputCellNames = gdsinstream.cells.keys()
    allBlackCellNames = []
    allWhiteCells = []
    for blackBasename, whiteFunction in black2whiteMap.items():
        blackCellGroup = __findBlackCells(blackBasename, allInputCellNames)
        allBlackCellNames += blackCellGroup

        cellsParams = __readPcellParams(blackCellGroup, gdsinstream,
            infolevel=infolevel)
        for cellname, parameters in cellsParams.items():
            try:
                whiteCell = whiteFunction(**parameters)
                allWhiteCells.append(whiteCell)
            except:
                print("Error: Can not generate whitebox {} with params {}".\
                    format(whiteFunction, parameters))
                #TODO: add empty/error cell to keep black and white list in sync
        if infolevel > 1:
            print('\nblackBasename = ', blackBasename)
            print('\nwhiteCellparams: ', allWhiteCells)

    if infolevel > 1:
        print('\nallBlackCells: {}'.format(allBlackCellNames))
        print('\nallWhiteCells: {}'.format(allWhiteCells))
        print('\n')

    #generate output cell names for all black and white cells using prefixes
    allBlackCellNamesOut = []
    allWhiteCellNames = []
    for name in allBlackCellNames:
        allBlackCellNamesOut.append(prefixb+name)
        allWhiteCellNames.append(prefixw+name)

    #map input black cellname to output black cellname (for gds debugging).
    blackmap = {old:new for old, new in zip(allBlackCellNames, allBlackCellNamesOut)}
    #map white cellname to black cellname
    whitemap = {white:black for white, black in zip(allWhiteCellNames, allBlackCellNames)}
    if infolevel > 0:
        print('mapping cell->black:')
        pprint.pprint(blackmap)
        print('mapping white->black cells:')
        pprint.pprint(whitemap)

    #Export all white cells into the <whitelibrary> gds file
    #TODO: use whitemap in export:
    nd.export_gds(allWhiteCells, filename=whitelibrary)
    # rename black cells:
    g = nd.GDSII_stream(whitelibrary, cellmap=blackmap)
    g.GDSII_write(whitelibrary)
    # rename white cells into original black cell names:
    g = nd.GDSII_stream(whitelibrary, cellmap=whitemap)
    g.GDSII_write(whitelibrary)

    #white to black cell name map:
    w2b = dict(zip(blackmap.keys(), blackmap.keys()))
    return whitelibrary, w2b


def replaceCells(gdsin, gdsout=None, PcellFunctionMap=None,
        ScellMapping=(None, None), infolevel=0):
    """Replace black with white cells in <gdsin> and export result to <gdsout>.

    Replacement is Pcell (parametric cell) or Scell (static cell) based.
    Only apply one mapping per function call, Pcell or Scell, to keep track of
    the mapping order. Sequential calls are possible.

    * Pcells: Needs a mapping of the Pcell basename to a cell function.
    * Scells: Needs a gds library <ScellFile> and a map <ScellMap>
        * ScellFile (str): filename of gds library
        * ScellMap (dict): {black_cellname: white_cellname}.

    Args:
        gdsin (str): input filename
        gdsout (str): optional output filename
        PcellMap (dict): black2white dictionary {black_cellbasename: white_Pcell_function}
        ScellMapping ((str, dict)): (<ScellFile>, <ScellMap>) for black2white mapping
        infolevel (int): amount of debug info printed

    Returns:
        str: gds output filename of file with replaced cells
    """

    if gdsout is None:
        timestr = time.strftime("%Y-%m-%d")
        gdsout = "{}_white_{}.gds".format(gdsin[:-4], timestr)

    if PcellFunctionMap is not None:
        whitePcellLibname, PcellMap = __createWhitePcellLib(gdsin,
            black2whiteMap=PcellFunctionMap, infolevel=infolevel)
        whitePcellstream = nd.GDSII_stream(whitePcellLibname)
        whitePcells = Counter(whitePcellstream.cells.keys())
        print("...created white Pcell gds library '{}'".format(whitePcellLibname))
        print("...cellmap: '{}'")
        pprint.pprint(PcellMap)
    else:
        whitePcells = Counter([])

    ScellFile, ScellMap = ScellMapping
    if ScellMap is None:
        ScellMap = []

    if ScellFile is not None:
        whiteScellstream = nd.GDSII_stream(ScellFile)
        whiteScells = Counter(whiteScellstream.cells.keys())
    else:
        whiteScells = Counter([])


    gdsinstream = nd.GDSII_stream(gdsin)
    gdsincells = Counter(gdsinstream.cells.keys())

    if len(whiteScells) > 0 and len(whitePcells) > 0:
        print("Use only Pcells or Scells in one call.")
        return None

    if len(whiteScells) > 0:
        whiteCells = whiteScells
        whiteCellstream = whiteScellstream
        cellMap = ScellMap
    elif len(whitePcells) > 0:
        whiteCells = whitePcells
        whiteCellstream = whitePcellstream
        cellMap = PcellMap


#==============================================================================
# Find which cells to replace and where needed rename
#==============================================================================
    replace = {}
    noreplace = {}
    for Bname, Wname in cellMap.items():
        if Bname in gdsincells:
            if Wname in whiteCells:
                #print("{} -> {}".format(Bname, Wname))
                replace[Bname] = Wname
            else:
                #print("{} -> no:{}".format(Bname, Wname))
                noreplace[Bname] = Wname
        else:
            #print("no:{}".format(Bname))
            noreplace[Bname] = Wname
    #list unmapped cells:
    stay = list(gdsincells - Counter(replace.keys()))
    #list any unmapped gdsin names in conflict with new cell names after map?
    stayNOTOK = list( (Counter(stay) & Counter(replace.values())).elements())
    stayOK = list(Counter(stay)-Counter(stayNOTOK))

    print('Validated mapping(s):')
    pprint.pprint(replace)
    if len(noreplace) > 0:
        print('Invalid mapping(s) removed:')
        pprint.pprint(noreplace)
    #print('gdsin cell(s) with name conflict:')
    #pprint.pprint(stayNOTOK)
    #print('untouched gdsin cell(s):')
    #pprint.pprint(stayOK)

    #create cellmap for gdsinstream and read gds again.
    #TODO: in mem
    cellmap = dict(replace)
    for cell in stayNOTOK:
        newname = cell+'$'
        cellmap[cell] = newname
        stayOK.append(newname)
    #print("cellmap:\n", cellmap)

    #print('untouched gdsin cells including renamed ones:')
    #pprint.pprint(stayOK)

    del(gdsinstream)
    gdsinstream = nd.GDSII_stream(gdsin, cellmap=cellmap)
    #print("topcells:")
    #pprint.pprint(gdsinstream.topcell())

#==============================================================================
# Replace cells
#==============================================================================
    #print('\nGenerate gdsout:')
    with open(gdsout, 'wb') as f:
        f.write(b''.join(rec.stream for rec in gdsinstream.header))

        #Walk tree for each top cell rather then linear loop.
        def gdsin_iter(name, level=0):
            """Iterator over cell <name>."""
            nonlocal stayOK
            if name in stayOK:
                yield name, 1
            else:
                yield name, 2
                return None
            for sub in gdsinstream.cells[name].snames:
                yield from gdsin_iter(sub, level+1)

        for top in gdsinstream.topcell():
            for cellname, source in gdsin_iter(top):
                if source == 1:
                    f.write(gdsinstream.cells[cellname].stream)
                elif source == 2:
                    f.write(whiteCellstream.GDSII_stream_cell(cellname))

        #for cellname in whitePcells:
        #    f.write(whitePcellstream.GDSII_stream_cell(cellname))

       # for Bname, Wname in replace.items():
       #     f.write(whiteScellstream.GDSII_stream_cell(Wname))

       # for cellname in stayOK:# gdsinstream.cells.keys():
       #     f.write(gdsinstream.cells[cellname].stream)

        f.write(gdsinstream.footer.stream)

    print("...Wrote white gds file '{}'".format(gdsout))

    return gdsout
