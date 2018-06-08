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
# (c) 2016 Xaveer Leijtens
# (c) 2017 Ronald Broeke

"""Module to read a GDSII file record by record.

Substitute cell names and/or layer names.
Keep cell constructions.

The hierarchy is:

-Cell
    -Element (polygon, polyline, annotation, cell-ref)
        -Records (layer, XY, etc.)

GDSII files in nazca are a subset of the GDSII spec:

gdsii:
    HEADER          \
    BGNLIB          | this is "header" in nazca
    LIBNAME         |
    UNITS           /
    {<structure>}*  > this is "cell" in nazca
    ENDLIB          > this is "footer" in nazca

<structure>:
    BGNSTR
    STRNAME
    {<element>}*
    ENDSTR

<element>:
    {<boundary>|<path>|<sref>|<aref>|<text>|<otherstuff>}
    ENDEL

 GDSII cells are a set of records:
 - start record
 - name record
 - various types of records that can also be references to other cells
 - end record

 All GDSII records are structured as:
 - 16 bit int RECORD_LENGTH (number of bytes in record)
 -  8 bit int RECORD_TYPE
 -  8 bit int DATA_TYPE
 -  the data: (RECORD_LENGTH - 4) / (DATA_TYPE length) elements
"""

import struct
import io
import sys
from collections import OrderedDict, defaultdict
from . import gds_base as gb


elm_open = {
    gb.gds_record.BOUNDARY,
    gb.gds_record.PATH,
    gb.gds_record.SREF,
    gb.gds_record.AREF,
    gb.gds_record.TEXT
}

elm_data = {
    gb.gds_record.LAYER,
    gb.gds_record.DATATYPE,
    gb.gds_record.XY,
    gb.gds_record.SNAME,
    gb.gds_record.ANGLE,
    gb.gds_record.TEXTTYPE,
    gb.gds_record.STRING
}


class GDSII_Error(Exception):
    pass


def DFS(graph):
    """Check if a cell is a DAG.

    GDSII cells form a sorted graph that should not contain loops. This
    routine does a Depth First Search of all cells and returns an
    topologically sorted list.
    """
    def DFS_visit(cell):
        for c in graph[cell]:
            if c in start:
                eprint("Error in GDS\nGraph:", graph)
                raise GDSII_Error("Cyclic reference at cell '{}'.".format(c))
            if c not in parent:
                parent[c] = cell
                start.add(c)
                DFS_visit(c)
                start.remove(c)
                order.append(c)
    parent = {}
    start = set()
    order = []
    cells = graph.keys()
    for c in cells:
        if c not in parent:
            parent[c] = None
            start.add(c)
            DFS_visit(c)
            start.remove(c)
            order.append(c)
#    order.reverse() # Start with top first.
    return order

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def unpack_uint8(byte1):
    return int(struct.unpack('>B', byte1)[0])

def unpack_uint16(byte2): # unsigned
    return int(struct.unpack('>H', byte2)[0])

def unpack_int16(byte2):
    return int(struct.unpack('>h', byte2)[0])

def unpack_int32(byte4):
    return int(struct.unpack('>l', byte4)[0])

# 4-byte real is not used.

def unpack_real8(byte8):
    if byte8 == b'00000000':
        return 0.0
    # value = M / 2**56 * 16**(E - 64)
    #       = M * 16**(E - 64 - 56/4)
    #       = M * 16**(E - 78)
    E, M6, M54, M3210 = struct.unpack(">BBHL",byte8) # 1, 1, 2, 4 bytes
    if E & 0x80:
        sign = -1
    else:
        sign = 1
    E = (E & 0x7f) - 78
    M = M6 * 0x1000000000000 + M54 * 0x100000000 + M3210
    return float(sign * M * (16**E))

def nx_uint8(strm):
    f = io.BytesIO(strm)
    while True: # Can be a stream with a number of 1-byte values
        byte1 = f.read(1)
        if byte1:
            yield unpack_uint8(byte1)
        else:
            break

def nx_int16(strm):
    f = io.BytesIO(strm)
    while True: # Can be a stream with a number of 2-byte values
        byte2 = f.read(2)
        if byte2:
            yield unpack_int16(byte2)
        else:
            break

def nx_int32(strm):
    f = io.BytesIO(strm)
    while True: # Can be a stream with a number of 4-byte values
        byte4 = f.read(4)
        if byte4:
            yield unpack_int32(byte4)
        else:
            break

def nx_real8(strm):
    f = io.BytesIO(strm)
    while True: # Can be a stream with a number of 4-byte values
        byte8 = f.read(8)
        if byte8:
            yield unpack_real8(byte8)
        else:
            break


class GDSII_element:
    """Class for GDS elements.

    A GDSii_element object stores a list of GDSII_records and a number of methods
    to read out the properties (layer and xy position) of each record.

    Elements are polylines, sref, aref, etc.
    Records are layer, XY, etc.
    """

    def __init__(self, records=None):
        """Initialize a GDS element.

        Args:
            records (list of GDSII_record):

        Returns:
            None
        """
        if records is None:
            records = []
        self.records = list(records) # force copy
        etype = None #element type


    def __str__(self):
        return '\n'.join(str(rec) for rec in self.records)


    def addrecord(self, record):
        """Add a GDSII-record to the GDSII_element.

        Args:
            record (GDSII_record):

        Returns:
            None
        """
        self.records.append(record)


    @property
    def etype(self):
        return self.records[0].rtype

    @property
    def annotation(self):
        """Return the properties of an annotation.

        Returns:
            [int, tuple, str]: layer, position (x, y), text
        """
        if self.etype != gb.gds_record.TEXT:
            return None
        lay, XY, text = None, None, None
        for r in self.records[1:]:
            if r.rtype == gb.gds_record.LAYER:
                lay = r.data[0]
            #TODO: add datatype
            if r.rtype == gb.gds_record.XY:
                XY = r.data
            if r.rtype == gb.gds_record.STRING:
                text = r.data
        return [lay, XY, text]


    @property
    def polyline(self):
        """Return the properties of a polyline.

        Returns:
            [int, tuple]: layer, position (x, y)
        """
        if self.etype != gb.gds_record.PATH:
            return None
        lay, XY = None, None
        for r in self.records[1:]:
            if r.rtype == gb.gds_record.LAYER:
                lay = r.data[0]
            #TODO: add datatype
            if r.rtype == gb.gds_record.XY:
                XY = r.data
        return [lay, XY]


    @property
    def polygon(self):
        """Return the properties of a polygon.

        Returns:
            [int, tuple]: layer, position (x, y)
        """
        if self.etype != gb.gds_record.BOUNDARY:
            return None
        lay, XY = None, None
        for r in self.records[1:]:
            if r.rtype == gb.gds_record.LAYER:
                lay = r.data[0]
            #TODO: add datatype
            if r.rtype == gb.gds_record.XY:
                XY = r.data
        return [lay, XY]


    @property
    def stream(self):
        """Create a stream of all records in the GDSII_element.

        Returns:
            bytearray: gds stream
        """
        return b''.join(rec.stream for rec in self.records)


class GDSII_cell:
    """Class for storing GDS cell content.

    A cell contains three attributes that constitute all cell content
        -header stream
        -list of element objects (sref, polygon, polyline, etc)
        -footer stream

    and a set of cell references <snames>
    """

    def __init__(self):
        """Initialize a GDSII_cell.

        Returns:
            None
        """
        self.snames = set()
        self.count_srefs = defaultdict(int) # keep track of number of references to cells
        self.header = []
        self.elements = []
        self.footer = GDSII_record(gb.gds_endstr())


    def __str__(self):
        return str('\n'.join(str(h) for h in self.header) + '\n' +
                '\n'.join(str(elem) for elem in self.elements) + '\n' +
                str(self.footer))


    def references(self, sname):
        """Add a cell reference <sname> to the cell.

        Args:
            sname (str): name of referenced cell

        Returns:
            None
        """
        self.snames.add(sname)
        self.count_srefs[sname] += 1

    def addelem(self, elem):
        """Add an GDSII-element to the cell.

        Args:
            elem (GDSII_element):

        Returns:
            None
        """
        self.elements.append(elem)


    @property
    def name(self):
        """
        Returns:
            str: name of the cell
        """
        return self.header[1].data


    @property
    def stream(self):
        """Create the full stream of the cell.

        Concatenate the list of streams of headers, elements and the footer.

        Returns:
            bytearray: stream of the cell
        """
        return  b''.join(h.stream for h in self.header) + \
                b''.join(e.stream for e in self.elements) + \
                self.footer.stream


class GDSII_record:
    """Class for storing a gdsii record in byte stream format.

    Store GDS records.
    """

    def __init__(self, strm, pos=0):
        """Construct a GDSII_record.

        Args:
            strm (bytearray): stream
            pos (int): start of record in the stream

        Returns:
            None
        """
        self.strm = strm # byte array
        self.pos = pos # Start of record


    def __str__(self):
        if self.rtype in elm_open:
            e = '┌'
            d = '│'
        elif self.rtype in elm_data:
            e = '├'
            d = '│'
        elif self.rtype == gb.gds_record.ENDEL:
            e = '└'
        else:
            e = ''
            d = ''
        s = '{}{}, {}, 4+{} bytes'.format(e,
            gb.gds_record.name[self.rtype],
            gb.gds_value.name[self.dtype],
            self.rlen-4)
        if self.dtype == gb.gds_value.ASCII:
            s += '\n{} "{}"'.format(d, self.data)
        elif self.dtype != gb.gds_value.NODATA:
            s += '\n{} {}'.format(d, self.data)
        return s


    @property
    # GDSII spec defines length as signed 2-byte int, but the sign does not
    # make sense and there are GDSII files with large (> 0x8000) record
    # lengths. By using unsigned here, we can also read those files.
    def rlen(self):
        return unpack_uint16(self.strm[self.pos:self.pos+2])


    @property
    def rtype(self):
        """Record type.

        Returns:
            uint8: record type value
        """
        return unpack_uint8(self.strm[self.pos+2:self.pos+3])


    @property
    def dtype(self):
        """ Data type.

        Returns:
            uint8: data type value
        """
        return unpack_uint8(self.strm[self.pos+3:self.pos+4])


    @property
    def data(self):
        """Convert the gds stream into proper data.

        Returns:
            None
        """
        if self.dtype == gb.gds_value.NODATA:
            return None
        elif self.dtype == gb.gds_value.BITARRAY:
            return 1
        elif self.dtype == gb.gds_value.INT16:
            return list(nx_int16(self.strm[self.pos+4:self.pos+self.rlen]))
        elif self.dtype == gb.gds_value.INT32:
            return list(nx_int32(self.strm[self.pos+4:self.pos+self.rlen]))
        elif self.dtype == gb.gds_value.REAL8:
            return list(nx_real8(self.strm[self.pos+4:self.pos+self.rlen]))
        elif self.dtype == gb.gds_value.ASCII:
            if self.strm[self.pos+self.rlen-1]: # Remove last 0 byte padding
                e = self.pos + self.rlen
            else:
                e = self.pos + self.rlen - 1
            return self.strm[self.pos+4:e].decode('utf-8')
        else:
            raise GDSII_Error('Unknown data type {}'.format(self.dtype))


    @property
    def stream(self):
        """Return the proper bytestream for this record.

        Return:
            bytearray: record
        """
        return self.strm[self.pos:self.pos + self.rlen]


class GDSII_stream:
    """Class to read, modify and write a GDSII stream.

    A stream consists of a header, elements with records and a footer.
    """

    def __init__(self, filename, cellmap=None, layermap=None):
        """Initiliaze a GDSII_stream.

        Args:
            filename (str): file name to read
            cellmap (dict): {celname_old:cellname_new} optional dictionary to map cell names in <filename> onto new
                names in the stream
            layermap (dict): {layer_number_old:layer_number_new} optional dictionary to map layer numbers in <filename> onto
                new layer numbers in the stream

        Returns:
            None
        """
        if cellmap is None:
            cellmap = {}
        if layermap is None:
            layermap = {}
        self.header = []
        self.cells = OrderedDict() # All cells: {cell_name: cell records}
        self.footer = GDSII_record(gb.gds_endlib())

        # Hold the GDSII file structure
        self.snames = set()  # Names of all referenced cells.
        self.graph = dict()  # Holds the DAG structure of the gds.
        self.ref_del = set() # Names of cells referenced from deleted cells.
        self.layer_remove = set() # Names of layers to be removed.
        self.cell_remove = set() # Names of cells to be removed.
        self.layermap = {}
        self.cellmap = {}
        self.filename = filename

        with open(filename, 'rb') as f:
            self.gds_stream = bytearray(f.read())
            self.stream_length = f.tell() # stream length in bytes.
            #eprint("Read '{}'".format(filename))
        for _layi, _layo in layermap.items():
            layi, layo = _layi, _layo
            if isinstance(layi, int):
                layi = (layi, 0)
            if isinstance(layo, int):
                layo = (layo, 0)
            if _layi is None or _layo is None:
                self.layer_remove.add(layi)
            else:
                assert isinstance(layi, tuple)
                assert isinstance(layo, tuple)
                self.layermap[layi] = layo
        for cell in cellmap:
            if cellmap[cell] is None:
                self.cell_remove.add(cell)
            else:
                self.cellmap[cell] = cellmap[cell]
        self.parse()
        for c in self.cells:
            # gather all referenced cells.
            self.snames |= self.cells[c].snames
            # Build the graph
            self.graph[c] = self.cells[c].snames
        # This will check for cyclic references as well
        self.order = DFS(self.graph)


    def count_srefs(self):
        """Print how many times cells have been instantiated."""
        for cellname, cell in sorted(self.cells.items()):
            print("{}:".format(cellname))
            for sref, count in sorted(cell.count_srefs.items()):
                print("   {:4d}x {}".format(count, sref))


    def GDSII_write(self, filename):
        """Write a GDSII stream to file."""
        with open(filename, 'wb') as f:
            f.write(b''.join(rec.stream for rec in self.header))
            if self.cell_remove == set(): # No cells were removed
                for cellname in self.cells:
                    f.write(self.cells[cellname].stream)
            else:
                for cellname in self.topcell():
                    # Write the topcell trees
                    self.GDSII_write_cell_and_under(f, cellname)
            f.write(self.footer.stream)


    def ASCII_write(self, filename=False):
        """Write the GDS in a human readable format."""
        buf = io.StringIO()
        for rec in self.header:
            buf.write(str(rec)+'\n')
        for cellname in self.cells:
            if cellname not in self.cell_remove:
                buf.write(str(self.cells[cellname])+'\n')
        buf.write(str(self.footer)+'\n')
        if filename is False: # Write to stdout
            sys.stdout.write(buf.getvalue())
        elif filename is None: # Don't write, only return the string.
            pass
        else:
            with open(filename, 'w') as f:
                f.write(buf.getvalue())
        return buf.getvalue()


    def GDSII_stream_cell_and_under(self, strm, cellname, init=True):
        """Obtain a GDS stream of cell named <cellname>, and its children.

        Iterative part of method CDSII_stream_cell.

        Args:
            strm (bytearrray): gds stream being build, starting as bytearray()
            cellname (str): cell name of top cell
            init (bool): flag to check is a top cell level (default:True)

        Returns:
            bytearray: gds stream
        """
        global done
        if init:
            done = []
        if cellname not in done:
            try:
                strm.extend(self.cells[cellname].stream)
            except:
                raise Exception("Error: Looking up a non-existing cellname '{}' "\
                    "in file '{}'. Valid cellnames are {}."\
                    .format(cellname, self.filename, list(self.cells.keys())))
            done.append(cellname)

        for subcellname in self.cells[cellname].snames:
            self.GDSII_stream_cell_and_under(strm, subcellname, init=False)
        return strm


    def GDSII_stream_cell(self, cellname):
        """Return the GDSII stream of the cell with name <cellname> and below.

        Args:
             cellname (str): name of the cell to get the gds stream off

        Returns:
            bytearray: gds stream
        """
        return self.GDSII_stream_cell_and_under(bytearray(), cellname,
            init=True)


    def GDSII_write_cell_and_under(self, f, cellname, init=True):
        """Write gds stream of cell named <cellname> to file <f>.

        Iterative part of method GDSII_write_cell.

        Args:
            f (file): file handle of file to write to
            cellname (str): name of the cell to write to file
            init (bool): flag to check if in top cell

        Returns:
            None
        """
        global done
        if init:
            done = []
        if cellname in self.cell_remove:
            return
        if cellname not in done:
            f.write(self.cells[cellname].stream)
            done.append(cellname)
        for subcellname in self.cells[cellname].snames:
            self.GDSII_write_cell_and_under(f, subcellname, False)


    def GDSII_write_cell(self, cellname, filename):
        """Write a cell (and the cells referenced by it).

        Args:
            cellname (str): name of the cell to write
            filename (str): name of the file the gds stream is saved in

        Returns:
            None
        """
        with open(filename, 'wb') as f:
            f.write(b''.join(rec.stream for rec in self.header))
            self.GDSII_write_cell_and_under(f, cellname)
            f.write(self.footer.stream)


    def topcell(self):
        """Return all topcells.

        Returns:
            list of str: list of all cell names that are top cells (are not referenced)
        """
        # This does not change if cells are deleted, which is used to write
        # out the GDSII with deep-deleted cells.
        cells = set(self.cells.keys())
        return cells - self.snames


    def cell_branch(self, cellname, cellnames=None, level=0):
        """Create a list of all cells in branch <cellname>.

        Returns
        """
        if level == 0:
            cellnames = set()
        cellnames.add(cellname)
        for subname in self.cells[cellname].snames:
            self.cell_branch(subname, cellnames, level+1)
        level -= 1
        if level == -1:
            return cellnames


    def print_structure(self, name, level=0):
        """Print the cell tree in ascii format.

        Args:
            name (str): cellname
            level (int): function internal recursive counter (default = 0)

        Returns:
            None
        """
        if level == 0:
            print("□", name)
        else:
            print("  {}└{}".format("│ " * (level-1), name))
        for sub in self.cells[name].snames:
            #TODO: sort alphabetically?
            self.print_structure(sub, level+1)


    # When removing a cell from the file, this is done by removing each
    # reference to that cell, but keeping the cell itself. That will make
    # it a separate topcell. When writing the file, this cell is then
    # not written. Everything below is also discarded, unless also
    # referenced by another cell.
    def gds_cell_iter(self, rec_iter):
        """Return an iterator over all cells.

        Args:
            rec (GDSII_record iterator): record

        Yields:
            GDSII_cell
        """
        for rec in rec_iter: # r is the next record of the record iterator rec.
            if rec.rtype == gb.gds_record.ENDLIB: # End of file.
                break
            cell = GDSII_cell()
            # Read BGNSTR
            if rec.rtype != gb.gds_record.BGNSTR:
                raise GDSII_Error("'BGNSTR' record expected. Got '{}'.".format(
                    gb.gds_record.name[rec.rtype]))
            cell.header.append(rec)
            # Read STRNAME
            rec = next(rec_iter)
            if rec.rtype != gb.gds_record.STRNAME:
                raise GDSII_Error("STRNAME record expected")
            cellname = rec.data
            if cellname in self.cellmap:
                if self.cellmap[cellname] is not None:
                    cellname = self.cellmap[cellname]
                    # Substitute new cellname: make new STRNAME record
                    rec = GDSII_record(gb.gds_strname(cellname))
            cell.header.append(rec)
            # The optional STRCLASS record is assumed to be absent.

            # Loop over this cell's elements
            for e in self.gds_elem_iter(rec_iter, cell):
                cell.addelem(e)
            yield cell


    def gds_elem_iter(self, rec_iter, cell):
        """Return an iterator over all the gds elements in <rec_iter> in <cell>.

        Args:
            rec (GDSII_record iterator): list of records
            cell (GDSII_cell): cell

        Yield:
            GDSII_element:
        """
        elem = GDSII_element()
        for rec in rec_iter:
            if rec.rtype == gb.gds_record.LAYER:
                lay = rec.data[0]
                continue
            elif rec.rtype == gb.gds_record.DATATYPE:
                dtype = rec.data[0]
                layID = (lay, dtype)
                if layID in self.layer_remove:
                    # Remove this element: read until ENDEL and delete.
                    for rec in rec_iter:
                        if rec.rtype == gb.gds_record.ENDEL:
                            elem = GDSII_element()
                            break
                    continue
                if layID in self.layermap: # Replace with new layer number.
                    rec1 = GDSII_record(gb.gds_layer(self.layermap[layID][0]))
                    rec2 = GDSII_record(gb.gds_datatype(self.layermap[layID][1]))
                else:
                    rec1 = GDSII_record(gb.gds_layer(lay))
                    rec2 = GDSII_record(gb.gds_datatype(dtype))
                elem.addrecord(rec1)
                elem.addrecord(rec2)
                lay = -1 #reset lay for next loop
                continue
            elif rec.rtype == gb.gds_record.TEXTTYPE:
                ttype = rec.data[0]
                layID = (lay, ttype)
                if layID in self.layer_remove:
                    # Remove this element: read until ENDEL and delete.
                    for rec in rec_iter:
                        if rec.rtype == gb.gds_record.ENDEL:
                            elem = GDSII_element()
                            break
                    continue
                if layID in self.layermap: # Replace with new layer number.
                    rec1 = GDSII_record(gb.gds_layer(self.layermap[layID][0]))
                    rec2 = GDSII_record(gb.gds_texttype(self.layermap[layID][1]))
                else:
                    rec1 = GDSII_record(gb.gds_layer(lay))
                    rec2 = GDSII_record(gb.gds_texttype(ttype))
                elem.addrecord(rec1)
                elem.addrecord(rec2)
                lay = -1 #reset lay for next loop
                continue
            elif rec.rtype == gb.gds_record.SNAME: # Reference to other cell.
                if rec.data in self.cell_remove: # Remove reference
                    cell.references(rec.data) # But still record the reference.
                    # Remove this element: read until ENDEL and delete.
                    for rec in rec_iter:
                        if rec.rtype == gb.gds_record.ENDEL:
                            elem = GDSII_element()
                            break
                    continue
                if rec.data in self.cellmap:
                    rec = GDSII_record(gb.gds_sname(self.cellmap[rec.data]))
                cell.references(rec.data)
            elif rec.rtype == gb.gds_record.ENDEL: # End of this iter.
                elem.addrecord(rec)
                yield elem
                elem = GDSII_element() # New element
                continue
            elif rec.rtype == gb.gds_record.ENDSTR: # Last element (cell footer)
                break # End of elements for this cell.
            elem.addrecord(rec)


    def gds_record_iter(self, strm, strmlen, pos=0):
        """Return an iterator of the cell records.

        Args:
            strm (byrearray):
            strmlen (int): length of the stream
            pos (int): position in the stream

        Yields:
            GDSII-record
        """
        while pos < strmlen:
            rec = GDSII_record(strm, pos)
            pos += rec.rlen
            yield rec


    def parse(self):
        """Find all the cells in a GDS and add them in an internal dictionary.

        Returns:
            None
        """
        rec_iter = self.gds_record_iter(self.gds_stream, self.stream_length)

        # Read the header: all records until (including) UNITS
        for rec in rec_iter:
            self.header.append(rec)
            if rec.rtype == gb.gds_record.UNITS:
                break # last record of header

        # Read all the cells, ends on ENDLIB record (end of file).
        for cell in self.gds_cell_iter(rec_iter):
            # Add cell reference
            self.cells[cell.name] = cell



    @property
    def libname(self):
        for rec in self.header:
            if rec.rtype == gb.gds_record.LIBNAME:
                return str(rec.data)
        else:
            return None # No LIBNAME record in GDS file header


    @property
    def gdsversion(self):
        # First header record should be HEADER with version number
        if self.header[0].rtype == gb.gds_record.HEADER:
            return self.header[0].data[0]
        else:
            return None
