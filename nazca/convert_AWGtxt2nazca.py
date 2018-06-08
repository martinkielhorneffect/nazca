# -*- coding: utf-8 -*-
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
"""
Convert BrightAWG coordinates.txt to nazca.txt

Reformat the AWG coordinates from the BrightAWG .txt outfile into Nazca format for all .txt. files in a specified directory.

    Rename port 'in' -> 'a'
    Rename port 'out' -> 'b'
    rotate ports 'in' 180 degree
    export as ...nazca.txt

Author: Ronald Broeke 2017(c)
"""

import os


def get_txt_files(dir):
    """Get all .txt. files in directory path 'dir'.

    Args:
        dir (str): root directory to scan

    Returns:
        list of str: list of filenames
    """
    root, dirs, files = next(os.walk(dir))
    txt = [os.path.join(root, f) for f in files if 's.txt' in f]
    return txt


def convert_file(file, xs, win, wout):
    """Reformat the AWG coordinates from the OD BrightAWG .txt file into Nazca format.

    Output format: port, x, y, a, xs, w

    A header row is included in the output
    OD port names are renamed to nazca pin names, e.g. in0 -> a0, out4 -> b4.
    Also, input port will undergo a 180 degree rotation to adhere to the Nazca
    phylosophy of "outward pointing" pins.

    Args:
        file (str): filename
        xs (float): xsection of connection
        win (float): waveguide width input side
        wout (float): waveguide width output side

    Returns:
        list: list containing table with pin data
    """
    skiprows = 3
    replacing = [('[', ''), (']', ''), (' ', ''), (':', ','), ('\n', '')]
    output = ['port, x, y, a, xs, w']
    with open(file, 'r') as f:
        for i, line in enumerate(f):
            if i >= skiprows:
                for mapping in replacing:
                    line = line.replace(mapping[0], mapping[1])
                if 'in' in line:
                    line = line.replace('in', 'a')
                    parts = line.split(',')
                    parts[3] = str(float(parts[3])+180)
                    line = ','.join(parts)
                    width = win
                if 'out' in line:
                    line = line.replace('out', 'b')
                    width = wout
                if 'fpri' in line or 'fpro' in line:
                    parts = line.split(',')
                    parts[3] = str(float(parts[3])+180)
                    line = ','.join(parts)
                    xs = None
                    width = 0.0
                line = '{},{},{}'.format(line, xs, width)
                output.append(line)
    return output


def convert(files, xs, win, wout):
    """Loop .txt input files and write converted nazca output files.

    Args:
        files (str | list of str): filename(s)
        xs (float): xsection of connection
        win (float): waveguide width input side
        wout (float): waveguide width output side

    Returns:
        None
    """

    if isinstance(files, str):
        files = [files]

    for name in files:
        if 'coordinates.txt' in name:
            newname = name[:-4]+'_nazca.txt'
            print(newname)
            updated = convert_file(name, xs, win, wout)
            with open(newname, 'w') as nf:
                for line in updated:
                    nf.write(line + '\n')
    return None


if __name__ == '__main__':
    path = '.'
    files = get_txt_files(path)
    xs = 'Deep'
    win = '1.5'
    wout = '2.0'
    convert(files, xs, win, wout)
