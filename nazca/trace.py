#!/usr/bin/env python3
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

# You should have received a copy of the GNU Affero General Public License
# along with Nazca.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------
#==============================================================================
# (c) 2017 Ronald Broeke
#==============================================================================

"""Module for tracing path lengths."""


from collections import defaultdict

trace_closed = None
trace_cnt = -1
trace_id_list = []
traces = defaultdict(list)

def trace_start(name=None):
    """Start recording a trace.

    If no trace_id is provided an ordinal number will be selected.

    Args:
        name (str | int): optional id of the trace to stop

    Returns:
        None
    """
    global trace_id, trace_cnt, trace_ids
    trace_cnt += 1
    if name is None:
        while trace_cnt in trace_id_list:
            trace_cnt += 1
        name = trace_cnt
    elif name in trace_id_list:
        print("Warning: starting already active trace_id = {} in trace_start.".\
           format(name))

        return
    trace_id_list.append(name)
    #print('tracelist:', trace_id_list)


def trace_stop(name=None):
    """Stop recording on trace.

    If no trace_id is given the last started trace is stopped.

    Args:
        trace_id (str | int): id of the trace to stop

    Returns:
        None
    """
    global trace_closed
    if name is None:
         name = trace_id_list[-1]
    #print('stop trace:', name)
    trace_closed = name
    if name == trace_id_list[-1]:
        del trace_id_list[-1]
        try:
            trace_id_list[-1]
        except:
            pass
    else:
        print("Warning: Trying to stop already stopped trace id={}.".\
            format(trace_id))


def trace_append(elm):
    traces[trace_id_list[-1]].append(elm)


def trace_length(name=None):
    if name is None:
        name = trace_closed
    length = 0
    #print('trace length id:', name)
    for elm in traces[name]:
        try:
            L = elm.cnode.cell.length_geo
            length += L
            #print(name, elm, L)
        except AttributeError:
            pass
    return length
