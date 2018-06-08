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
#
# Write the klayout layer properties file
#
# (c) 2016 Xaveer Leijtens
# 
def lypwrite(df, filename):
    lyp_fields = ["frame-color", "fill-color", "frame-brightness",
            "fill-brightness", "dither-pattern", "valid", "visible",
            "transparent", "width", "marked", "animation"]

    with open(filename, 'w') as lyp:
        lyp.write('<?xml version="1.0" encoding="utf-8"?>\n')
        lyp.write('<layer-properties>\n')
        for index, row in df.iterrows():
            lyp.write(' <properties>\n')
            for field in lyp_fields:
                val = row[field]
                lyp.write('  <%s>%s</%s>\n' % (field, val, field))
            lyp.write('  <name>%s-%s</name>\n' % (row['Number'], row['Name']))
            lyp.write('  <source>%s/%s@1</source>\n' % \
                (row['Number'],row['Datatype']))
            lyp.write(' </properties>\n')
        lyp.write(' <name>COBRA</name>\n</layer-properties>\n')
