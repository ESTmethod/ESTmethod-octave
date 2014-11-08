##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-07-05
## Last modified: Andres Lahe, 2013-07-08
##                Mattias Põldaru, 2014-11-06
## Copyright (c)  2010 by Tallinn University of Technology
##                Department of Mechanics
##                http://www.ttu.ee/
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software Foundation, Inc.
## 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
##=========================================================================

function SpTM2x2 = SpTeisendusMaatriks2x2(node_count, element_count, eid, coordinates, element_properties)
##disp(' The matrix for transformation the vector [N, Q]'' ')
##disp(' from local to global coordinates.')
##disp(' OUTPUT: SpTM2x2 -- the transformation Matrix as sparse matrix. ')

## node_count - the number of frame nodes
## element_count - the number of elements
## coordinates - the nodal coordinates
## element_properties - the topology
## eid - the number of the element

if nargin != 5
    error(' Function SpTM2x2=SpTeisendusMaatriks2x2 has wrong number of input arguments.')
end

end_nid = element_properties(eid, 16);
start_nid = element_properties(eid, 17);
dX = coordinates(end_nid, 1) - coordinates(start_nid, 1);
dZ = coordinates(end_nid, 2) - coordinates(start_nid, 2);
len = sqrt(dX^2 + dZ^2);
cosAlpha = dX / len;
cosBeta = dZ / len;
# The direction cosines of element
dir_cosine = zeros(2, 2);
dir_cosine(1, 1) = cosAlpha;
dir_cosine(1, 2) = -cosBeta;
dir_cosine(2, 1) = cosBeta;
dir_cosine(2, 2) = cosAlpha;
SpTM2x2 = sparse(dir_cosine);
endfunction

