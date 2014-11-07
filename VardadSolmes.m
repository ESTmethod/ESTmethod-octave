##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2013-06-26
## Last modified: Andres Lahe, 2013-06-27
## Copyright (c)  2013 by Tallinn University of Technology
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
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
##=========================================================================

function Sv = VardadSolmes(NSARV, NEARV, Solm, AB, ABB)
#
# Sort the rows of matrix ABB according to the order of the (column=8) node
# number. Returned matrix with the order of the row reversed.
#
if nargin != 5
    error(' function VardadSolmes has wrong number of input arguments!')
end
## NSARV - the number of frame_nodes
## NEARV - the number of elements
## Solm - the node number
## AB - the elements connected at the nodes
## ABB - matrix[Element number, u w fi N Q M, Node, hinges N Q M]

i = Solm;
k = 1;
for j = 1 : size(AB)
    if AB(j) == i
        ABB2(k, :) = ABB(j, :);
        k += 1;
    endif
endfor

[!, idx] = sortrows(ABB2(:, 11));
Sv = ABB2(idx, :);
endfunction

