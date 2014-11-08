##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-06-26
## Last modified: Andres Lahe, 2013-07-08
##                Mattias Põldaru, 2014-11-08
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
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
##=========================================================================

%======================================================================
%> @file elementLengths.m
%> @brief Calculate lengths of all elements.
%======================================================================
%> @brief Calculate lengths of all elements.
%>
%> Calculate lengths for all elements and return result as a vector.
%>
%> @param coordinates Node coordinates in a table.
%> @param elements Elements topology and properties in a huge matrix.
%>
%> @retval element_lengths Vector of lengths of all elements.
%======================================================================
function element_lengths = elementLengths(coordinates, elements)
if nargin != 2
    error('Function VardaPikkus has wrong number of input arguments.')
end

element_count = size(elements)(1);
for i = 1:element_count
    end_nid = elements(i, 16);
    start_nid = elements(i, 17);
    dX = coordinates(end_nid, 1) - coordinates(start_nid, 1);
    dZ = coordinates(end_nid, 2) - coordinates(start_nid, 2);
    element_lengths(i, 1) = sqrt(dX^2 + dZ^2);
endfor
endfunction

