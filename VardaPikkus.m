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
%> @file VardaPikkus.m
%> @brief Calculate lengths of all elements (deprecated, use elementLengths()).
%======================================================================
%> @brief Calculate lengths of all elements (deprecated).
%>
%> Calculate lengths for all elements and return result as a vector.
%> Deprecated, use elementLengths() instead.
%>
%> @param node_count Total number of nodes (unused).
%> @param element_count Total number of elements (unused).
%> @param coordinates Node coordinates in a table.
%> @param elements Elements topology and properties in a huge matrix.
%>
%> @retval element_lengths Vector of lengths of all elements.
%======================================================================
function element_lengths = VardaPikkus(node_count, element_count, coordinates, elements)

if nargin != 4
    error('Function VardaPikkus has wrong number of input arguments.')
end
element_lengths = elementLengths(coordinates, elements)
warning('VardaPikkus() is deprecated. Use elementLengths() instead.')
endfunction

