##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2013-06-26
## Last modified: Andres Lahe, 2013-06-27
##                Mattias Põldaru, 2014-11-08
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

%======================================================================
%> @file VardadSolmes.m
%> @brief Find and sort elements connected to a given node (deprecated, use elementsInNode()).
%======================================================================
%> @brief Find and sort elements connected to a given node (deprecated).
%>
%> Find all elements which are connected to a given node, sort elements by
%> nids (ABB column 8), flip the results and sort by column 11.
%> This function is deprecated, use elementsInNode() instead.
%>
%> @param node_count Total number of nodes (unused).
%> @param element_count Total number of elements (unused).
%> @param nid The original matrix.
%> @param AB A matrix indicating which elements are connected to nodes.
%> @param ABB A element properties matrix.
%>            [eid, u w fi N Q M, nid, hinges N Q M]
%>
%> @retval elements All elements which are connected to node nid.
%======================================================================
function elements = VardadSolmes(node_count, element_count, nid, AB, ABB)
if nargin != 5
    error(' function Vardadnides has wrong number of input arguments!')
end

elements = elementsInNode(nid, AB, ABB)
warning('VardadSolmes() is deprecated. Use elementsInNode() instead.')
endfunction
