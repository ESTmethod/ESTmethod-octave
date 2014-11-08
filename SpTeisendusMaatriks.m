##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-07-16
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

%======================================================================
%> @file SpTeisendusMaatriks.m
%> @brief Transformation 3x3 matrix from directional cosines
%> (deprecated, use spTransformation()).
%======================================================================
%> @brief Transformation matrix from directional cosines (deprecated).
%>
%> Generate a transformation matrix from directional cosines.
%> Used to transform vector [N, Q, M] from local to global coordinates.
%> It is deprecated, use spTransformation() with first argument 3 instead.
%>
%> @param node_count Total number of nodes (not used).
%> @param element_count Total number of elements (not used).
%> @param eid Element id.
%> @param coordinates Matrix containing all node coordinates.
%> @param element_properties Topology matrix containing all element properties,
%>            used to determine the beginning and ending nodes of element eid.
%>
%> @retval 3x3 transformation sparse matrix.
%======================================================================
function transform = SpTeisendusMaatriks(node_count, element_count, eid, coordinates, element_properties)
if nargin != 5
    error(' Function SpTeisendusMaatriks has wrong number of input arguments.')
end
transform = spTransformation(3, eid, coordinates, element_properties)
warning('SpTeisendusMaatriks() is deprecated, use spTransformation(3, ...) instead.')
endfunction

