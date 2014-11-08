##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-07-16
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
## along with this program; if not, write to the Free Software Foundation, Inc.
## 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
##=========================================================================

%======================================================================
%> @file spTransformation.m
%> @brief Transformation matrix of directional cosines.
%======================================================================
%> @brief Transformation matrix of directional cosines.
%>
%> Generate a transformation matrix of directional cosines.
%> Used to transform vector [N, Q, M] from local to global coordinates.
%>
%> @param dim Dimension of output matrix. Can be either 2 or 3.
%> @param eid Element id.
%> @param coordinates Matrix containing all node coordinates.
%> @param element_properties Topology matrix containing all element properties,
%>            used to determine the beginning and ending nodes of element eid.
%>
%> @retval Newly created tranformational sparse matrix, either 2x2 or 3x3.
%======================================================================
function transform = spTransformation(dim, eid, coordinates, element_properties)

if nargin != 4
    error(' Function spTransformation() has wrong number of input arguments.')
end

if (dim < 2) or (dim > 3)
    error(' Function spTransformation() first argument dimension has to be either 2 or 3.')
end

end_nid = element_properties(eid, 16);
start_nid = element_properties(eid, 17);
dX = coordinates(end_nid, 1) - coordinates(start_nid, 1);
dZ = coordinates(end_nid, 2) - coordinates(start_nid, 2);
len = sqrt(dX^2 + dZ^2);
cosAlpha = dX / len;
cosBeta = dZ / len;
# The direction cosines of element
dir_cosine = zeros(dim, dim);
dir_cosine(1, 1) = cosAlpha;
dir_cosine(1, 2) = -cosBeta;
dir_cosine(2, 1) = cosBeta;
dir_cosine(2, 2) = cosAlpha;
if dim == 3
    dir_cosine(3, 3) = 1.0;
endif
transform = sparse(dir_cosine);
endfunction

