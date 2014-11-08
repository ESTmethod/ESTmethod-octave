##=========================================================================
## Copyright (c)  2013 by Tallinn University of Technology
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2013-06-26
## Last modified: Andres Lahe, 2013-06-27
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
%> @file SpTeisendusUhikMaatriks.m
%> @brief Create a sparse identity matrix 3x3 (deprecated,
%> use speye(3)).
%======================================================================
%> @brief Create a sparse identity matrix 3x3 (deprecated).
%>
%> Create a sparse identity matrix 3x3. It has been deprecated in
%> favour of built-in speye() function, so use speye(3) instead.
%>
%> @retval New 3x3 identity matrix.
%======================================================================
function identity = SpTeisendusUhikMaatriks()
##   The sparse identity matrix of reaction at the node
identity = speye(3);
warning('SpTeisendusUhikMaatriks() is deprecated. Use speye(3) instead.')
endfunction

