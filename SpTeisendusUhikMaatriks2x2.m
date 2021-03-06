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
%> @file SpTeisendusUhikMaatriks2x2.m
%> @brief Create a sparse identity matrix 2x2 (deprecated,
%> use speye(2)).
%======================================================================
%> @brief Create a sparse identity matrix 2x2 (deprecated).
%>
%> Create a sparse identity matrix 2x2. It has been deprecated in
%> favour of built-in speye() function, so use speye(2) instead.
%>
%> @retval New 2x2 identity matrix.
%======================================================================
function identity = SpTeisendusUhikMaatriks2x2()
##
##   The sparse identity matrix of reaction at the node
##
identity = speye(2)
warning('SpTeisendusUhikMaatriks2x2() is deprecated. Use speye(2) instead.')
endfunction

