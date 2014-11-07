##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-06-26
## Last modified: Andres Lahe, 2013-07-08
## Copyright (c)  2010 by Tallinn University of Technology
##                Department of Mechanics
##                http://www.ttu.ee/
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
##=========================================================================

function lvarras = VardaPikkus(NSARV, NEARV, krdn, selem)

## NSARV - the number of frame nodes
## NEARV - the number of elements
## krdn - the nodal coordinates
## selem - the topology
## lvarras - the length of the element
##
if nargin != 4
    error(' Function VardaPikkus has wrong number of input arguments.')
end

NEARV = size(selem)(1);
for i = 1:NEARV
    LkoordN = selem(i, 16);
    AkoordN = selem(i, 17);
    dX = krdn(LkoordN, 1) - krdn(AkoordN, 1);
    dZ = krdn(LkoordN, 2) - krdn(AkoordN, 2);
    lvarras(i, 1) = sqrt(dX^2 + dZ^2);
endfor
endfunction

