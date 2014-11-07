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
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
##=========================================================================

function SpTUv = SpToeSiirdeUvektor(NSARV, NEARV, n, krdn, selem)
##disp(' The vector for transformation the vector [u, w, Fi]'' ')
##disp(' from local to Ux in global coordinates. ')
##disp(' OUTPUT: SpTUv -- the transformation vector as sparse vector. ')
##
## NSARV - the number of frame nodes
## NEARV - the number of elements
## krdn - the nodal coordinates
## selem - the topology
## n - the number of the element

if nargin != 5
    error(' function SpToeSiirdeUvektor has wrong number of input arguments!')
end

LkoordN = selem(n, 16);
AkoordN = selem(n, 17);
dX = krdn(LkoordN, 1) - krdn(AkoordN, 1);
dZ = krdn(LkoordN, 2) - krdn(AkoordN, 2);
len = sqrt(dX^2 + dZ ^ 2);
# --------- The direction cosines of element -------
suunakosin = zeros(1, 3);
suunakosin(1, 1) = dX / len; # Cosalpha
suunakosin(1, 2) = dZ / len; # Cosbeta
suunakosin(1, 3) = 0.0;
SpTUv = sparse(suunakosin);
endfunction

