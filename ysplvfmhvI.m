##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2009-02-13
## Last modified: Andres Lahe, 2010-03-12
## Copyright (c)  2009 by Tallinn University of Technology
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

function spvF = ysplvfmhvI(baasi0, x, l, EA, GAr, EJ)
# The basic sparse matrix spvF (6,12) for a frame element with the transfer matrix U.
#     spvF=Zp-U*Zv, Zp'=[up wp fip Np Qp Mp]
# The Sign Convention II for normal-, shear force and bending moment
# ---------------------------------------------
# W. D. Pilkey, W. Wunderlich. Mechanics of Structures. Variational and
# Computational Methods. CRC Press, Boca Raton, Ann Arbor, London, Tokyo, 1994.
# 4.3.1 Transfer Matrix (the basic equations) p. 174
#----------------------------------------------
#  EA - the extension stiffness of the element
#  GAr - the shear stiffness of the element
#  EJ -  the bending stiffness of the element
#  x - coordinate
#  l - the length of the element
#  baasi0 - scaling multiplier for the displacements (io= EJo/lo)

if nargin != 6
    error(' function ysplvfmhvI has wrong number of input arguments!')
end

spvF = sparse(6, 12);
spvI = speye(6);
spvF = spInsertBtoA(spvF, 1, 1, spvI);
spvU1 = ysplfhlin(baasi0, x, EA, GAr, EJ);
spvU = spvU1 .* (-1);
spvF = spInsertBtoA(spvF, 1, 7, spvU);

endfunction

