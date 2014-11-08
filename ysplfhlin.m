##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2009-02-14
## Last modified: Andres Lahe, 2010-03-12
## Copyright (c)  2009 by Tallinn University of Technology
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

function spF = ysplfhlin(scale, x, EA, GAr, EJ)
#
# The sparse matrix spF (6,6) for linear equations of frame element with the transfer matrix U.
#  Zx=U*Zo+lq
#  Zo'=[uo wo fio No Qo Mo]
#  spF=sparse(U)
# The Sign Convention II for normal-, shear force and bending moment
#  http://digi.lib.ttu.ee/opik_eme/Ehitusmehaanika.pdf#page=397   page=693
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
#  scale - scaling multiplier for the displacements (io= EJo/lo)
#

if nargin != 5
    error('Function ysplfhlin has wrong number of input arguments!')
end

F = zeros(6, 6);
F(1, 1) = 1.0;
F(2, 2) = 1.0;
F(3, 3) = 1.0;
F(1, 4) = -scale * x / EA;
F(2, 3) = -x;
F(2, 5) = scale * x^3 / (6 * EJ);
F(2, 6) = scale * x^2 / (2 * EJ);
F(3, 5) = -scale * x^2 / (2 * EJ);
F(3, 6) = -scale * x / EJ;
F(4, 4) = -1;
F(5, 5) = -1;
F(6, 5) = -x;
F(6, 6) = -1;
spF = sparse(F);
endfunction

