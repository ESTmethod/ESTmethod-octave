##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 1998-05-20
## Last modified: Andres Lahe, 2009-02-13
## Copyright (c)  2004 by Tallinn University of Technology
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

function Fzx = yzfzv(scale, x, a, Fx, Fz, EA, EJ)
#  Fzx -The load vector  F (projections Fx and Fz)
#  in the tranfer matrix method   ZL - U*ZA =Fzx (the II sign convention).
#  http://digi.lib.ttu.ee/opik_eme/Ehitusmehaanika.pdf#page=396  page=694
#  EA - the extension stiffness of the element
#  EJ -  the bending stiffness of the element
#  x - coordinate
#  a - the point of load
#  qx - the uniformly distributed load in x direction
#  qz - the uniformly distributed load in z direction
#  scale - scaling multiplier for the displacements (io= EJo/lo)
#

if nargin != 7
    error('Function yzfzv has wrong number of input arguments!')
end

Fzx = zeros(6, 1);
xp = x - a;
if xp >= 0
    Fzx(1, 1) = - scale * Fx * xp / EA;
    Fzx(2, 1) = scale * Fz * xp^3 / (6 * EJ);
    Fzx(3, 1) = - scale * Fz * xp^2 / (2 * EJ);
    Fzx(4, 1) = - Fx;
    Fzx(5, 1) = - Fz;
    Fzx(6, 1) = - Fz * xp;
endif
endfunction

