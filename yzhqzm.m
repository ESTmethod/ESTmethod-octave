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
## along with this program; if not, write to the Free Software Foundation, Inc.
## 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
##=========================================================================

function Zq0 = yzhqzm(baasi0, x, a, qx, qz, EA, EJ)
#  Zq0 -The load vector for uniformly distriputed load q (projections qx and qz)
#  in the tranfer matrix method    ZL - U*ZA =Zq0 (the II sign convention).
#  http://digi.lib.ttu.ee/opik_eme/Ehitusmehaanika.pdf#page=396  page=694
#  EA - the extension stiffness of the element
#  EJ -  the bending stiffness of the element
#  x - coordinate
#  qx - the uniformly distributed load in x direction
#  qz - the uniformly distributed load in z direction
#  baasi0 - scaling multiplier for the displacements (io= EJo/lo)

if nargin != 7
    error(' Function yzhqz has wrong number of input arguments!')
end

i0 = baasi0;

if (x - a) >= 0
    Zq0(1, 1) = - i0 * qx * x^2 / (2 * EA);
    Zq0(2, 1) = i0 * qz * x^4 / (24 * EJ);
    Zq0(3, 1) = - i0 * qz * x^3 / (6 * EJ);
    Zq0(4, 1) = - qx * x;
    Zq0(5, 1) = - qz * x;
    Zq0(6, 1) = - qz * x^2 / 2;
else
    Zq0 = zeros(6, 1)
endif
endfunction

