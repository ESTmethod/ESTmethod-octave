##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2006-02-25
## Last modified: Andres Lahe, 2010-08-05
##                Mattias Põldaru, 2014-11-06
## Copyright (c)  2006 by Tallinn University of Technology
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

function inner_forces = ESTFrKrmus(baasi0, xx ,Li, Fjoud, qkoormus, EA, EI)

#  computed a loading vector for uniformly distributed load
#  and point load (q + F) for a frame
#  EA -- axial stiffness of the element
#  EJ -- bending stiffness of the element
#  baasi0 - scaling multiplier for the displacements (i = EJ/l)

if nargin != 7
    error('Function ESTFrKrmus has wrong number of input arguments.')
end

#L # length of member
#Fjoud # point load and distance of a load point from start point
#qkoormus # uniformly distributed load
## %%%%%%%%%%%%%%%%
# ==========
# Point load and distance of a load point from start point
# Fz Fx aF
# ==========
#Fjoud = [0.0 0.0 0.0; # force 1
# 0.0 0.0 0.0; # force 2
# 0.0 0.0 0.0; # force 3
# 0.0 0.0 0.0; # force 4
# 0.0 0.0 0.0]; # force 5
# ==========
# ==========
# Uniformly distributed load and distances of start point and end point
# qz qx qA qL
# ==========
#qkoormus = [0.0 0.0 0.0 0.0 0.0; # uniformly distributed load qz1
# 0.0 0.0 0.0 0.0 0.0; # uniformly distributed load qz2
# 0.0 0.0 0.0 0.0 0.0]; # uniformly distributed load qz3
# ==========

% Initialize vector for inner forces.
inner_forces = zeros(6, 1);

% Find forces which are not zero.
nonzero_forces = find(Fjoud(:, 1));

for i = nonzero_forces
    [Fz Fx aLx] = num2cell(Fjoud(i, 1:3)){:};
    inner_forces += yzfzv(baasi0, xx, aLx, Fx, Fz, EA, EI);
end

% Find loads which are not zero.
nonzero_loads = find(qkoormus(:, 1));

for i = nonzero_loads
    [qz qx Li] = num2cell(qkoormus(i, 1:3)){:};
    inner_forces += yzhqzm(baasi0, xx, Li, qx, qz, EA, EI);
    #Li = qkoormus(1, 4); # aqL - kui koormus lõppeb enne lõppu
    #Li = qkoormus(1, 5); # varda lõpp (vaja kohendada)
    #Zq0(1, 1) = - qx * x;
    #Zq0(2, 1) = - qz * x;
    #Zq0(3, 1) = - qz * x^2 / 2;
    # vaja lisada juht, kui koormus lõppeb enne lõppu
end
endfunction

