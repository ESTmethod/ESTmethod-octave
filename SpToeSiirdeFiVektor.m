##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-07-05
## Last modified: Andres Lahe, 2013-07-08
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

function SpTFiV = SpToeSiirdeFiVektor()
##disp(' The vector for transformation the vector [u, w, Fi]'' ')
##disp(' from local to Fiy in global coordinates.')
##disp(' OUTPUT: SpTFiV -- the transformation vector as sparse vector. ')

# This returns a sparse [0 0 1] matrix.
SpTFiV = sparse(1, 3, 1);
endfunction

