##=========================================================================
## Copyright (c)  2002 by Tallinn University of Technology
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2002-05-15
## Last modified: Andres Lahe, 2009-01-03
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

function sA = spInsertBtoA(spA, IM, JN, spB)
#
# Inserted sparse matrix spB  into sparse matrix spA,
# starting at row index IM and column index JN.
# When matrix spA and matrix spB elements overlaped
# then they added together.

#% 1 row, then   iiB = iiB( : )'; jjB = jjB( : )' ;
#% colunm, then   iiB = iiB( : ); jjB = jjB( : );
#%%%%%%%%%%%%%%%%%
if nargin != 4
    error(' function spInsertBtoA has wrong number of input arguments!')
end
#
OnHore = issparse(spA);
if OnHore == 0
    disp('!!!!!')
    spA
    disp('!!!!!')
    disp('The matrix spA(M,N) must be sparse!')
    error('The sparse matrix spA(M,N) must be instantiated!')
    ## return
end
# %%%%%%%%%%%%%%%%%
#
IMB = IM - 1;
JNB = JN - 1;
#
[iiA, jjA, aaA] = find(spA);
[iiB, jjB, aaB] = find(spB);
iiB = iiB .+ IMB;
jjB = jjB .+ JNB;
iiB = iiB(:);
jjB = jjB(:);
aaB = aaB(:);
iis = [iiA; iiB];
jjs = [jjA; jjB];
aas = [aaA; aaB];

Auus = sparse(iis, jjs, aas);
#
sA = Auus;
#
endfunction

