## SpTeisendusUhikMaatriks2x2.m
function SpTM2x2xz=SpTeisendusUhikMaatriks2x2(VarrasN)
##
##   The sparse identity matrix of reaction at the node
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
##
## NSARV - the number of frame nodes
## NEARV - the number of elements
## VarrasN - the number of the element
##
#if nargin != 1
if ~(nargin==1)
error(' function SpTeisendusUhikMaatriks2x2 have wrong number of input arguments!')
end
#
#%%http://www.network-theory.co.uk/docs/octave3/octave_87.html
VarrasN=1;
#if (AF1-eps1) > 0
switch (VarrasN);
case{1}
i=VarrasN;
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=1.0;
        suunakosin(1,2)=0.0;
        suunakosin(2,1)=0.0;
        suunakosin(2,2)=1.0;
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#
#%%%%%%%%%%%

endswitch


SpTM2x2xz=SpTM2k2;
#
endfunction
