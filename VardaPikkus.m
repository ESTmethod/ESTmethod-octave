## VardaPikkus.m
function lvarras=VardaPikkus(NSARV,NEARV,krdn,selem)
##
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
## NSARV - the number of frame nodes
## NEARV - the number of elements
## krdn - the nodal coordinates
## selem - the topology
## lvarras - the length of the element
##
#if nargin != 4
if ~(nargin==4)
error(' function VardaPikkus have wrong number of input arguments!')
end
#
EARV=size(selem);
NEARV=EARV(1,1);
for i=1:NEARV
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
endfor
#
#Varda pikkus
#lvarras
for i=1:NEARV
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
endfor
#  ---------The direction cosines of element -------

#for i=1:NEARV
#cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
#cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#endfor
#

lvarras;

#
endfunction
