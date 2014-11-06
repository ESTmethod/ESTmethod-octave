## SpToeSiirdeFiVektor.m
function SpTFiV=SpToeSiirdeFiVektor(VarrasN)
##function SpTFiV=SpToeSiirdeFiVektor(NSARV,NEARV,VarrasN,krdn,selem)
##
##disp('==================================================================')
##disp(' SpToeSiirdeFiVektor.m     ')
##disp(' The vector for transformation the vector [u, w, Fi]'' ')
##disp(' from local to Fiy in global coordinates.')
##disp(' OUTPUT: SpTFiV -- the transformation vector as sparse vector. ')
##
## GNU octave version >= 3.0.x
##=========================================================================
## This Program is writtwn by Andres Lahe,   2010-07-05
##                    e-mail: andres.lahe@ttu.ee
## LAST MODIFIED: Andres Lahe,   2013-07-08
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
##
## NSARV - the number of frame nodes
## NEARV - the number of elements
## krdn - the nodal coordinates
## selem - the topology
## VarrasN - the number of the element
##
#if nargin != 1
if ~(nargin==1)
error(' function SpToeSiirdeFiVektor have wrong number of input arguments!')
end
#
#%selem=[selemjl(:,1:23)];
#EARV=size(selem);
#NEARV=EARV(1,1);
#%%http://www.network-theory.co.uk/docs/octave3/octave_87.html
##if (AF1-eps1) > 0
#switch (VarrasN);
#case{1}
#i=VarrasN;
#
#LkoordN=selem(i,16);
#AkoordN=selem(i,17);
#DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
#DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
#VGRx(i,1)=krdn(AkoordN,1);
#VGRx(i,2)=krdn(LkoordN,1);
#VGRz(i,1)=krdn(AkoordN,2);
#VGRz(i,2)=krdn(LkoordN,2);
#
## lvarras -- the length of the element
#
#lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
#cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
#cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#---------The direction cosines of element -------
suunakosin= zeros(1,3);
        suunakosin(1,1)=0.0;
        suunakosin(1,2)=0.0;
        suunakosin(1,3)=1.0;
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
#
#%%%%%%%%%%%

#endswitch


SpTFiV=SpTM3k3;
#
endfunction
