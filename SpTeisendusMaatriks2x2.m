## SpTeisendusMaatriks2x2.m
function SpTM2x2=SpTeisendusMaatriks2x2(NSARV,NEARV,VarrasN,krdn,selem)
##
##
##disp('==================================================================')
##disp(' SpTeisendusMaatriks2x2.m     ')
##disp(' The matrix for transformation the vector [N, Q]'' ')
##disp(' from local to global coordinates.')
##disp(' OUTPUT: SpTM2x2 -- the transformation Matrix as sparse matrix. ')
##
## GNU octave version >= 3.0.x
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
## along with this program; if not, write to the Free Software Foundation, Inc.
## 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
##=========================================================================
##
## NSARV - the number of frame nodes
## NEARV - the number of elements
## krdn - the nodal coordinates
## selem - the topology
## VarrasN - the number of the element
##
#if nargin != 5
if ~(nargin==5)
error(' function SpTM2x2=SpTeisendusMaatriks2x2 have wrong number of input arguments!')
end
#
##selem=[selemjl(:,1:23)];
EARV=size(selem);
NEARV=EARV(1,1);
#%%http://www.network-theory.co.uk/docs/octave3/octave_87.html
#if (AF1-eps1) > 0
switch (VarrasN);
case{1}
i=VarrasN;
%
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
%
# lvarras -- the length of the element
%
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
%
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
%
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

%suunakosin
%TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
%
#%%%%%%%%%%%
case{2}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#
#%%%%%%%%%%%


case{3}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#
#%%%%%%%%%%%


case{4}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#
#%%%%%%%%%%%


case{5}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#
#%%%%%%%%%%%


case{6}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#
#%%%%%%%%%%%


case{7}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#
#%%%%%%%%%%%

case{8}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#
#%%%%%%%%%%%

case{9}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#
#%%%%%%%%%%%

case{10}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#
#%%%%%%%%%%%

case{11}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{12}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{13}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{14}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{15}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{16}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{17}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{18}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{19}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{20}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{21}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{22}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{23}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{24}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{25}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{26}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{27}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{28}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{29}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{30}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{31}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{32}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{33}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{34}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{35}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{36}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{37}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{38}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{39}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{40}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{41}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{42}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{43}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{44}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{45}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{46}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{47}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{48}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{49}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#

case{50}
i=VarrasN;
#
LkoordN=selem(i,16);
AkoordN=selem(i,17);
DeltaX(i)=krdn(LkoordN,1)-krdn(AkoordN,1);
DeltaZ(i)=krdn(LkoordN,2)-krdn(AkoordN,2);
VGRx(i,1)=krdn(AkoordN,1);
VGRx(i,2)=krdn(LkoordN,1);
VGRz(i,1)=krdn(AkoordN,2);
VGRz(i,2)=krdn(LkoordN,2);
#
# lvarras -- the length of the element
#
lvarras(i,1)=sqrt(DeltaX(i)^2+DeltaZ(i)^2);
#
cosAlpha(i,1)=DeltaX(i)/lvarras(i,1);
cosBeta(i,1)=DeltaZ(i)/lvarras(i,1);
#
suunakosin= zeros(2,2);
        suunakosin(1,1)=cosAlpha(i,1);
        suunakosin(1,2)=-cosBeta(i,1);

        suunakosin(2,1)=cosBeta(i,1);
        suunakosin(2,2)=cosAlpha(i,1);

#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM2x2=suunakosin;
SpTM2k2=sparse(TM2x2);
#


case{51}
case{52}
case{53}
case{54}
case{55}
case{56}
case{57}
case{58}
case{59}
case{60}

endswitch


SpTM2x2=SpTM2k2;
#
endfunction
