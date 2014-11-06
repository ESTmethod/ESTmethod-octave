## SpToeSiirdeWvektor.m
function SpTWv=SpToeSiirdeWvektor(NSARV,NEARV,VarrasN,krdn,selem)
##
## 
##disp('==================================================================')
##disp(' SpToeSiirdeWvektor.m     ')
##disp(' The vector for transformation the vector [u, w, Fi]'' ')
##disp(' from local to Wz in global coordinates. ')
##disp(' OUTPUT: SpTWv -- the transformation vector as sparse vector. ')
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
#if nargin != 5
if ~(nargin==5)
error(' function SpToeSiirdeWvektor have wrong number of input arguments!') 
end
#
#%selem=[selemjl(:,1:23)];
EARV=size(selem);
NEARV=EARV(1,1);
#%%http://www.network-theory.co.uk/docs/octave3/octave_87.html
#%kui tahad, et arvud oleksid alati 1 reas, siis   iiB = iiB( : )'; jjB = jjB( : )' ;
#%kui tahad, et arvud oleksid alati tulbas, siis   iiB = iiB( : ); jjB = jjB( : );

%if (AF1-eps1) > 0
switch (VarrasN);
case{1}
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
#---------The direction cosines of element -------
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
#
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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
suunakosin= zeros(1,3);
         
         
                 
        suunakosin(1,1)=-cosBeta(i,1);  
        suunakosin(1,2)=cosAlpha(i,1);  
        suunakosin(1,3)=0.0;        
          
           
                      
#suunakosin
#TJ=suunakosin(NEARV,cosAlpha(i,1),cosBeta(i,2));

TM3x3=suunakosin;
SpTM3k3=sparse(TM3x3);
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


SpTWv=SpTM3k3;
#
endfunction
