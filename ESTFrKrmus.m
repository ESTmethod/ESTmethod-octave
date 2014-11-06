# ESTFrKrmus.m
function Sj=ESTFrKrmus(baasi0,xx,Li,Fjoud,qkoormus,EA,EI)
#
#  computed a loading vector for uniformly distributed load 
#  and point load (q + F) for a frame
#  EA -- axial stiffness of the element
#  EJ -- bending stiffness of the element
#  baasi0 - scaling multiplier for the displacements (i = EJ/l)
#
##=========================================================================
## This Program is written by Andres Lahe,    2006-02-25
##                    e-mail: andres.lahe@ttu.ee
## LAST MODIFIED: Andres Lahe,    2010-08-05
## Copyright (c)  2006  by Tallinn University of Technology 
##                Department of Mechanics
##                http://www.ttu.ee/
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
##========================================================================= 
#
#if nargin != 7
if ~(nargin==7)
error(' function spInsertBtoA have wrong number of input arguments!') 
end
#
#
#L         # length of member
#Fjoud     # point load and distance of a load point  from start point
#qkoormus  # uniformly distributed load 
##%%%%%%%%%%%%%%%%
#l=Li;
Li=Li;
xx=xx;
#==========
# Point load and distance of a load point  from start point
#       Fz    Fx     aF
#==========
#Fjoud=[0.0    0.0   0.0;	# force 1
#	 0.0    0.0   0.0;	# force 2
#	 0.0    0.0   0.0;	# force 3
#	 0.0    0.0   0.0;	# force 4
#	 0.0    0.0   0.0];	# force 5
#==========
#==========
# Uniformly distributed load and distances of start point and end point
#       qz     qx     qA      qL
#==========
#qkoormus=[0.0   0.0   0.0  0.0   0.0;	# uniformly distributed load qz1
#          0.0   0.0   0.0  0.0   0.0;	# uniformly distributed load qz2
#          0.0   0.0   0.0  0.0   0.0];	# uniformly distributed load qz3
#==========
Fjoud;
qkoormus;
 eps1=0.000001;
FARV=size(Fjoud);
NFARV=FARV(1,1);
#
F1=Fjoud(1,1);
AF1=abs(F1);
FjouduS=size(Fjoud);
Fjoudu=FjouduS(1,1);
FZZZ=Fjoud(:,1);
FZZ=find(FZZZ);
FZY=size(FZZ);
Fjoudu=FZY(1,1);
#%%http://www.network-theory.co.uk/docs/octave3/octave_87.html
if (AF1-eps1) > 0
switch (Fjoudu);
case{1}
Fz1=Fjoud(1,1);
Fx1=Fjoud(1,2);
aLx1=Fjoud(1,3);
vFz1=yzfzv(baasi0,xx,aLx1,Fx1,Fz1,EA,EI);
vFz=vFz1;
case{2}
Fz1=Fjoud(1,1);
Fx1=Fjoud(1,2);
aLx1=Fjoud(1,3);
vFz1=yzfzv(baasi0,xx,aLx1,Fx1,Fz1,EA,EI);
Fz2=Fjoud(2,1);
Fx2=Fjoud(2,2);
aLx2=Fjoud(2,3);
vFz2=yzfzv(baasi0,xx,aLx2,Fx2,Fz2,EA,EI);
vFz=vFz1+vFz2;
case{3}
Fz1=Fjoud(1,1);
Fx1=Fjoud(1,2);
aLx1=Fjoud(1,3);
vFz1=yzfzv(baasi0,xx,aLx1,Fx1,Fz1,EA,EI);
Fz2=Fjoud(2,1);
Fx2=Fjoud(2,2);
aLx2=Fjoud(2,3);
vFz2=yzfzv(baasi0,xx,aLx2,Fx2,Fz2,EA,EI);
Fz3=Fjoud(3,1);
Fx3=Fjoud(3,2);
aLx3=Fjoud(3,3);
vFz3=yzfzv(baasi0,xx,aLx3,Fx3,Fz3,EA,EI);
vFz=vFz1+vFz2+vFz3;
case{4}
Fz1=Fjoud(1,1);
Fx1=Fjoud(1,2);
aLx1=Fjoud(1,3);
vFz1=yzfzv(baasi0,xx,aLx1,Fx1,Fz1,EA,EI);
Fz2=Fjoud(2,1);
Fx2=Fjoud(2,2);
aLx2=Fjoud(2,3);
vFz2=yzfzv(baasi0,xx,aLx2,Fx2,Fz2,EA,EI);
Fz3=Fjoud(3,1);
Fx3=Fjoud(3,2);
aLx3=Fjoud(3,3);
vFz3=yzfzv(baasi0,xx,aLx3,Fx3,Fz3,EA,EI);
Fz4=Fjoud(4,1);
Fx4=Fjoud(4,2);
aLx4=Fjoud(4,3);
vFz4=yzfzv(baasi0,xx,aLx4,Fx4,Fz4,EA,EI);
vFz=vFz1+vFz2+vFz3+vFz4;
case{5}
Fz1=Fjoud(1,1);
Fx1=Fjoud(1,2);
aLx1=Fjoud(1,3);
vFz1=yzfzv(baasi0,xx,aLx1,Fx1,Fz1,EA,EI);
Fz2=Fjoud(2,1);
Fx2=Fjoud(2,2);
aLx2=Fjoud(2,3);
vFz2=yzfzv(baasi0,xx,aLx2,Fx2,Fz2,EA,EI);
Fz3=Fjoud(3,1);
Fx3=Fjoud(3,2);
aLx3=Fjoud(3,3);
vFz3=yzfzv(baasi0,xx,aLx3,Fx3,Fz3,EA,EI);
Fz4=Fjoud(4,1);
Fx4=Fjoud(4,2);
aLx4=Fjoud(4,3);
vFz4=yzfzv(baasi0,xx,aLx4,Fx4,Fz4,EA,EI);
Fz5=Fjoud(5,1);
Fx5=Fjoud(5,2);
aLx5=Fjoud(5,3);
vFz5=yzfzv(baasi0,xx,aLx5,Fx5,Fz5,EA,EI);
vFz=vFz1+vFz2+vFz3+vFz4+vFz5;
#%?otherwise
endswitch
else
# disp(' Koondatud force puudub  '); 
vFz=[0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
#vFz=vFz'
endif

#%%%%%%%%%%%%%%%%%%%%%
qzv=qkoormus(1,1);
Aq1=abs(qzv);
qZZZ=qkoormus(:,1);
qZZ=find(qZZZ);
qZY=size(qZZ);
qARV=qZY(1,1);
#
if (Aq1-eps1) > 0
switch (qARV)
case{1}
qz1=qkoormus(1,1);
qx1=qkoormus(1,2);
Li1=qkoormus(1,3);
#Li=qkoormus(1,4); # aqL -kui koormus lõppeb enne lõppu
#Li=qkoormus(1,5); # varda lõpp (vaja kohendada)
#if (xx-aLi1) >= eps
#xarv=xx-aqA;
#disp(' olen krms... rida176  ');
vB1=yzhqzm(baasi0,xx,Li1,qx1,qz1,EA,EI);
vB=vB1;
#disp(' olen krms... rida179  ');
#else
#vB(1,1)=0.0;
#vB(2,1)=0.0;
#vB(3,1)=0.0;
#endif
#Zq0(1,1)=-qx*x;
#Zq0(2,1)=-qz*x;
#Zq0(3,1)=-qz*x^2/2;
#  vaja lisada juht kui koormus lõppeb enne lõppu
#
case{2}
qz1=qkoormus(1,1);
qx1=qkoormus(1,2);
Li1=qkoormus(1,3);
vB1=yzhqzm(baasi0,xx,Li1,qx1,qz1,EA,EI);
qz2=qkoormus(2,1);
qx2=qkoormus(2,2);
Li2=qkoormus(2,3);
vB2=yzhqzm(baasi0,xx,Li2,qx2,qz2,EA,EI);
vB=vB1+vB2;
case{3}
qz1=qkoormus(1,1);
qx1=qkoormus(1,2);
Li1=qkoormus(1,3);
vB1=yzhqzm(baasi0,xx,Li1,qx1,qz1,EA,EI);
qz2=qkoormus(2,1);
qx2=qkoormus(2,2);
Li2=qkoormus(2,3);
vB2=yzhqzm(baasi0,xx,Li2,qx2,qz2,EA,EI);
qz3=qkoormus(3,1);
qx3=qkoormus(3,2);
Li3=qkoormus(3,3);
vB3=yzhqzm(baasi0,xx,Li3,qx3,qz3,EA,EI);
vB=vB1+vB2+vB3;
case{4}
qz1=qkoormus(1,1);
qx1=qkoormus(1,2);
Li1=qkoormus(1,3);
vB1=yzhqzm(baasi0,xx,Li1,qx1,qz1,EA,EI);
qz2=qkoormus(2,1);
qx2=qkoormus(2,2);
Li2=qkoormus(2,3);
vB2=yzhqzm(baasi0,xx,Li2,qx2,qz2,EA,EI);
qz3=qkoormus(3,1);
qx3=qkoormus(3,2);
Li3=qkoormus(3,3);
vB3=yzhqzm(baasi0,xx,Li3,qx3,qz3,EA,EI);
qz4=qkoormus(4,1);
qx4=qkoormus(4,2);
Li4=qkoormus(4,3);
vB4=yzhqzm(baasi0,xx,Li4,qx4,qz4,EA,EI);
vB=vB1+vB2+vB3+vB4;
#%?otherwise
endswitch
else
# disp(' Jaotatud koormus puudub  '); 
vB=[0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
#vB=vB;
endif
#%%%%%%%%%%%%%%
#%vB=yzShqz(x,qx,qz)  #vB=yzhqz(baasi0,Li,qx,qz,EA,EI);    
#%vFz=yzSfzv(xx,aLx,Fx,Fz) #yzTfzv(baasi0,xx,aLx,Fz,EI); #yzhqzm(baasi0,xx,Li,qx,qz,EA,EI);
#vB=vB+vFz;
#%%%%%%%%%%%%%%
vBs=vB+vFz;
Sj=vBs;
#
endfunction
