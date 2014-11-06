## SisejoudPunktism.m
function SF=SisejoudPunktism(VardaNr,X,AlgPar,lvarras,selem,esQkoormus,esFjoud,suurused)
## The forces of frame element, ElementNr, at x = X
##=========================================================================
## This Program is written by Andres Lahe,   2010-07-16
##                    e-mail: andres.lahe@ttu.ee
## LAST MODIFIED: Andres Lahe,   2010-08-03
## Copyright (c)   2010 by Tallinn University of Technology 
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
##===========================================================================
#
#if nargin != 8
if ~(nargin==8)
error(' function SisejoudPunktis have wrong number of input arguments!') 
endif
#%%%%%%%%%%%%%%%%%
i=VardaNr;
xx=X;
baasi0=1.0;
%%%%%%%%%%
%vF=zeros(6,12);
    EI=selem(i,13);   %  from topology
    EA=selem(i,14);   %    "      "
    GAr=selem(i,15);  %    "      "
Li=lvarras(i,1);

#%qx=qxZ(i,1);
#%qz=qzZ(i,1);
#%aLx=aLXx(i,1);
#%Fz=FZz(i,1);
#%Fx=FZx(i,1);

#Fjoud=esFjoud(:,1:2,i);
Fjoud=esFjoud(:,1:3,i);
#%       Fz    Fx     a
#%Fjoud=[0.0    0.0   0.0;
#qkoormus=esQkoormus(:,1:3,i); 
qkoormus=esQkoormus(:,1:4,i);
#%       qz     qx     qA      qL
#%qkoormus=[0.0   0.0   0.0  0.0   0.0]; 

AP=AlgPar(i,:)';
#---------The transfer matrix equation --------
vvF=ylfhlin(1.0,xx,EA,GAr,EI);
#%vvB=yzhqz(1.0,xx,qx,qz,EA,EI);
#%vvFz=yzfzv(1.0,xx,aLx,Fx,Fz,EA,EI);
Sj=ESTFrKrmus(baasi0,xx,Li,Fjoud,qkoormus,EA,EI); 
 vB=Sj;
        Fvv(:,1)=vvF*AP+vB;

#%Fvv(:,1)=vvF*AP+vvB+vvFz;
##  %9.2f %9.2e
disp(sprintf('%15s %2i %9s  %4.7f   ', 'Forces of element',VardaNr,' at x =', X))
%
for i=1:3 
disp(sprintf('%14s %9.5e   ',suurused(i,:), Fvv(i,1)))
endfor
%
for i=4:6 
disp(sprintf('%14s %9.5f   ',suurused(i,:), Fvv(i,1)))
endfor
%
%disp('------------------')
%
SF=Fvv(:,1);
%
endfunction
