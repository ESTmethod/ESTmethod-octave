## spESTframe1LaheDefWFI.m
close all
clear
##
#dbwhere(175)
##
disp('==================================================================')
disp(' spESTframe1LaheDefWFI.m   ')
disp(' Analysis of a plane frame using the EST method. ')
disp(' The EST method finds the displacements and forces. ')
disp('---  ')
disp(' For problem: http://digi.lib.ttu.ee/opik_eme/ylesanded.pdf#page=39 ')
disp('---  ')
# For problem: http://digi.lib.ttu.ee/opik_eme/ylesanded.pdf#page=39
disp(' The equations matrix are non-symmetric sparse matrix. ')
## http://digi.lib.ttu.ee/opik_eme/slaidid/raamEST1en.pdf
##disp('  ')
#disp('---  ')
##disp(' http://digi.lib.ttu.ee/opik_eme/ylesanded.pdf#page=39 ')
#disp(' http://digi.lib.ttu.ee/opik_eme/Ehitusmehaanika.pdf#page=419 ')
## GNU octave version >= 3.0.x
##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-07-16
## Last modified: Andres Lahe, 2013-07-21
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
##siireVardaA
 eps1=0.000001;
 epsgraf=0.01;
 global spA
  global X
 global siireVardaA
disp('---  ')
disp(' daf -- the displacements and the forces at the ends of element. ')
disp(' The element of plane frame have 12 daf''s (6+6). ')
disp(' The number of unknowns are 12*number_of_elements + number_of_support_reactions. ')
disp(' The matrices are assembled into compressed column sparse matrices. ')
disp('---  ')
disp('  ')
#---------  Number of frame nodes ----------
Number_of_frame_nodes=8  # 6
                            SolmedeArv=Number_of_frame_nodes;
#---------  Number of elements ----------
Number_of_elements=7   #  5
                            ElementideArv=Number_of_elements;
#---------  Number of supports ----------
Number_of_support_reactions=8
                            Ntoerkts=Number_of_support_reactions;
disp('  ')

spNNK=12*ElementideArv+Ntoerkts;
#==========
   The_number_of_unknowns_are=spNNK
disp('---  ')
#==========
disp('--------- Element properties ------------')
#E=2.1E+11  # Elastsusmoodul (E=210GPa) E=2.1E+11 [Pa]
#E=2.1E+05 # kPa
#% I-40  Hristloige=400mm=0.4m   I=1902cm4=1.902E+03cm4=1.902E-05m4
#
EIp=20000
#EIp=4.6*10^15
#EIr=4.6*10^13
#EIr=40000
#EAp=4.6*10^6;
EAp=4.6*10^15;
#EAr=6.8*10^6;
EAr=6.8*10^15;
GAp=0.4*EAp;
#GAp=0.4*4.6*10^6;
GAr=0.4*EAr;
disp('--------- Scaling multiplier for displacements----------')
#baasi0=EIp/4.0  # scaling multiplier for displacements
#%baasi0=1.0;
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' Displacements and forces calculated on parts of the element ''Nmitmeks''.  ')
Nmitmeks=4
Lp=8.0; # graphics axis
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##Fjoud=esFjoud(:,1:2,i);
#Fjoud=esFjoud(:,1:3,i);
##qkoormus=esQkoormus(:,1:3,i);
#qkoormus=esQkoormus(:,1:4,i);
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## ---- load variants -----
load_variant=1;
#load_variant=2;
#load_variant=3;
#load_variant=4;
#load_variant=5;
#load_variant=6;
#load_variant=7;
#load_variant=8;
#load_variant=9;
#load_variant=10;
        koormusvariant=load_variant;
###
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (koormusvariant)
case{1}
#
disp(' Load variant 1 ')
#
disp(' Initial data  ')
l=6.0
h=3.0
Isuhe=2  % Isuhe=I1/I2
p1=12.0
p2=0.0
p3=0.0
F1=50.0
F2=0.0
F3=0.0
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr=Isuhe*EIp
#
baasi0=EIp/h
xi=0.6;
l1=l;
l2=l;
h08=0.8*h
deltah=(1.0*0.2*h)/l1;
#% sqrt(h^2+l1^2);
L1=h08;
##L2=l1;
L2=sqrt(deltah^2+1.0^2); # length of the element
L3=sqrt((0.2*h)^2+l1^2); # length of the element
L4=h;
L5=sqrt((0.2*h)^2+l2^2); # length of the element
L6=L1;
L7=L2;
#sqrt(h^2+l4^2);
Lp=l1+l2;   # graphics axis
sinA3=(0.2*h)/L3;
cosA3=l1/L3;
sinA2=sinA3;
cosA2=cosA3;
sinA5=(0.2*h)/L5;
cosA5=l2/L5;
#
h1=xi*h;
h2=(1-xi)*h;
d=h1;
H=h;
#
q1=p3;
#
  if (p1 > 0)
    p12=p1
  else
    p12=p2
   endif
#
q2=p12
q3=p12
q4=0.0;
q5=0.0;
q6=0.0;
q7=0.0;
#
qz3v=q3*l1/L3;   # load/length of the element
qz3=qz3v*cosA3;  # projection onto z-axis
qx3=-qz3v*sinA3; # projection onto x-axis
qz2v=q2*1.0/L2;   # load/length of the element
qz2=-qz2v*cosA2;  # projection onto z-axis
qx2=qz2v*sinA2; # projection onto x-axis


qz2=-qz3;
qx2=-qx3;
#
Lpunkt1=0.4*h;
Lpunkt2=0.3*h;
Lpunkt3=0.4*L5;
#
#Fz1=0.0;
#Fz2=0.0;
Fz4=-F1;
aF4=Lpunkt1;
##Fx3=0.0;
Fz5=F3*cosA5;   % Fz2 projection onto z-axis
Fx5=F3*sinA5;  % Fx2 projection onto x-axis
aF5=Lpunkt3;
#
Fz6=F2;
aF6=Lpunkt2;
#Fz7=0.0;
#
#disp(' Element load  in local coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN/m^2
LoadsqONelement=4;
esQkoormu=zeros(LoadsqONelement,4,ElementideArv);
esQkoormus(1,1:4,1)=[q1 0.0 0.0 L1];
esQkoormus(1,1:4,2)=[qz2 qx2 0.0 L2];
esQkoormus(1,1:4,3)=[qz3 qx3 0.0 L3];
esQkoormus(1,1:4,4)=[0.0 0.0 0.0 L4];
esQkoormus(1,1:4,5)=[0.0 0.0 0.0 L5];
esQkoormus(1,1:4,6)=[0.0 0.0 0.0 L6];
esQkoormus(1,1:4,7)=[0.0 0.0 0.0 L7];
#esQkoormus
#
# Point load in local coordinate z and x direction  kN
#  Fz, Fx, aF - coordinate of the point of force application
LoadsF_on_Element=5;
esFjoud=zeros(LoadsF_on_Element,2,ElementideArv);
esFjoud(1,1:3,1)=[0.0 0.0 L1];
esFjoud(1,1:3,2)=[0.0 0.0 L2];
esFjoud(1,1:3,3)=[0.0 0.0 L3];
esFjoud(1,1:3,4)=[Fz4 0.0 aF4];
esFjoud(1,1:3,5)=[Fz5 Fx5 aF5];
esFjoud(1,1:3,6)=[Fz6 0.0 aF6];
esFjoud(1,1:3,7)=[0.0 0.0 L7];
#esFjoud

#disp(' Node forces in global  coordinates ')
# sSolmF(forces,1,nodes); forces=[Fx; Fz; My]
sSolmF = zeros(3,1,SolmedeArv);
#sSolmF = zeros(3,1,SolmedeArv);
#sSolmF(:,1,1)= 0.0
#sSolmF(:,1,2)= 0.0
#sSolmF(:,1,3)= 0.0
#sSolmF(:,1,4)= 0.0
#sSolmF(:,1,5)= 0.0

#%%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3,1,SolmedeArv);
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
##%%%%
case{2}
disp(' Load variant 2 ')
#
disp(' Initial data  ')
l=8.0
h=4.0
Isuhe=2  % Isuhe=I1/I2
p1=0.0
p2=14.0
p3=0.0
F1=0.0
F2=50.0
F3=0.0
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr=Isuhe*EIp
#
baasi0=EIp/h
xi=0.6;
l1=l;
l2=l;
h08=0.8*h
deltah=(1.0*0.2*h)/l1;
#% sqrt(h^2+l1^2);
L1=h08;
##L2=l1;
L2=sqrt(deltah^2+1.0^2); # length of the element
L3=sqrt((0.2*h)^2+l1^2); # length of the element
L4=h;
L5=sqrt((0.2*h)^2+l2^2); # length of the element
L6=L1;
L7=L2;
#sqrt(h^2+l4^2);
Lp=l1+l2;   # graphics axis
sinA3=(0.2*h)/L3;
cosA3=l1/L3;
sinA2=sinA3;
cosA2=cosA3;
sinA5=(0.2*h)/L5;
cosA5=l2/L5;
#
h1=xi*h;
h2=(1-xi)*h;
d=h1;
H=h;
#
q1=p3;
#
  if (p1 > 0)
    p12=p1
  else
    p12=p2
   endif
#
q2=p12
q3=p12
q4=0.0;
q5=0.0;
q6=0.0;
q7=0.0;
#
qz3v=q3*l1/L3;   # load/length of the element
qz3=qz3v*cosA3;  # projection onto z-axis
qx3=-qz3v*sinA3; # projection onto x-axis
qz2v=q2*1.0/L2;   # load/length of the element
qz2=-qz2v*cosA2;  # projection onto z-axis
qx2=qz2v*sinA2; # projection onto x-axis


qz2=-qz3;
qx2=-qx3;
#
Lpunkt1=0.4*h;
Lpunkt2=0.3*h;
Lpunkt3=0.4*L5;
#
#Fz1=0.0;
#Fz2=0.0;
Fz4=-F1;
aF4=Lpunkt1;
##Fx3=0.0;
Fz5=F3*cosA5;   % Fz2 projection onto z-axis
Fx5=F3*sinA5;  % Fx2 projection onto x-axis
aF5=Lpunkt3;
#
Fz6=F2;
aF6=Lpunkt2;
#Fz7=0.0;
#
#disp(' Element load  in local coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN/m^2
LoadsqONelement=4;
esQkoormu=zeros(LoadsqONelement,4,ElementideArv);
esQkoormus(1,1:4,1)=[q1 0.0 0.0 L1];
esQkoormus(1,1:4,2)=[qz2 qx2 0.0 L2];
esQkoormus(1,1:4,3)=[qz3 qx3 0.0 L3];
esQkoormus(1,1:4,4)=[0.0 0.0 0.0 L4];
esQkoormus(1,1:4,5)=[0.0 0.0 0.0 L5];
esQkoormus(1,1:4,6)=[0.0 0.0 0.0 L6];
esQkoormus(1,1:4,7)=[0.0 0.0 0.0 L7];
#esQkoormus
#
# Point load in local coordinate z and x direction  kN
#  Fz, Fx, aF - coordinate of the point of force application
LoadsF_on_Element=5;
esFjoud=zeros(LoadsF_on_Element,2,ElementideArv);
esFjoud(1,1:3,1)=[0.0 0.0 L1];
esFjoud(1,1:3,2)=[0.0 0.0 L2];
esFjoud(1,1:3,3)=[0.0 0.0 L3];
esFjoud(1,1:3,4)=[Fz4 0.0 aF4];
esFjoud(1,1:3,5)=[Fz5 Fx5 aF5];
esFjoud(1,1:3,6)=[Fz6 0.0 aF6];
esFjoud(1,1:3,7)=[0.0 0.0 L7];
#esFjoud

#disp(' Node forces in global  coordinates ')
# sSolmF(forces,1,nodes); forces=[Fx; Fz; My]
sSolmF = zeros(3,1,SolmedeArv);
#sSolmF = zeros(3,1,SolmedeArv);
#sSolmF(:,1,1)= 0.0
#sSolmF(:,1,2)= 0.0
#sSolmF(:,1,3)= 0.0
#sSolmF(:,1,4)= 0.0
#sSolmF(:,1,5)= 0.0

#%%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3,1,SolmedeArv);
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
#%%%%%
case{3}
disp(' Load variant 3 ')
#
disp(' Initial data  ')
l=10.0
h=4.0
Isuhe=2  % Isuhe=I1/I2
p1=0.0
p2=0.0
p3=8.0
F1=0.0
F2=0.0
F3=50.0
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr=Isuhe*EIp
#
baasi0=EIp/h
xi=0.6;
l1=l;
l2=l;
h08=0.8*h
deltah=(1.0*0.2*h)/l1;
#% sqrt(h^2+l1^2);
L1=h08;
##L2=l1;
L2=sqrt(deltah^2+1.0^2); # length of the element
L3=sqrt((0.2*h)^2+l1^2); # length of the element
L4=h;
L5=sqrt((0.2*h)^2+l2^2); # length of the element
L6=L1;
L7=L2;
#sqrt(h^2+l4^2);
Lp=l1+l2;   # graphics axis
sinA3=(0.2*h)/L3;
cosA3=l1/L3;
sinA2=sinA3;
cosA2=cosA3;
sinA5=(0.2*h)/L5;
cosA5=l2/L5;
#
h1=xi*h;
h2=(1-xi)*h;
d=h1;
H=h;
#
q1=p3;
#
  if (p1 > 0)
    p12=p1
  else
    p12=p2
   endif
#
q2=p12
q3=p12
q4=0.0;
q5=0.0;
q6=0.0;
q7=0.0;
#
qz3v=q3*l1/L3;   # load/length of the element
qz3=qz3v*cosA3;  # projection onto z-axis
qx3=-qz3v*sinA3; # projection onto x-axis
qz2v=q2*1.0/L2;   # load/length of the element
qz2=-qz2v*cosA2;  # projection onto z-axis
qx2=qz2v*sinA2; # projection onto x-axis


qz2=-qz3;
qx2=-qx3;
#
Lpunkt1=0.4*h;
Lpunkt2=0.3*h;
Lpunkt3=0.4*L5;
#
#Fz1=0.0;
#Fz2=0.0;
Fz4=-F1;
aF4=Lpunkt1;
##Fx3=0.0;
Fz5=F3*cosA5;   % Fz2 projection onto z-axis
Fx5=F3*sinA5;  % Fx2 projection onto x-axis
aF5=Lpunkt3;
#
Fz6=F2;
aF6=Lpunkt2;
#Fz7=0.0;
#
#disp(' Element load  in local coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN/m^2
LoadsqONelement=4;
esQkoormu=zeros(LoadsqONelement,4,ElementideArv);
esQkoormus(1,1:4,1)=[q1 0.0 0.0 L1];
esQkoormus(1,1:4,2)=[qz2 qx2 0.0 L2];
esQkoormus(1,1:4,3)=[qz3 qx3 0.0 L3];
esQkoormus(1,1:4,4)=[0.0 0.0 0.0 L4];
esQkoormus(1,1:4,5)=[0.0 0.0 0.0 L5];
esQkoormus(1,1:4,6)=[0.0 0.0 0.0 L6];
esQkoormus(1,1:4,7)=[0.0 0.0 0.0 L7];
#esQkoormus
#
# Point load in local coordinate z and x direction  kN
#  Fz, Fx, aF - coordinate of the point of force application
LoadsF_on_Element=5;
esFjoud=zeros(LoadsF_on_Element,2,ElementideArv);
esFjoud(1,1:3,1)=[0.0 0.0 L1];
esFjoud(1,1:3,2)=[0.0 0.0 L2];
esFjoud(1,1:3,3)=[0.0 0.0 L3];
esFjoud(1,1:3,4)=[Fz4 0.0 aF4];
esFjoud(1,1:3,5)=[Fz5 Fx5 aF5];
esFjoud(1,1:3,6)=[Fz6 0.0 aF6];
esFjoud(1,1:3,7)=[0.0 0.0 L7];
#esFjoud

#disp(' Node forces in global  coordinates ')
# sSolmF(forces,1,nodes); forces=[Fx; Fz; My]
sSolmF = zeros(3,1,SolmedeArv);
#sSolmF = zeros(3,1,SolmedeArv);
#sSolmF(:,1,1)= 0.0
#sSolmF(:,1,2)= 0.0
#sSolmF(:,1,3)= 0.0
#sSolmF(:,1,4)= 0.0
#sSolmF(:,1,5)= 0.0

#%%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3,1,SolmedeArv);
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
#%%%%%
case{4}
disp(' Load variant 4 ')
#
disp(' Initial data  ')
l=6.0
h=4.0
Isuhe=2  % Isuhe=I1/I2
p1=14.0
p2=0.0
p3=0.0
F1=40.0
F2=0.0
F3=0.0
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr=Isuhe*EIp
#
baasi0=EIp/h
xi=0.6;
l1=l;
l2=l;
h08=0.8*h
deltah=(1.0*0.2*h)/l1;
#% sqrt(h^2+l1^2);
L1=h08;
##L2=l1;
L2=sqrt(deltah^2+1.0^2); # length of the element
L3=sqrt((0.2*h)^2+l1^2); # length of the element
L4=h;
L5=sqrt((0.2*h)^2+l2^2); # length of the element
L6=L1;
L7=L2;
#sqrt(h^2+l4^2);
Lp=l1+l2;   # graphics axis
sinA3=(0.2*h)/L3;
cosA3=l1/L3;
sinA2=sinA3;
cosA2=cosA3;
sinA5=(0.2*h)/L5;
cosA5=l2/L5;
#
h1=xi*h;
h2=(1-xi)*h;
d=h1;
H=h;
#
q1=p3;
#
  if (p1 > 0)
    p12=p1
  else
    p12=p2
   endif
#
q2=p12
q3=p12
q4=0.0;
q5=0.0;
q6=0.0;
q7=0.0;
#
qz3v=q3*l1/L3;   # load/length of the element
qz3=qz3v*cosA3;  # projection onto z-axis
qx3=-qz3v*sinA3; # projection onto x-axis
qz2v=q2*1.0/L2;   # load/length of the element
qz2=-qz2v*cosA2;  # projection onto z-axis
qx2=qz2v*sinA2; # projection onto x-axis


qz2=-qz3;
qx2=-qx3;
#
Lpunkt1=0.4*h;
Lpunkt2=0.3*h;
Lpunkt3=0.4*L5;
#
#Fz1=0.0;
#Fz2=0.0;
Fz4=-F1;
aF4=Lpunkt1;
##Fx3=0.0;
Fz5=F3*cosA5;   % Fz2 projection onto z-axis
Fx5=F3*sinA5;  % Fx2 projection onto x-axis
aF5=Lpunkt3;
#
Fz6=F2;
aF6=Lpunkt2;
#Fz7=0.0;
#
#disp(' Element load  in local coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN/m^2
LoadsqONelement=4;
esQkoormu=zeros(LoadsqONelement,4,ElementideArv);
esQkoormus(1,1:4,1)=[q1 0.0 0.0 L1];
esQkoormus(1,1:4,2)=[qz2 qx2 0.0 L2];
esQkoormus(1,1:4,3)=[qz3 qx3 0.0 L3];
esQkoormus(1,1:4,4)=[0.0 0.0 0.0 L4];
esQkoormus(1,1:4,5)=[0.0 0.0 0.0 L5];
esQkoormus(1,1:4,6)=[0.0 0.0 0.0 L6];
esQkoormus(1,1:4,7)=[0.0 0.0 0.0 L7];
#esQkoormus
#
# Point load in local coordinate z and x direction  kN
#  Fz, Fx, aF - coordinate of the point of force application
LoadsF_on_Element=5;
esFjoud=zeros(LoadsF_on_Element,2,ElementideArv);
esFjoud(1,1:3,1)=[0.0 0.0 L1];
esFjoud(1,1:3,2)=[0.0 0.0 L2];
esFjoud(1,1:3,3)=[0.0 0.0 L3];
esFjoud(1,1:3,4)=[Fz4 0.0 aF4];
esFjoud(1,1:3,5)=[Fz5 Fx5 aF5];
esFjoud(1,1:3,6)=[Fz6 0.0 aF6];
esFjoud(1,1:3,7)=[0.0 0.0 L7];
#esFjoud

#disp(' Node forces in global  coordinates ')
# sSolmF(forces,1,nodes); forces=[Fx; Fz; My]
sSolmF = zeros(3,1,SolmedeArv);
#sSolmF = zeros(3,1,SolmedeArv);
#sSolmF(:,1,1)= 0.0
#sSolmF(:,1,2)= 0.0
#sSolmF(:,1,3)= 0.0
#sSolmF(:,1,4)= 0.0
#sSolmF(:,1,5)= 0.0

#%%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3,1,SolmedeArv);
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
#%%%%%
case{5}
disp(' Load variant 5 ')
#
disp(' Initial data  ')
l=8.0
h=5.0
Isuhe=2  % Isuhe=I1/I2
p1=0.0
p2=8.0
p3=0.0
F1=0.0
F2=40.0
F3=0.0
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr=Isuhe*EIp
#
baasi0=EIp/h
xi=0.6;
l1=l;
l2=l;
h08=0.8*h
deltah=(1.0*0.2*h)/l1;
#% sqrt(h^2+l1^2);
L1=h08;
##L2=l1;
L2=sqrt(deltah^2+1.0^2); # length of the element
L3=sqrt((0.2*h)^2+l1^2); # length of the element
L4=h;
L5=sqrt((0.2*h)^2+l2^2); # length of the element
L6=L1;
L7=L2;
#sqrt(h^2+l4^2);
Lp=l1+l2;   # graphics axis
sinA3=(0.2*h)/L3;
cosA3=l1/L3;
sinA2=sinA3;
cosA2=cosA3;
sinA5=(0.2*h)/L5;
cosA5=l2/L5;
#
h1=xi*h;
h2=(1-xi)*h;
d=h1;
H=h;
#
q1=p3;
#
  if (p1 > 0)
    p12=p1
  else
    p12=p2
   endif
#
q2=p12
q3=p12
q4=0.0;
q5=0.0;
q6=0.0;
q7=0.0;
#
qz3v=q3*l1/L3;   # load/length of the element
qz3=qz3v*cosA3;  # projection onto z-axis
qx3=-qz3v*sinA3; # projection onto x-axis
qz2v=q2*1.0/L2;   # load/length of the element
qz2=-qz2v*cosA2;  # projection onto z-axis
qx2=qz2v*sinA2; # projection onto x-axis


qz2=-qz3;
qx2=-qx3;
#
Lpunkt1=0.4*h;
Lpunkt2=0.3*h;
Lpunkt3=0.4*L5;
#
#Fz1=0.0;
#Fz2=0.0;
Fz4=-F1;
aF4=Lpunkt1;
##Fx3=0.0;
Fz5=F3*cosA5;   % Fz2 projection onto z-axis
Fx5=F3*sinA5;  % Fx2 projection onto x-axis
aF5=Lpunkt3;
#
Fz6=F2;
aF6=Lpunkt2;
#Fz7=0.0;
#
#disp(' Element load  in local coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN/m^2
LoadsqONelement=4;
esQkoormu=zeros(LoadsqONelement,4,ElementideArv);
esQkoormus(1,1:4,1)=[q1 0.0 0.0 L1];
esQkoormus(1,1:4,2)=[qz2 qx2 0.0 L2];
esQkoormus(1,1:4,3)=[qz3 qx3 0.0 L3];
esQkoormus(1,1:4,4)=[0.0 0.0 0.0 L4];
esQkoormus(1,1:4,5)=[0.0 0.0 0.0 L5];
esQkoormus(1,1:4,6)=[0.0 0.0 0.0 L6];
esQkoormus(1,1:4,7)=[0.0 0.0 0.0 L7];
#esQkoormus
#
# Point load in local coordinate z and x direction  kN
#  Fz, Fx, aF - coordinate of the point of force application
LoadsF_on_Element=5;
esFjoud=zeros(LoadsF_on_Element,2,ElementideArv);
esFjoud(1,1:3,1)=[0.0 0.0 L1];
esFjoud(1,1:3,2)=[0.0 0.0 L2];
esFjoud(1,1:3,3)=[0.0 0.0 L3];
esFjoud(1,1:3,4)=[Fz4 0.0 aF4];
esFjoud(1,1:3,5)=[Fz5 Fx5 aF5];
esFjoud(1,1:3,6)=[Fz6 0.0 aF6];
esFjoud(1,1:3,7)=[0.0 0.0 L7];
#esFjoud

#disp(' Node forces in global  coordinates ')
# sSolmF(forces,1,nodes); forces=[Fx; Fz; My]
sSolmF = zeros(3,1,SolmedeArv);
#sSolmF = zeros(3,1,SolmedeArv);
#sSolmF(:,1,1)= 0.0
#sSolmF(:,1,2)= 0.0
#sSolmF(:,1,3)= 0.0
#sSolmF(:,1,4)= 0.0
#sSolmF(:,1,5)= 0.0

#%%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3,1,SolmedeArv);
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
#%%%%%
case{6}
disp(' Load variant 6 ')
#
disp(' Initial data  ')
l=10.0
h=5.0
Isuhe=3  % Isuhe=I1/I2
p1=0.0
p2=0.0
p3=12.0
F1=0.0
F2=0.0
F3=40.0
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr=Isuhe*EIp
#
baasi0=EIp/h
xi=0.6;
l1=l;
l2=l;
h08=0.8*h
deltah=(1.0*0.2*h)/l1;
#% sqrt(h^2+l1^2);
L1=h08;
##L2=l1;
L2=sqrt(deltah^2+1.0^2); # length of the element
L3=sqrt((0.2*h)^2+l1^2); # length of the element
L4=h;
L5=sqrt((0.2*h)^2+l2^2); # length of the element
L6=L1;
L7=L2;
#sqrt(h^2+l4^2);
Lp=l1+l2;   # graphics axis
sinA3=(0.2*h)/L3;
cosA3=l1/L3;
sinA2=sinA3;
cosA2=cosA3;
sinA5=(0.2*h)/L5;
cosA5=l2/L5;
#
h1=xi*h;
h2=(1-xi)*h;
d=h1;
H=h;
#
q1=p3;
#
  if (p1 > 0)
    p12=p1
  else
    p12=p2
   endif
#
q2=p12
q3=p12
q4=0.0;
q5=0.0;
q6=0.0;
q7=0.0;
#
qz3v=q3*l1/L3;   # load/length of the element
qz3=qz3v*cosA3;  # projection onto z-axis
qx3=-qz3v*sinA3; # projection onto x-axis
qz2v=q2*1.0/L2;   # load/length of the element
qz2=-qz2v*cosA2;  # projection onto z-axis
qx2=qz2v*sinA2; # projection onto x-axis


qz2=-qz3;
qx2=-qx3;
#
Lpunkt1=0.4*h;
Lpunkt2=0.3*h;
Lpunkt3=0.4*L5;
#
#Fz1=0.0;
#Fz2=0.0;
Fz4=-F1;
aF4=Lpunkt1;
##Fx3=0.0;
Fz5=F3*cosA5;   % Fz2 projection onto z-axis
Fx5=F3*sinA5;  % Fx2 projection onto x-axis
aF5=Lpunkt3;
#
Fz6=F2;
aF6=Lpunkt2;
#Fz7=0.0;
#
#disp(' Element load  in local coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN/m^2
LoadsqONelement=4;
esQkoormu=zeros(LoadsqONelement,4,ElementideArv);
esQkoormus(1,1:4,1)=[q1 0.0 0.0 L1];
esQkoormus(1,1:4,2)=[qz2 qx2 0.0 L2];
esQkoormus(1,1:4,3)=[qz3 qx3 0.0 L3];
esQkoormus(1,1:4,4)=[0.0 0.0 0.0 L4];
esQkoormus(1,1:4,5)=[0.0 0.0 0.0 L5];
esQkoormus(1,1:4,6)=[0.0 0.0 0.0 L6];
esQkoormus(1,1:4,7)=[0.0 0.0 0.0 L7];
#esQkoormus
#
# Point load in local coordinate z and x direction  kN
#  Fz, Fx, aF - coordinate of the point of force application
LoadsF_on_Element=5;
esFjoud=zeros(LoadsF_on_Element,2,ElementideArv);
esFjoud(1,1:3,1)=[0.0 0.0 L1];
esFjoud(1,1:3,2)=[0.0 0.0 L2];
esFjoud(1,1:3,3)=[0.0 0.0 L3];
esFjoud(1,1:3,4)=[Fz4 0.0 aF4];
esFjoud(1,1:3,5)=[Fz5 Fx5 aF5];
esFjoud(1,1:3,6)=[Fz6 0.0 aF6];
esFjoud(1,1:3,7)=[0.0 0.0 L7];
#esFjoud

#disp(' Node forces in global  coordinates ')
# sSolmF(forces,1,nodes); forces=[Fx; Fz; My]
sSolmF = zeros(3,1,SolmedeArv);
#sSolmF = zeros(3,1,SolmedeArv);
#sSolmF(:,1,1)= 0.0
#sSolmF(:,1,2)= 0.0
#sSolmF(:,1,3)= 0.0
#sSolmF(:,1,4)= 0.0
#sSolmF(:,1,5)= 0.0

#%%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3,1,SolmedeArv);
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
#%%%%%
case{7}
disp(' Load variant 7 ')
#
disp(' Initial data  ')
l=6.0
h=3.0
Isuhe=3  % Isuhe=I1/I2
p1=8.0
p2=0.0
p3=0.0
F1=30.0
F2=0.0
F3=0.0
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr=Isuhe*EIp
#
baasi0=EIp/h
xi=0.6;
l1=l;
l2=l;
h08=0.8*h
deltah=(1.0*0.2*h)/l1;
#% sqrt(h^2+l1^2);
L1=h08;
##L2=l1;
L2=sqrt(deltah^2+1.0^2); # length of the element
L3=sqrt((0.2*h)^2+l1^2); # length of the element
L4=h;
L5=sqrt((0.2*h)^2+l2^2); # length of the element
L6=L1;
L7=L2;
#sqrt(h^2+l4^2);
Lp=l1+l2;   # graphics axis
sinA3=(0.2*h)/L3;
cosA3=l1/L3;
sinA2=sinA3;
cosA2=cosA3;
sinA5=(0.2*h)/L5;
cosA5=l2/L5;
#
h1=xi*h;
h2=(1-xi)*h;
d=h1;
H=h;
#
q1=p3;
#
  if (p1 > 0)
    p12=p1
  else
    p12=p2
   endif
#
q2=p12
q3=p12
q4=0.0;
q5=0.0;
q6=0.0;
q7=0.0;
#
qz3v=q3*l1/L3;   # load/length of the element
qz3=qz3v*cosA3;  # projection onto z-axis
qx3=-qz3v*sinA3; # projection onto x-axis
qz2v=q2*1.0/L2;   # load/length of the element
qz2=-qz2v*cosA2;  # projection onto z-axis
qx2=qz2v*sinA2; # projection onto x-axis


qz2=-qz3;
qx2=-qx3;
#
Lpunkt1=0.4*h;
Lpunkt2=0.3*h;
Lpunkt3=0.4*L5;
#
#Fz1=0.0;
#Fz2=0.0;
Fz4=-F1;
aF4=Lpunkt1;
##Fx3=0.0;
Fz5=F3*cosA5;   % Fz2 projection onto z-axis
Fx5=F3*sinA5;  % Fx2 projection onto x-axis
aF5=Lpunkt3;
#
Fz6=F2;
aF6=Lpunkt2;
#Fz7=0.0;
#
#disp(' Element load  in local coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN/m^2
LoadsqONelement=4;
esQkoormu=zeros(LoadsqONelement,4,ElementideArv);
esQkoormus(1,1:4,1)=[q1 0.0 0.0 L1];
esQkoormus(1,1:4,2)=[qz2 qx2 0.0 L2];
esQkoormus(1,1:4,3)=[qz3 qx3 0.0 L3];
esQkoormus(1,1:4,4)=[0.0 0.0 0.0 L4];
esQkoormus(1,1:4,5)=[0.0 0.0 0.0 L5];
esQkoormus(1,1:4,6)=[0.0 0.0 0.0 L6];
esQkoormus(1,1:4,7)=[0.0 0.0 0.0 L7];
#esQkoormus
#
# Point load in local coordinate z and x direction  kN
#  Fz, Fx, aF - coordinate of the point of force application
LoadsF_on_Element=5;
esFjoud=zeros(LoadsF_on_Element,2,ElementideArv);
esFjoud(1,1:3,1)=[0.0 0.0 L1];
esFjoud(1,1:3,2)=[0.0 0.0 L2];
esFjoud(1,1:3,3)=[0.0 0.0 L3];
esFjoud(1,1:3,4)=[Fz4 0.0 aF4];
esFjoud(1,1:3,5)=[Fz5 Fx5 aF5];
esFjoud(1,1:3,6)=[Fz6 0.0 aF6];
esFjoud(1,1:3,7)=[0.0 0.0 L7];
#esFjoud

#disp(' Node forces in global  coordinates ')
# sSolmF(forces,1,nodes); forces=[Fx; Fz; My]
sSolmF = zeros(3,1,SolmedeArv);
#sSolmF = zeros(3,1,SolmedeArv);
#sSolmF(:,1,1)= 0.0
#sSolmF(:,1,2)= 0.0
#sSolmF(:,1,3)= 0.0
#sSolmF(:,1,4)= 0.0
#sSolmF(:,1,5)= 0.0

#%%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3,1,SolmedeArv);
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
#%%%%%
case{8}
disp(' Load variant 8 ')
#
disp(' Initial data  ')
l=8.0
h=4.0
Isuhe=3  % Isuhe=I1/I2
p1=0.0
p2=16.0
p3=0.0
F1=0.0
F2=30.0
F3=0.0
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr=Isuhe*EIp
#
baasi0=EIp/h
xi=0.6;
l1=l;
l2=l;
h08=0.8*h
deltah=(1.0*0.2*h)/l1;
#% sqrt(h^2+l1^2);
L1=h08;
##L2=l1;
L2=sqrt(deltah^2+1.0^2); # length of the element
L3=sqrt((0.2*h)^2+l1^2); # length of the element
L4=h;
L5=sqrt((0.2*h)^2+l2^2); # length of the element
L6=L1;
L7=L2;
#sqrt(h^2+l4^2);
Lp=l1+l2;   # graphics axis
sinA3=(0.2*h)/L3;
cosA3=l1/L3;
sinA2=sinA3;
cosA2=cosA3;
sinA5=(0.2*h)/L5;
cosA5=l2/L5;
#
h1=xi*h;
h2=(1-xi)*h;
d=h1;
H=h;
#
q1=p3;
#
  if (p1 > 0)
    p12=p1
  else
    p12=p2
   endif
#
q2=p12
q3=p12
q4=0.0;
q5=0.0;
q6=0.0;
q7=0.0;
#
qz3v=q3*l1/L3;   # load/length of the element
qz3=qz3v*cosA3;  # projection onto z-axis
qx3=-qz3v*sinA3; # projection onto x-axis
qz2v=q2*1.0/L2;   # load/length of the element
qz2=-qz2v*cosA2;  # projection onto z-axis
qx2=qz2v*sinA2; # projection onto x-axis


qz2=-qz3;
qx2=-qx3;
#
Lpunkt1=0.4*h;
Lpunkt2=0.3*h;
Lpunkt3=0.4*L5;
#
#Fz1=0.0;
#Fz2=0.0;
Fz4=-F1;
aF4=Lpunkt1;
##Fx3=0.0;
Fz5=F3*cosA5;   % Fz2 projection onto z-axis
Fx5=F3*sinA5;  % Fx2 projection onto x-axis
aF5=Lpunkt3;
#
Fz6=F2;
aF6=Lpunkt2;
#Fz7=0.0;
#
#disp(' Element load  in local coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN/m^2
LoadsqONelement=4;
esQkoormu=zeros(LoadsqONelement,4,ElementideArv);
esQkoormus(1,1:4,1)=[q1 0.0 0.0 L1];
esQkoormus(1,1:4,2)=[qz2 qx2 0.0 L2];
esQkoormus(1,1:4,3)=[qz3 qx3 0.0 L3];
esQkoormus(1,1:4,4)=[0.0 0.0 0.0 L4];
esQkoormus(1,1:4,5)=[0.0 0.0 0.0 L5];
esQkoormus(1,1:4,6)=[0.0 0.0 0.0 L6];
esQkoormus(1,1:4,7)=[0.0 0.0 0.0 L7];
#esQkoormus
#
# Point load in local coordinate z and x direction  kN
#  Fz, Fx, aF - coordinate of the point of force application
LoadsF_on_Element=5;
esFjoud=zeros(LoadsF_on_Element,2,ElementideArv);
esFjoud(1,1:3,1)=[0.0 0.0 L1];
esFjoud(1,1:3,2)=[0.0 0.0 L2];
esFjoud(1,1:3,3)=[0.0 0.0 L3];
esFjoud(1,1:3,4)=[Fz4 0.0 aF4];
esFjoud(1,1:3,5)=[Fz5 Fx5 aF5];
esFjoud(1,1:3,6)=[Fz6 0.0 aF6];
esFjoud(1,1:3,7)=[0.0 0.0 L7];
#esFjoud

#disp(' Node forces in global  coordinates ')
# sSolmF(forces,1,nodes); forces=[Fx; Fz; My]
sSolmF = zeros(3,1,SolmedeArv);
#sSolmF = zeros(3,1,SolmedeArv);
#sSolmF(:,1,1)= 0.0
#sSolmF(:,1,2)= 0.0
#sSolmF(:,1,3)= 0.0
#sSolmF(:,1,4)= 0.0
#sSolmF(:,1,5)= 0.0

#%%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3,1,SolmedeArv);
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
#%%%%%
case{9}
disp(' Load variant 9 ')
#
disp(' Initial data  ')
l=10.0
h=4.0
Isuhe=3  % Isuhe=I1/I2
p1=0.0
p2=0.0
p3=14.0
F1=0.0
F2=0.0
F3=30.0
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr=Isuhe*EIp
#
baasi0=EIp/h
xi=0.6;
l1=l;
l2=l;
h08=0.8*h
deltah=(1.0*0.2*h)/l1;
#% sqrt(h^2+l1^2);
L1=h08;
##L2=l1;
L2=sqrt(deltah^2+1.0^2); # length of the element
L3=sqrt((0.2*h)^2+l1^2); # length of the element
L4=h;
L5=sqrt((0.2*h)^2+l2^2); # length of the element
L6=L1;
L7=L2;
#sqrt(h^2+l4^2);
Lp=l1+l2;   # graphics axis
sinA3=(0.2*h)/L3;
cosA3=l1/L3;
sinA2=sinA3;
cosA2=cosA3;
sinA5=(0.2*h)/L5;
cosA5=l2/L5;
#
h1=xi*h;
h2=(1-xi)*h;
d=h1;
H=h;
#
q1=p3;
#
  if (p1 > 0)
    p12=p1
  else
    p12=p2
   endif
#
q2=p12
q3=p12
q4=0.0;
q5=0.0;
q6=0.0;
q7=0.0;
#
qz3v=q3*l1/L3;   # load/length of the element
qz3=qz3v*cosA3;  # projection onto z-axis
qx3=-qz3v*sinA3; # projection onto x-axis
qz2v=q2*1.0/L2;   # load/length of the element
qz2=-qz2v*cosA2;  # projection onto z-axis
qx2=qz2v*sinA2; # projection onto x-axis


qz2=-qz3;
qx2=-qx3;
#
Lpunkt1=0.4*h;
Lpunkt2=0.3*h;
Lpunkt3=0.4*L5;
#
#Fz1=0.0;
#Fz2=0.0;
Fz4=-F1;
aF4=Lpunkt1;
##Fx3=0.0;
Fz5=F3*cosA5;   % Fz2 projection onto z-axis
Fx5=F3*sinA5;  % Fx2 projection onto x-axis
aF5=Lpunkt3;
#
Fz6=F2;
aF6=Lpunkt2;
#Fz7=0.0;
#
#disp(' Element load  in local coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN/m^2
LoadsqONelement=4;
esQkoormu=zeros(LoadsqONelement,4,ElementideArv);
esQkoormus(1,1:4,1)=[q1 0.0 0.0 L1];
esQkoormus(1,1:4,2)=[qz2 qx2 0.0 L2];
esQkoormus(1,1:4,3)=[qz3 qx3 0.0 L3];
esQkoormus(1,1:4,4)=[0.0 0.0 0.0 L4];
esQkoormus(1,1:4,5)=[0.0 0.0 0.0 L5];
esQkoormus(1,1:4,6)=[0.0 0.0 0.0 L6];
esQkoormus(1,1:4,7)=[0.0 0.0 0.0 L7];
#esQkoormus
#
# Point load in local coordinate z and x direction  kN
#  Fz, Fx, aF - coordinate of the point of force application
LoadsF_on_Element=5;
esFjoud=zeros(LoadsF_on_Element,2,ElementideArv);
esFjoud(1,1:3,1)=[0.0 0.0 L1];
esFjoud(1,1:3,2)=[0.0 0.0 L2];
esFjoud(1,1:3,3)=[0.0 0.0 L3];
esFjoud(1,1:3,4)=[Fz4 0.0 aF4];
esFjoud(1,1:3,5)=[Fz5 Fx5 aF5];
esFjoud(1,1:3,6)=[Fz6 0.0 aF6];
esFjoud(1,1:3,7)=[0.0 0.0 L7];
#esFjoud

#disp(' Node forces in global  coordinates ')
# sSolmF(forces,1,nodes); forces=[Fx; Fz; My]
sSolmF = zeros(3,1,SolmedeArv);
#sSolmF = zeros(3,1,SolmedeArv);
#sSolmF(:,1,1)= 0.0
#sSolmF(:,1,2)= 0.0
#sSolmF(:,1,3)= 0.0
#sSolmF(:,1,4)= 0.0
#sSolmF(:,1,5)= 0.0

#%%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3,1,SolmedeArv);
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
#%%%%%
case{10}
disp(' Load variant 10 ')
#
disp(' Initial data  ')
l=12.0
h=5.0
Isuhe=3  % Isuhe=I1/I2
p1=10.0
p2=0.0
p3=0.0
F1=20.0
F2=0.0
F3=0.0
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr=Isuhe*EIp
#
baasi0=EIp/h
xi=0.6;
l1=l;
l2=l;
h08=0.8*h
deltah=(1.0*0.2*h)/l1;
#% sqrt(h^2+l1^2);
L1=h08;
##L2=l1;
L2=sqrt(deltah^2+1.0^2); # length of the element
L3=sqrt((0.2*h)^2+l1^2); # length of the element
L4=h;
L5=sqrt((0.2*h)^2+l2^2); # length of the element
L6=L1;
L7=L2;
#sqrt(h^2+l4^2);
Lp=l1+l2;   # graphics axis
sinA3=(0.2*h)/L3;
cosA3=l1/L3;
sinA2=sinA3;
cosA2=cosA3;
sinA5=(0.2*h)/L5;
cosA5=l2/L5;
#
h1=xi*h;
h2=(1-xi)*h;
d=h1;
H=h;
#
q1=p3;
#
  if (p1 > 0)
    p12=p1
  else
    p12=p2
   endif
#
q2=p12
q3=p12
q4=0.0;
q5=0.0;
q6=0.0;
q7=0.0;
#
qz3v=q3*l1/L3;   # load/length of the element
qz3=qz3v*cosA3;  # projection onto z-axis
qx3=-qz3v*sinA3; # projection onto x-axis
qz2v=q2*1.0/L2;   # load/length of the element
qz2=-qz2v*cosA2;  # projection onto z-axis
qx2=qz2v*sinA2; # projection onto x-axis


qz2=-qz3;
qx2=-qx3;
#
Lpunkt1=0.4*h;
Lpunkt2=0.3*h;
Lpunkt3=0.4*L5;
#
#Fz1=0.0;
#Fz2=0.0;
Fz4=-F1;
aF4=Lpunkt1;
##Fx3=0.0;
Fz5=F3*cosA5;   % Fz2 projection onto z-axis
Fx5=F3*sinA5;  % Fx2 projection onto x-axis
aF5=Lpunkt3;
#
Fz6=F2;
aF6=Lpunkt2;
#Fz7=0.0;
#
#disp(' Element load  in local coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN/m^2
LoadsqONelement=4;
esQkoormu=zeros(LoadsqONelement,4,ElementideArv);
esQkoormus(1,1:4,1)=[q1 0.0 0.0 L1];
esQkoormus(1,1:4,2)=[qz2 qx2 0.0 L2];
esQkoormus(1,1:4,3)=[qz3 qx3 0.0 L3];
esQkoormus(1,1:4,4)=[0.0 0.0 0.0 L4];
esQkoormus(1,1:4,5)=[0.0 0.0 0.0 L5];
esQkoormus(1,1:4,6)=[0.0 0.0 0.0 L6];
esQkoormus(1,1:4,7)=[0.0 0.0 0.0 L7];
#esQkoormus
#
# Point load in local coordinate z and x direction  kN
#  Fz, Fx, aF - coordinate of the point of force application
LoadsF_on_Element=5;
esFjoud=zeros(LoadsF_on_Element,2,ElementideArv);
esFjoud(1,1:3,1)=[0.0 0.0 L1];
esFjoud(1,1:3,2)=[0.0 0.0 L2];
esFjoud(1,1:3,3)=[0.0 0.0 L3];
esFjoud(1,1:3,4)=[Fz4 0.0 aF4];
esFjoud(1,1:3,5)=[Fz5 Fx5 aF5];
esFjoud(1,1:3,6)=[Fz6 0.0 aF6];
esFjoud(1,1:3,7)=[0.0 0.0 L7];
#esFjoud

#disp(' Node forces in global  coordinates ')
# sSolmF(forces,1,nodes); forces=[Fx; Fz; My]
sSolmF = zeros(3,1,SolmedeArv);
#sSolmF = zeros(3,1,SolmedeArv);
#sSolmF(:,1,1)= 0.0
#sSolmF(:,1,2)= 0.0
#sSolmF(:,1,3)= 0.0
#sSolmF(:,1,4)= 0.0
#sSolmF(:,1,5)= 0.0

#%%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3,1,SolmedeArv);
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
#%%%%%
#%%%%%
otherwise
disp(' No load variant cases ')
endswitch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 disp(' Elements loads q in local coordinates ')
#%%%qz=qzZ'
#%%%qx=qxZ'
disp(' qz, qx, qA, qL  ')
disp('  For the element load q = 0,  qL = L - the end of element  ')
esQkoormus
disp(' Elements forces F  in local coordinates  ')
#%%%Fx=FZx'
#%%%Fz=FZz'
##esFjoud
disp(' Fz, Fx, aF -- coordinate of the point of force F  application in local coordinates ')
#disp(' x = aL -- coordinate of the point of force F  application in local coordinates ')
disp(' For the element force F = 0,  aL = L - the end of element  ')
#%aL=aLXx'
esFjoud
#%%%
#disp(' Node 2 forces in global  coordinates ')
#sSolmF(:,:,2)
#%s2F'
#disp(' Node 4 forces in global  coordinates ')
#sSolmF(:,:,4)
#%s4F'
#disp(' Node 6 forces in global  coordinates ')
#sSolmF(:,:,6)
#s6F'
disp('------  ')
disp(' Node forces in global  coordinates ')
disp(' Fx Fz My; : ; Node number; ')
disp('------  ')
sSolmF(:,1,:)=sSolmF(:,1,:)
#
#==========
# Nodal coordinates
#==========
deltah=(1.0*0.2*h)/l1;
krdn=[#	 x                z
        0.0               0.0;   % node 1
        0.0              -h08;   % node 2
       -1.0       -h08+deltah;   % node 3
	l1                 -h;   % node 4
        l1                0.0;   % node 5
        l1+l2            -h08;   % node 6
        l1+l2             0.0;   % node 7
        l1+l2+1.0 -h08+deltah];   % node 8
#==========
 # disp(' No      toel  u  w  fi  N   Q   M ')
#==========
#  Restrtictions on  support displacements.
# Support  No   u   w  fi   hold on - 1, open - 0
#==========
tsolm=[1       1 1 1;   % node 1
       5       1 1 1;   % node 3
       7       1 1 0];  % node 5
#==========
# Support shift multiplied by scaling multiplier
#tSiire(1,1,1)= 0.0
#tSiire(2,1,1)= 0.01*baasi0
#tSiire(1,1,3)= 0.05
#tSiire(2,1,3)= 0.0
#tSiire(3,1,3)= 0.0
#tSiire(1,1,5)= 0.0
#tSiire(2,1,5)= 0.0
#==========
 # ------------- Element properties, topology and hinges ---------
#-----------------------
elasts=[# Element properties,
#                     n2 - end of the element,
#                          n1 - beginning of the element,
#                               N, Q, M hinges at end of the element
#                                        N, Q, M hinges at beginning of the element.
	EIp EAp GAp   2     1    0 0 1    0 0 0;    % element 1
	EIr EAr GAr   3     2    0 0 0    0 0 0;    % element 2	
	EIr EAr GAr   4     2    0 0 0    0 0 0;    % element 3
	EIp EAp GAp   5     4    0 0 0    0 0 0;    % element 4
	EIr EAr GAr   6     4    0 0 0    0 0 0;    % element 5
	EIp EAp GAp   7     6    0 0 1    0 0 0;    % element 6
	EIr EAr GAr   8     6    0 0 0    0 0 0];   % element 7
# 1-hinge 'true' (axial, shear, moment hinge)
#
#==========

#==========
SARV=size(krdn);
NSARV=SARV(1,1);
#%%%%%%%%%%%%
for i=1:NSARV
sJrN(i)=i;
end
#
sJrNT=sJrN';
#%%%%%%%%%%%
#disp('-----------------------------')
disp('===============================')
disp(' Nodal coordinates ')
disp('  Node     X        Z  ')
disp('-----------------------------')
#disp([krdn])
for i=1:NSARV
disp(sprintf('  %2i    %7.4f  %7.4f',sJrNT(i),    krdn(i,1), krdn(i,2)))
endfor
disp('-----------------------------')
#%%%
if ((SolmedeArv-NSARV) != 0)
disp(' ')
disp(' ------- ')
disp(' ')
Number_of_frame_nodes
Error_by_nodes_input=sprintf('Number of nodes must be, but is: %d, %d', SolmedeArv, NSARV)
error(' Error by nodes input.')
disp(' ')
endif
#%%%

disp('===============================================')
disp('        Topology and  hinges. ')
disp('  End of the element, beginning of the element,  ')
disp('  axial-, shear-, moment hinge; 1 - hinge ''true'' ')
disp('-----------------------------------------------')
elemkinni=elasts(:,4:11);
Topology_and_and_hinges=elemkinni
EARV=size(elemkinni);
NEARV=EARV(1,1);
NEARV;
#
#%%%%%%%%%%%%%%%%%%%%
# --------The displacements and the forces numbers at the end and beginning of element.  ----- -----
slmnr=0;
for i=1:NEARV
 for j=1:12
   slmnr=slmnr+1;
   selemjln(i,j)=slmnr;
 endfor
endfor
#
disp('====================================================')
disp('  Displacements and  forces numbers (DaFs) ')
disp('  at the end and at the beginning of elements.  ')
disp('----------------------------------------------------')
selemjln;
Displacements_and_the_forces_numbers=selemjln
#
selemjl=[selemjln elasts];
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#selemjl=[1  2  3  4  5  6  7  8  9 10 11 12 EIp EAp GAp 2 1 ;  # element 1
#	13 14 15 16 17 18 19 20 21 22 23 24 EIr EAr GAr 4 2 ;	# element 2
#	25 26 27 28 29 30 31 32 33 34 35 36 EIp EAp GAp 4 3 ;	# element 3
#	37 38 39 40 41 42 43 44 45 46 47 48 EIr EAr GAr 6 4 ;	# element 4
#	49 50 51 52 53 54 55 56 57 58 59 60 EIp EAp GAp 5 6 ];	# element 5
#
selem=[selemjl(:,1:23)];
#EARV=size(selem);
#NEARV=EARV(1,1);
#%%%%%%%%%%%%
for i=1:NEARV
JrN(i,1)=i;
endfor
#
JrNT=JrN';
#%%%%%%%%%%%
disp('========================================================================================================================================')
disp(' Topology of elements                                                                                                       Moment hinge')
disp('                                                                                                               Node number   ''true''=1 ')
disp(' No.  At the end:  u, w, fi, N, Q, M.  At the beginning: u, w, fi, N, Q, M.      EI         EA          GA    End Beginning  M    M')
disp('----------------------------------------------------------------------------------------------------------------------------------------')
#
for i=1:NEARV
disp(sprintf('  %2i     %2i   %2i   %2i   %2i   %2i   %2i     %2i   %2i   %2i    %2i    %2i    %2i     %9.3e   %9.3e   %9.3e  %2i  %3i    %4i   %2i',JrN(i,1),    selemjl(i,1), selemjl(i,2), selemjl(i,3),selemjl(i,4), ...
selemjl(i,5), selemjl(i,6), selemjl(i,7),selemjl(i,8), ...
selemjl(i,9), selemjl(i,10), selemjl(i,11),selemjl(i,12), ...
selemjl(i,13), selemjl(i,14), selemjl(i,15), selemjl(i,16),    selemjl(i,17),   selemjl(i,20),  selemjl(i,23)))
endfor
disp('----------------------------------------------------------------------------------------------------------------------------------------')
#%%%
EARV=size(selem);
NEARV=EARV(1,1);
#%%%
if ((ElementideArv-NEARV)  !=0 )
disp(' ')
disp(' ------- ')
disp(' ')
Number_of_elements
Error_by_elements_input=sprintf(' Number of elements must be, and is  %d, %d',ElementideArv,NEARV)
disp(' ')
error(' Error inserting the topology of elements.')
endif
#%%%
#==========
TSARV=size(tsolm);
NTSARV=TSARV(1,1);
#%%%%%%%%%%%%
for i=1:NTSARV
sJrN(i)=i;
end
#
sTJrNT=sJrN';
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Restrtictions on  support displacements.
# disp(' No      toel  u  w  fi  N   Q   M ')
%disp('--------------------------------')
disp('==================================================')
disp(' Restrtictions on  support displacements. ')
disp(' No    Node No   u  w  fi;  hold on-1, open-0')
disp('--------------------------------------------------')
#disp([tsolm])
for i=1:NTSARV
disp(sprintf('  %2i     %5i   %2i %2i %2i %2i %2i %2i',sTJrNT(i),    tsolm(i,1), tsolm(i,2), tsolm(i,3), tsolm(i,4)))
endfor
disp('--------------------------------------------------')
%%%
TOEARV=size(tsolm);
TOARV=TOEARV(1,1);
##%%%%%%%%%%%
ToeSidemeteArv=0;
for i=1:NTSARV
  for j=2:4
     if (tsolm(i,j) == 1)
       ToeSidemeteArv=ToeSidemeteArv+1;
     endif
  endfor
endfor
#
#%%%
#ToeSidemeteArv
if ((Ntoerkts-ToeSidemeteArv)  !=0 )
disp(' ')
disp(' ------- ')
disp(' ')
Number_of_support_reactions
Error_by_inserting_support_reactions=sprintf('Number of support reactions must be, and is  %d, %d',Ntoerkts,ToeSidemeteArv)
disp(' ')
error(' Error inserting the restrtictions on the support displacements.')
endif
#
#%%%
#
lvarras=VardaPikkus(NSARV,NEARV,krdn,selem);
#lvarras
#%%%%%%%%%%%
disp('   ')
NEARV;
NNK=12*NEARV+Ntoerkts;

lvarras;
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
for i=1:NEARV
elemendiN(i,1)=i;
siireVardaA(i,1:3)=selemjl(i,7:9);
siireVardaL(i,1:3)=selemjl(i,1:3);
joudVardaA(i,1:3)=selemjl(i,10:12);
joudVardaL(i,1:3)=selemjl(i,4:6);
endfor
#
siireVardaA;
siireVardaL;
joudVardaA;
joudVardaL;
#%%%%%%%%%%%%


%%%%%ABB=ABB


#--------- Define a matrix spA to be sparse -----------
disp('-------- Sparse matrix instantiation --------  ')
disp('  spA=sparse(NNK,NNK)  ')
spA=sparse(NNK,NNK) % nullistan võrrandisüsteemi kordajate maatriksi
#--------- Define a vector B ------------
disp('-------- Right-hand side of the equations (RHS). --------  ')
disp('  B=zeros(NNK,1);  ')
B=zeros(NNK,1);   % nullistan võrrandisüsteemi vabaliikmed
disp('   ')
%%%%%%%%%%%%%%%%%%%%%%%%%

disp('  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' function AlgParm=LaheFrameDFIm(baasi0,Ntoerkts,esQkoormus,esFjoud,sSolmF,tsolm,tSiire,krdn,selem) ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp('  ')

AlgParm=LaheFrameDFIm(baasi0,Ntoerkts,esQkoormus,esFjoud,sSolmF,tsolm,tSiire,krdn,selem);
disp('  ')
AlgPar=AlgParm;
disp('  ')


disp('  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Element displacements and forces  determined by transfer matrix ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
load_variant
disp('  ')
#-------------The element displacements and forces----------
#
  suurused=[' displacement u - ';
            ' displacement w - ';
            ' rotation    fi - ';
            ' normal force N - ';
            ' shear force  Q - ';
            ' moment force M - '];
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#i=1
for i=1:NEARV
   krda=i;
   vF=zeros(6,12);
      EI=selem(i,13);   % from topology
      EA=selem(i,14);   %   "      "
      GAr=selem(i,15);  %   "      "
   Li=lvarras(i,1);
#    SII=0.0;
#    qx=qxZ(i,1);
#    qz=qzZ(i,1);
#    aLx=aLXx(i,1);
#    Fz=FZz(i,1);
#    Fx=FZx(i,1);
baasi0=1.0;
#Fjoud=esFjoud(:,1:2,i);
Fjoud=esFjoud(:,1:3,i);
#%       Fz    Fx     a
#%Fjoud=[0.0    0.0   0.0;
#qkoormus=esQkoormus(:,1:3,i);
qkoormus=esQkoormus(:,1:4,i);
#%       qz     qx     qA      qL
#%qkoormus=[0.0   0.0   0.0  0.0   0.0];

   xsamm=Li/Nmitmeks; % element is divided into 4-th parts
   xx=0;
   AP=AlgPar(i,:)';
#---------The transfer matrix equation --------
     for ij=1:Nmitmeks+1  # 5 - displacements and forces at x=0.0
       vvF=ylfhlin(1.0,xx,EA,GAr,EI);
#       vvB=yzhqz(1.0,xx,qx,qz,EA,EI);
#       vvFz=yzfzv(1.0,xx,aLx,Fx,Fz,EA,EI);
 Sj=ESTFrKrmus(1.0,xx,Li,Fjoud,qkoormus,EA,EI);
 vB=Sj;
   Fvv(:,ij)=vvF*AP+vB;
#       Fvv(:,ij)=vvF*AP+vvB+vvFz;
       xx=xx+xsamm;
     endfor
#-------------The element displacements and forces----------
   #
   VardaNr=i;
   disp(sprintf('%13s %2i %15s %5.3f %2s', 'Displacements and forces of element no',VardaNr,'  length of element ',Li,'m,  '))
    disp('    element is divided into 4 th parts ')
   %disp('  ')
#
     for i=1:3
       disp(sprintf('%14s %9.5e   %9.5e   %9.5e   %9.5e  %9.5e',suurused(i,:), Fvv(i,1),   Fvv(i,2),   Fvv(i,3),  Fvv(i,4),  Fvv(i,5)))
     endfor
#
     for i=4:6
       disp(sprintf('%14s %9.5f   %9.5f   %9.5f   %9.5f  %9.5f',suurused(i,:), Fvv(i,1),   Fvv(i,2),   Fvv(i,3),  Fvv(i,4),  Fvv(i,5)))
     endfor
#
   #disp('  ')
#
      disp('----------------------------------------------------------------------------')
#
endfor
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##SF=SisejoudPunktism(VardaNr,X,AlgPar,lvarras,selem,esQkoormus,esFjoud,suurused)
SF=SisejoudPunktism(4,Lpunkt1-eps1,AlgPar,lvarras,selem,esQkoormus,esFjoud,suurused);
disp('-----------------------------------')
SF=SisejoudPunktism(4,Lpunkt1,AlgPar,lvarras,selem,esQkoormus,esFjoud,suurused);
disp('-----------------------------------')
SF=SisejoudPunktism(5,Lpunkt3-eps1,AlgPar,lvarras,selem,esQkoormus,esFjoud,suurused);
disp('-----------------------------------')
SF=SisejoudPunktism(5,Lpunkt3,AlgPar,lvarras,selem,esQkoormus,esFjoud,suurused);
disp('-----------------------------------')
SF=SisejoudPunktism(6,Lpunkt2-eps1,AlgPar,lvarras,selem,esQkoormus,esFjoud,suurused);
disp('-----------------------------------')
SF=SisejoudPunktism(6,Lpunkt2,AlgPar,lvarras,selem,esQkoormus,esFjoud,suurused);
disp('-----------------------------------')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%
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
#%%%%%%%%%%%%
disp('  ')
disp(' Testing a static equilibrium:  ')
disp('  ')
disp(' Support_reactions=[X(85) X(86) X(87) X(88) X(89) X(90) X(91) X(92)] ')
Support_reactions=[X(85) X(86) X(87) X(88) X(89) X(90) X(91) X(92)]
disp('The_lengths=[l1 l2 h08 h1  h] ')
The_lengths=[l1 l2 h08 h1  h]
## q2=p12; 61+24=85  -- 92
disp(' q3=p1 or q3=p2  ')
q3=q3
p3=p3
disp('  ')
disp(' sumZ=q3*(l1+1.0)+F3+X(86)+X(89)+X(92) ')
sumZ=q3*(l1+1.0)+F3+X(86)+X(89)+X(92)
disp('  ')
disp(' sumX=p3*h08+F1-F2+X(85)+X(88)+X(91)')
sumX=p3*h08+F1-F2+X(85)+X(88)+X(91)
disp('  ')
disp(' Sum of the moments acting  about point 1: ')
disp(' sumM1=-p3*h08*h08/2-q3*(l1+1.0)*((l1+1.0)/2-1.0)-F1*0.6*h... ')
disp('       -F3*(l1+0.4*l2)+F2*0.5*h+X(87)+X(90)-l1*X(89)-(l1+l2)*X(92) ')
sumM1=-p3*h08*h08/2-q3*(l1+1.0)*((l1+1.0)/2-1.0)-F1*h1-F3*(l1+0.4*l2)+...
F2*0.5*h+X(87)+X(90)-l1*X(89)-(l1+l2)*X(92)
disp('  ')
disp('  Sum of the moments acting  about point 5:  ')
disp(' sumM7=-p3*h08*h08/2+q3*(l1+1.0)*((l1+1.0)/2+l2)-F1*0.6*h... ')
disp('       +F3*(0.6*l2)+F2*0.5*h+X(87)+X(90)+l1*X(89)+(l1+l2)*X(86) ')
sumM7=-p3*h08*h08/2+q3*(l1+1.0)*((l1+1.0)/2+l2)-F1*0.6*h+F3*(0.6*l2)+F2*0.5*h+X(87)+X(90)+l1*X(89)+(l1+l2)*X(86)
disp('  ')
disp(' Calculations verified the static equilibrium of the frame   ')
disp('  ')
#%%%%%%%%%%%%
figure(1)
   ax = gca ();
   #axy=[ 0   0   1   0.55];
   #joonepaksus=0.50000
   #joonepaksus=0.75000
   #joonepaksus=3
   #set (ax, "outerposition", axy);
   #set (ax, "outerposition", [ 0   0   1   0.55])
   #ax = gca ();
   #set (ax, "linewidth", joonepaksus);
   #
   grid off
   #
   Lpikk=Lp;
   d=4.0;
   #
     for j=1:NEARV
         IR=j;
         #
         #   axis([-d Lpikk+d -H-2.0 1.0], "ij")
#         axis([-d Lpikk+d  -Lpikk-0.2*d d], "ij")
         axis([-d Lpikk+d  -H-8.0 0.5*d], "ij")
         #
         plot(VGRx(IR,:),VGRz(IR,:),"3")
         hold on
     endfor
   #
   jaotT=2;
   #
   title('spESTframe1DefWFI',"fontsize",12)
   ax=gca();
      xlabel('x',"fontsize",12)
      ylabel('z',"fontsize",12)
   xticks(1,:)=-d:jaotT:Lpikk+d;
   set (ax, "xtick", xticks)
   #
   #%Sõlmede numbrid
   TTT=' ';
   TST=' ';
   TT1=' ';
   TT1=' ';
   #%%%%%%%%%
    for i=1:NSARV;
      #
      IR=i;
      TTT=eval(sprintf('%i',IR));
      TST=int2str(TTT);

      #text(krdn(IR,1), krdn(IR,2)+0.17, TST)
      text(krdn(IR,1), krdn(IR,2)+0.25, TST)
      #
     endfor
#
#%% Varraste numbrid
for j=1:NEARV;
   #
   IR=j;
   TTT=eval(sprintf('%i',IR));
   TST=int2str(TTT);
   text((VGRx(IR,1)+VGRx(IR,2))/2+0.1, (VGRz(IR,1)+VGRz(IR,2))/2, TST)
   #
endfor
#text(L/2-d/2, -H-d/2, 'Raamskeem')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xvalg=-1.0;
#yvalg=-15.5;
yvalg=-H-6.0;
sammuga=0.5;
#sammuga=0.8;
text(-1.0, yvalg-1.5*sammuga, 'Numeration of displacements and forces ')
#TSA=[TSA1 Tyhk TSA2 Tyhk TSA3] ;
text(-2.0, yvalg, 'u w fi N Q M at the beginning')
#%xvalg=xvalg+sammuga;
yvalg=yvalg+sammuga;

for i=1:NEARV;
   #siireVardaA  joudVardaA(i,1:3)
   TTA1=eval(sprintf('%3i',siireVardaA(i,1)));
   TTA2=eval(sprintf('%3i',siireVardaA(i,2)));
   TTA3=eval(sprintf('%3i',siireVardaA(i,3)));
   TTA4=eval(sprintf('%3i',joudVardaA(i,1)));
   TTA5=eval(sprintf('%3i',joudVardaA(i,2)));
   TTA6=eval(sprintf('%3i',joudVardaA(i,3)));
      TSA1=int2str(TTA1);
      TSA2=int2str(TTA2);
      TSA3=int2str(TTA3);
      TSA4=int2str(TTA4);
      TSA5=int2str(TTA5);
      TSA6=int2str(TTA6);
      Tyhk='  ';
   TSA=[TSA1 Tyhk TSA2 Tyhk TSA3 Tyhk TSA4 Tyhk TSA5 Tyhk TSA6] ;
   text(xvalg, yvalg, TSA)
   #%xvalg=xvalg+sammuga;
   yvalg=yvalg+sammuga;
endfor
yvalg=yvalg+sammuga;

text(-1.0, yvalg, 'Support reactions: 85 86 87 88 89 90 91 92 ')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xvalg=7.0;
#xvalg=10.0;
#yvalg=-15.5;
yvalg=-H-6.0;
sammuga=0.5;
#sammuga=0.8;
#TSA=[TSA1 Tyhk TSA2 Tyhk TSA3] ;
text(6.0, yvalg, 'u w fi N Q M at the end')
#text(9.5, yvalg, 'u w fi N Q M at the end')
#%xvalg=xvalg+sammuga;
yvalg=yvalg+sammuga;
for i=1:NEARV;
   #siireVardaL  joudVardaL(i,1:3)
   TTA1=eval(sprintf('%3i',siireVardaL(i,1)));
   TTA2=eval(sprintf('%3i',siireVardaL(i,2)));
   TTA3=eval(sprintf('%3i',siireVardaL(i,3)));
   TTA4=eval(sprintf('%3i',joudVardaL(i,1)));
   TTA5=eval(sprintf('%3i',joudVardaL(i,2)));
   TTA6=eval(sprintf('%3i',joudVardaL(i,3)));
     TSA1=int2str(TTA1);
     TSA2=int2str(TTA2);
     TSA3=int2str(TTA3);
     TSA4=int2str(TTA4);
     TSA5=int2str(TTA5);
     TSA6=int2str(TTA6);
     Tyhk='  ';
   TSA=[TSA1 Tyhk TSA2 Tyhk TSA3 Tyhk TSA4 Tyhk TSA5 Tyhk TSA6] ;
   text(xvalg, yvalg, TSA)
   #%xvalg=xvalg+sammuga;
   yvalg=yvalg+sammuga;
endfor
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
###### hold off
#
#Vaata: http://www.gnu.org/software/octave/doc/interpreter/Printing-Plots.html
#
#print('-deps','spESTframe1DefWFI.eps');
#print('-dpng','spESTframe1DefWFI.png');
print('-dfig','-landscape','-mono','-solid','spESTframe1DefWFI.fig');
refresh
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
spy(spA)
title('spy(spA,14) - the sparse matrix spA(92,92) non zero elements [3%] ')
refresh
#
#print('-deps','spESTframe1DefWFI_h6remaatriks.eps');
#print('-dpng','spESTframe1DefWFI_h6remaatriks.png');
print('-dfig','-landscape','spESTframe1DefWFI_sparse_matrix.fig');
%%print('-dfig','-landscape','-mono','-solid','spESTframe1DefWFI_sparse_matrix.fig');
refresh
#
# octave
disp('-----------------------------------')
disp('  ')
disp(' The END ')
#%%%%%%%%%%%%%%%