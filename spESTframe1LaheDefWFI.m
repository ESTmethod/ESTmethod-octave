##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-07-16
## Last modified: Andres Lahe, 2013-07-21
##                Mattias Põldaru, 2013-07-07
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

close all
clear
disp('==================================================================')
disp(' spESTframe1LaheDefWFI.m')
disp(' Analysis of a plane frame using the EST method.')
disp(' The EST method finds the displacements and forces.')
disp(' The problem: http://digi.lib.ttu.ee/opik_eme/ylesanded.pdf#page=39')
disp('---')
disp(' The equations matrix is non-symmetric sparse matrix.')
disp('---')

# The original problem:
# http://digi.lib.ttu.ee/opik_eme/ylesanded.pdf#page=39
# http://digi.lib.ttu.ee/opik_eme/slaidid/raamEST1en.pdf
# http://digi.lib.ttu.ee/opik_eme/Ehitusmehaanika.pdf#page=419

global spA
global X
global siireVardaA

disp(' daf -- the displacements and the forces at the ends of element.')
disp(' The element of plane frame have 12 daf''s (6+6).')
disp(' The number of unknowns are 12*number_of_elements + number_of_support_reactions.')
disp(' The matrixes are assembled into compressed column sparse matrices.')
disp('')

# Number nodes, elements, support reactions.
number_of_nodes = 8;
number_of_elements = 7;
number_of_support_reactions = 8;
number_of_unknowns = 12 * number_of_elements + number_of_support_reactions;

# Element properties
EIp = 20000;
EAp = 4.6 * 10^15;
EAr = 6.8 * 10^15;
GAp = 0.4 * EAp;
GAr = 0.4 * EAr;
divisions = 4;

# ---- Load variants -----
load_variant = 1;

switch (load_variant)
    case{1}
        disp(' Load variant 1 ')
        l = 6.0
        h = 3.0
        Isuhe = 2 % Isuhe=I1/I2
        p1 = 12.0
        p2 = 0.0
        p3 = 0.0
        F1 = 50.0
        F2 = 0.0
        F3 = 0.0
    case{2}
        disp(' Load variant 2 ')
        l = 8.0
        h = 4.0
        Isuhe = 2 % Isuhe=I1/I2
        p1 = 0.0
        p2 = 14.0
        p3 = 0.0
        F1 = 0.0
        F2 = 50.0
        F3 = 0.0
    case{3}
        disp(' Load variant 3 ')
        l = 10.0
        h = 4.0
        Isuhe = 2 % Isuhe=I1/I2
        p1 = 0.0
        p2 = 0.0
        p3 = 8.0
        F1 = 0.0
        F2 = 0.0
        F3 = 50.0
    case{4}
        disp(' Load variant 4 ')
        l = 6.0
        h = 4.0
        Isuhe = 2 % Isuhe=I1/I2
        p1 = 14.0
        p2 = 0.0
        p3 = 0.0
        F1 = 40.0
        F2 = 0.0
        F3 = 0.0
    case{5}
        disp(' Load variant 5 ')
        l = 8.0
        h = 5.0
        Isuhe = 2 % Isuhe=I1/I2
        p1 = 0.0
        p2 = 8.0
        p3 = 0.0
        F1 = 0.0
        F2 = 40.0
        F3 = 0.0
    case{6}
        disp(' Load variant 6 ')
        l = 10.0
        h = 5.0
        Isuhe = 3 % Isuhe=I1/I2
        p1 = 0.0
        p2 = 0.0
        p3 = 12.0
        F1 = 0.0
        F2 = 0.0
        F3 = 40.0
    case{7}
        disp(' Load variant 7 ')
        l = 6.0
        h = 3.0
        Isuhe = 3 % Isuhe=I1/I2
        p1 = 8.0
        p2 = 0.0
        p3 = 0.0
        F1 = 30.0
        F2 = 0.0
        F3 = 0.0
    case{8}
        disp(' Load variant 8 ')
        l = 8.0
        h = 4.0
        Isuhe = 3 % Isuhe=I1/I2
        p1 = 0.0
        p2 = 16.0
        p3 = 0.0
        F1 = 0.0
        F2 = 30.0
        F3 = 0.0
    case{9}
        disp(' Load variant 9 ')
        l = 10.0
        h = 4.0
        Isuhe = 3 % Isuhe=I1/I2
        p1 = 0.0
        p2 = 0.0
        p3 = 14.0
        F1 = 0.0
        F2 = 0.0
        F3 = 30.0
    case{10}
        disp(' Load variant 10 ')
        l = 12.0
        h = 5.0
        Isuhe = 3 % Isuhe=I1/I2
        p1 = 10.0
        p2 = 0.0
        p3 = 0.0
        F1 = 20.0
        F2 = 0.0
        F3 = 0.0
endswitch

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%
EIr = Isuhe * EIp;
baasi0 = EIp / h;
xi = 0.6;
h08 = 0.8 * h;
deltah = (1.0 * 0.2 * h) / l;
L1 = h08;
L2 = sqrt(deltah^2 + 1.0 ^ 2); # length of the element
L3 = sqrt((0.2 * h)^2 + l ^ 2); # length of the element
L4 = h;
L5 = sqrt((0.2 * h)^2 + l ^ 2); # length of the element
L6 = L1;
L7 = L2;

Lp = 2 * l; # graphics axis
sinA3 = (0.2 * h) / L3;
cosA3 = l / L3;
sinA2 = sinA3;
cosA2 = cosA3;
sinA5 = (0.2 * h) / L5;
cosA5 = l / L5;

h1 = xi * h;
h2 = (1 - xi) * h;
d = h1;
H = h;

q1 = p3;
q = max(p1, p2)

qz3v = q * l / L3; # load / length of the element
qz3 = qz3v * cosA3; # projection onto z - axis
qx3 = - qz3v * sinA3; # projection onto x - axis
qz2v = q * 1.0 / L2; # load / length of the element
qz2 = - qz2v * cosA2; # projection onto z - axis
qx2 = qz2v * sinA2; # projection onto x - axis

qz2 = - qz3;
qx2 = - qx3;

Lpunkt1 = 0.4 * h;
Lpunkt2 = 0.3 * h;
Lpunkt3 = 0.4 * L5;

Fz4 = - F1;
aF4 = Lpunkt1;
Fz5 = F3 * cosA5; % Fz2 projection onto z-axis
Fx5 = F3 * sinA5; % Fx2 projection onto x-axis
aF5 = Lpunkt3;
Fz6 = F2;
aF6 = Lpunkt2;

#disp(' Element load  in local node_coordinates ')
#disp('    qz     qx     qA      qL ')
# Uniformly distributed load in local coordinate z and x direction kN / m^2
esQkoormus = zeros(1, 4, number_of_elements);
esQkoormus(1, 1:4, 1) = [q1 0.0 0.0 L1];
esQkoormus(1, 1:4, 2) = [qz2 qx2 0.0 L2];
esQkoormus(1, 1:4, 3) = [qz3 qx3 0.0 L3];

# Point load in local coordinate z and x direction kN
# Fz, Fx, aF - coordinate of the point of force application
esFjoud = zeros(1, 3, number_of_elements);
esFjoud(1, 1:3, 1) = [0.0 0.0 L1];
esFjoud(1, 1:3, 2) = [0.0 0.0 L2];
esFjoud(1, 1:3, 3) = [0.0 0.0 L3];
esFjoud(1, 1:3, 4) = [Fz4 0.0 aF4];
esFjoud(1, 1:3, 5) = [Fz5 Fx5 aF5];
esFjoud(1, 1:3, 6) = [Fz6 0.0 aF6];
esFjoud(1, 1:3, 7) = [0.0 0.0 L7];

#disp(' Node forces in global node_coordinates ')
# sSolmF(forces, 1, nodes); forces = [Fx; Fz; My]
sSolmF = zeros(3, 1, number_of_nodes);

# %%%%%
# Support shift - tSiire#
# Support shift multiplied by scaling multiplier
tSiire = zeros(3, 1, number_of_nodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' Elements loads q in local node_coordinates ')
disp(' qz, qx, qA, qL  ')
disp('  For the element load q = 0,  qL = L - the end of element  ')
esQkoormus
disp(' Elements forces F in local node_coordinates  ')
disp(' Fz, Fx, aF -- coordinate of the point of force F  application in local node_coordinates ')
#disp(' x = aL -- coordinate of the point of force F  application in local node_coordinates ')
disp(' For the element force F = 0,  aL = L - the end of element  ')
esFjoud

disp('------  ')
disp(' Node forces in global  node_coordinates ')
disp(' Fx Fz My; : ; Node number; ')
disp('------  ')
sSolmF(:, 1, :)

# ==========
# Node coordinates
# ==========
deltah = (1.0 * 0.2 * h) / l;

node_coordinates = [# x  z
          0.0       0.0;          % node 1
          0.0       -h08;         % node 2
          -1.0      -h08+deltah;  % node 3
          l         -h;           % node 4
          l         0.0;          % node 5
          l*2       -h08;         % node 6
          l*2       0.0;          % node 7
          l*2+1.0   -h08+deltah]; % node 8
node_count = size(node_coordinates)(1);
# ==========
# disp(' No      toel  u  w  fi  N   Q   M ')
# ==========
# Restrictions on support displacements.
# Support node u w fi hold on - 1, open - 0
# ==========
support_nodes = [ 1 1 1 1;
                  5 1 1 1;
                  7 1 1 0];

# ------------- Element properties, topology and hinges ---------
element_properties = [# Element properties,
# n2 - end of the element,
# n1 - beginning of the element,
# N, Q, M hinges at end of the element
# N, Q, M hinges at beginning of the element.
EIp EAp GAp 2 1 0 0 1 0 0 0; % element 1
EIr EAr GAr 3 2 0 0 0 0 0 0; % element 2
EIr EAr GAr 4 2 0 0 0 0 0 0; % element 3
EIp EAp GAp 5 4 0 0 0 0 0 0; % element 4
EIr EAr GAr 6 4 0 0 0 0 0 0; % element 5
EIp EAp GAp 7 6 0 0 1 0 0 0; % element 6
EIr EAr GAr 8 6 0 0 0 0 0 0];% element 7
# 1 - hinge 'true' (axial, shear, moment hinge)
# ==========

disp('=============================')
disp(' Nodal node_coordinates ')
disp('  Node      X        Z  ')
disp('-----------------------------')
for i = 1:node_count
    disp(sprintf('  %2i     %7.4f  %7.4f', i, node_coordinates(i, 1:2)))
endfor
disp('-----------------------------')
# %%%
if number_of_nodes != node_count
    display(sprintf('Number of nodes must be %d, but is %d', number_of_nodes, node_count))
    error(' Faulty node data.')
endif

disp('===============================================')
disp('        Topology and  hinges. ')
disp('  End of the element, beginning of the element,  ')
disp('  axial-, shear-, moment hinge; 1 - hinge ''true'' ')
disp('-----------------------------------------------')
topology_and_hinges = element_properties(:, 4:11)

# -------- Displacements and force numbers at the end and the beginning of element. ----------
element_count = size(element_properties)(1);
for i = 1 : element_count
    displacements_forces_numbers(i, :) = [i*12 - 11 : i*12];
endfor
#
disp('====================================================')
disp('  Displacements and  forces numbers (DaFs) ')
disp('  at the end and at the beginning of elements.  ')
disp('----------------------------------------------------')
displacements_forces_numbers;
#
disp_force_numbers_properties = [displacements_forces_numbers element_properties];
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#disp_force_numbers_properties = [
#  1  2  3  4  5  6  7  8  9 10 11 12 EIp EAp GAp 2 1; # element 1
# 13 14 15 16 17 18 19 20 21 22 23 24 EIr EAr GAr 4 2; # element 2
# 25 26 27 28 29 30 31 32 33 34 35 36 EIp EAp GAp 4 3; # element 3
#
disp('========================================================================================================================================')
disp(' Topology of elements                                                                                                       Moment hinge')
disp('                                                                                                               Node number   ''true''=1 ')
disp(' No  At the end:  u, w, fi, N, Q, M.  At the beginning: u, w, fi, N, Q, M.      EI         EA          GA    End Beginning  M    M')
disp('----------------------------------------------------------------------------------------------------------------------------------------')
#
for i = 1:element_count
    disp(sprintf('%2i  %2i  %2i  %2i  %2i %2i  %2i  %2i  %2i  %2i  %2i  %2i  %2i  %9.2e %9.2e %9.2e %2i %2i %2i %2i', i, disp_force_numbers_properties(i, 1:17), disp_force_numbers_properties(i, 20), disp_force_numbers_properties(i, 23)))
endfor
disp('----------------------------------------------------------------------------------------------------------------------------------------')


# Checking data
element_count = size(element_properties)(1)
if number_of_elements != element_count
    disp(sprintf('Number of elements must  %d, but is %d', number_of_elements, element_count))
    error('Faulty elements data')
endif

disp('==================================================')
disp(' Restrictions on support displacements. ')
disp(' No    Node No   u  w  fi;  hold on-1, open-0')
disp('--------------------------------------------------')
n = size(support_nodes)(1);
for i = 1 : n
    disp(sprintf('  %2i     %5i   %2i %2i %2i %2i %2i', i, support_nodes(i, 1:4)))
endfor
disp('--------------------------------------------------')

# Checking data
support_reactions_count = sum(sum(support_nodes(:, 2:4)))
if number_of_support_reactions != support_reactions_count
    disp(sprintf('Number of support reactions must be %d, but is %d', number_of_support_reactions, support_reactions_count))
    error('Faulty support reactions data')
endif

lvarras = VardaPikkus(node_count, element_count, node_coordinates, disp_force_numbers_properties);

for i = 1:element_count
    elemendiN(i, 1) = i;
    siireVardaA(i, 1:3) = disp_force_numbers_properties(i, 7:9);
    siireVardaL(i, 1:3) = disp_force_numbers_properties(i, 1:3);
    joudVardaA(i, 1:3) = disp_force_numbers_properties(i, 10:12);
    joudVardaL(i, 1:3) = disp_force_numbers_properties(i, 4:6);
endfor

# Size of equation system matrix
NNK = 12 * element_count + number_of_support_reactions;

# Do the actual calculation.
AlgPar = LaheFrameDFIm(baasi0, number_of_support_reactions, esQkoormus, esFjoud, sSolmF, support_nodes, tSiire, node_coordinates, disp_force_numbers_properties);

disp('  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' Element displacements and forces  determined by transfer matrix ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
load_variant
disp('  ')
# ------------- The element displacements and forces ----------
row_names = [
    'displacement u';
    'displacement w';
    'rotation    fi';
    'normal force N';
    'shear force  Q';
    'moment force M'];
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#i = 1
for i = 1:element_count
    EI = disp_force_numbers_properties(i, 13); % from topology
    EA = disp_force_numbers_properties(i, 14);
    GAr = disp_force_numbers_properties(i, 15);
    Li = lvarras(i, 1);
    baasi0 = 1.0;
    Fjoud = esFjoud(:, 1:3, i);
    # %       Fz    Fx     a
    # %Fjoud=[0.0    0.0   0.0;
    #qkoormus = esQkoormus(:, 1:3, i);
    qkoormus = esQkoormus(:, 1:4, i);
    # %       qz     qx     qA      qL
    # %qkoormus=[0.0   0.0   0.0  0.0   0.0];
 
    xsamm = Li / divisions; % element is divided into 4-th parts
    xx = 0;
    AP = AlgPar(i, :)';
    # --------- The transfer matrix equation --------
    for ij = 1:divisions + 1 # 5 - displacements and forces at x = 0.0
        vvF = ylfhlin(1.0, xx, EA, GAr, EI);
        # vvB = yzhqz(1.0, xx, qx, qz, EA, EI);
        # vvFz = yzfzv(1.0, xx, aLx, Fx, Fz, EA, EI);
        Sj = ESTFrKrmus(1.0, xx, Li, Fjoud, qkoormus, EA, EI);
        vB = Sj;
        Fvv(:, ij) = vvF * AP + vB;
        # Fvv(:, ij) = vvF * AP + vvB + vvFz;
        xx += xsamm;
    endfor
    # ------------- Element displacements and forces ----------
    disp(sprintf('Element %i (l=%5.3f m)', i, Li))

    for i = 1:3
        disp(sprintf(' %s   %9.2e   %9.2e   %9.2e   %9.2e  %9.2e', row_names(i), Fvv(i, 1:5)))
    endfor

    for i = 4:6
        disp(sprintf(' %s   %9.5f   %9.5f   %9.5f   %9.5f  %9.5f', row_names(i), Fvv(i, 1:5)))
    endfor
    disp('--------------------------------------------------------------------------')
endfor

for i = 1:element_count
    LkoordN = disp_force_numbers_properties(i, 16);
    AkoordN = disp_force_numbers_properties(i, 17);
    DeltaX(i) = node_coordinates(LkoordN, 1) - node_coordinates(AkoordN, 1);
    DeltaZ(i) = node_coordinates(LkoordN, 2) - node_coordinates(AkoordN, 2);
    VGRx(i, 1) = node_coordinates(AkoordN, 1);
    VGRx(i, 2) = node_coordinates(LkoordN, 1);
    VGRz(i, 1) = node_coordinates(AkoordN, 2);
    VGRz(i, 2) = node_coordinates(LkoordN, 2);
endfor

disp('  ')
disp(' Testing for static equilibrium:  ')
disp('  ')
disp(' Support_reactions=[X(85) X(86) X(87) X(88) X(89) X(90) X(91) X(92)] ')
Support_reactions = [X(85) X(86) X(87) X(88) X(89) X(90) X(91) X(92)]
disp('The_lengths=[l1 l2 h08 h1  h] ')
The_lengths = [l l h08 h1 h]
disp(' q3=p1 or q3=p2  ')
q3 = q
p3 = p3
disp('  ')
disp(' sumZ=q3*(l1+1.0)+F3+X(86)+X(89)+X(92) ')
sumZ = q3 * (l + 1.0) + F3 + X(86) + X(89) + X(92)
disp('  ')
disp(' sumX=p3*h08+F1-F2+X(85)+X(88)+X(91)')
sumX = p3 * h08 + F1 - F2 + X(85) + X(88) + X(91)
disp('  ')
disp(' Sum of the moments acting  about point 1: ')
disp(' sumM1=-p3*h08*h08/2-q3*(l1+1.0)*((l1+1.0)/2-1.0)-F1*0.6*h... ')
  disp('       -F3*(l1+0.4*l2)+F2*0.5*h+X(87)+X(90)-l1*X(89)-(l1+l2)*X(92) ')
sumM1 = -p3 * h08 * h08 / 2 - q3 * (l + 1.0) * ((l + 1.0) / 2 - 1.0) ...
  - F1 * h1 - F3 * (l + 0.4 * l) + ...
  F2 * 0.5 * h + X(87) + X(90) - l * X(89) - (l*2) * X(92)
disp('  ')
disp('  Sum of the moments acting  about point 5:  ')
disp(' sumM5=-p3*h08*h08/2+q3*(l1+1.0)*((l1+1.0)/2+l2)-F1*0.6*h... ')
  disp('       +F3*(0.6*l2)+F2*0.5*h+X(87)+X(90)+l1*X(89)+(l1+l2)*X(86) ')
sumM7 = - p3 * h08 * h08 / 2 + q3 * (l + 1.0) * ((l + 1.0) / 2 + l) ...
  - F1 * 0.6 * h + F3 * (0.6 * l) + F2 * 0.5 * h + X(87) + X(90) + l * X(89) + (l + l) * X(86)
disp('  ')
disp(' Calculations verified the static equilibrium of the frame   ')
disp('  ')


exit
# %%%%%%%%%%%%
figure(1)
ax = gca ();
#axy = [0 0 1 0.55];
#joonepaksus = 0.50000
#joonepaksus = 0.75000
#joonepaksus = 3
#set (ax, "outerposition", axy);
#set (ax, "outerposition", [0 0 1 0.55])
#ax = gca ();
#set (ax, "linewidth", joonepaksus);
#
grid off
#
Lpikk = Lp;
d = 4.0;
#
for j = 1:element_count
    IR = j;
    axis([-d Lpikk+d -H-8.0 0.5*d], "ij")
    plot(VGRx(IR, :), VGRz(IR, :), "3")
    hold on
endfor
#
jaotT = 2;
#
title('spESTframe1DefWFI', "fontsize", 12)
ax = gca();
xlabel('x', "fontsize", 12)
ylabel('z', "fontsize", 12)
xticks(1, :) = - d:jaotT:Lpikk + d;
set (ax, "xtick", xticks)
#
# %Sõlmede numbrid
TTT = ' ';
TST = ' ';
TT1 = ' ';
TT1 = ' ';
# %%%%%%%%%
for i = 1:node_count;
    #
    IR = i;
    TTT = eval(sprintf('%i', IR));
    TST = int2str(TTT);
 
    #text(node_coordinates(IR, 1), node_coordinates(IR, 2) + 0.17, TST)
    text(node_coordinates(IR, 1), node_coordinates(IR, 2) + 0.25, TST)
    #
endfor
#
# %% Varraste numbrid
for j = 1:element_count;
    #
    IR = j;
    TTT = eval(sprintf('%i', IR));
    TST = int2str(TTT);
    text((VGRx(IR, 1) + VGRx(IR, 2)) / 2 + 0.1, (VGRz(IR, 1) + VGRz(IR, 2)) / 2, TST)
    #
endfor
#text(L / 2 - d / 2, - H - d / 2, 'Raamskeem')
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xvalg = - 1.0;
#yvalg = - 15.5;
yvalg = - H - 6.0;
sammuga = 0.5;
#sammuga = 0.8;
text(- 1.0, yvalg - 1.5 * sammuga, 'Numeration of displacements and forces ')
#TSA = [TSA1 Tyhk TSA2 Tyhk TSA3];
text(- 2.0, yvalg, 'u w fi N Q M at the beginning')
# %xvalg=xvalg+sammuga;
yvalg = yvalg + sammuga;

for i = 1:element_count;
    #siireVardaA joudVardaA(i, 1:3)
    TTA1 = eval(sprintf('%3i', siireVardaA(i, 1)));
    TTA2 = eval(sprintf('%3i', siireVardaA(i, 2)));
    TTA3 = eval(sprintf('%3i', siireVardaA(i, 3)));
    TTA4 = eval(sprintf('%3i', joudVardaA(i, 1)));
    TTA5 = eval(sprintf('%3i', joudVardaA(i, 2)));
    TTA6 = eval(sprintf('%3i', joudVardaA(i, 3)));
    TSA1 = int2str(TTA1);
    TSA2 = int2str(TTA2);
    TSA3 = int2str(TTA3);
    TSA4 = int2str(TTA4);
    TSA5 = int2str(TTA5);
    TSA6 = int2str(TTA6);
    Tyhk = '  ';
    TSA = [TSA1 Tyhk TSA2 Tyhk TSA3 Tyhk TSA4 Tyhk TSA5 Tyhk TSA6];
    text(xvalg, yvalg, TSA)
    # %xvalg=xvalg+sammuga;
    yvalg = yvalg + sammuga;
endfor
yvalg = yvalg + sammuga;

text(- 1.0, yvalg, 'Support reactions: 85 86 87 88 89 90 91 92 ')
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xvalg = 7.0;
#xvalg = 10.0;
#yvalg = - 15.5;
yvalg = - H - 6.0;
sammuga = 0.5;
#sammuga = 0.8;
#TSA = [TSA1 Tyhk TSA2 Tyhk TSA3];
text(6.0, yvalg, 'u w fi N Q M at the end')
#text(9.5, yvalg, 'u w fi N Q M at the end')
# %xvalg=xvalg+sammuga;
yvalg = yvalg + sammuga;
for i = 1:element_count;
    #siireVardaL joudVardaL(i, 1:3)
    TTA1 = eval(sprintf('%3i', siireVardaL(i, 1)));
    TTA2 = eval(sprintf('%3i', siireVardaL(i, 2)));
    TTA3 = eval(sprintf('%3i', siireVardaL(i, 3)));
    TTA4 = eval(sprintf('%3i', joudVardaL(i, 1)));
    TTA5 = eval(sprintf('%3i', joudVardaL(i, 2)));
    TTA6 = eval(sprintf('%3i', joudVardaL(i, 3)));
    TSA1 = int2str(TTA1);
    TSA2 = int2str(TTA2);
    TSA3 = int2str(TTA3);
    TSA4 = int2str(TTA4);
    TSA5 = int2str(TTA5);
    TSA6 = int2str(TTA6);
    Tyhk = '  ';
    TSA = [TSA1 Tyhk TSA2 Tyhk TSA3 Tyhk TSA4 Tyhk TSA5 Tyhk TSA6];
    text(xvalg, yvalg, TSA)
    # %xvalg=xvalg+sammuga;
    yvalg = yvalg + sammuga;
endfor
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
###### hold off
#
#Vaata: http://www.gnu.org/software/octave/doc/interpreter/Printing-Plots.html
#print('-deps', 'spESTframe1DefWFI.eps');
#print('-dpng', 'spESTframe1DefWFI.png');
print('-dfig', '-landscape', '-mono', '-solid', 'spESTframe1DefWFI.fig');
refresh
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
spy(spA)
title('spy(spA,14) - the sparse matrix spA(92,92) non zero elements [3%] ')
refresh
#print('-deps', 'spESTframe1DefWFI_h6remaatriks.eps');
#print('-dpng', 'spESTframe1DefWFI_h6remaatriks.png');
print('-dfig', '-landscape', 'spESTframe1DefWFI_sparse_matrix.fig');
%%print('-dfig','-landscape','-mono','-solid','spESTframe1DefWFI_sparse_matrix.fig');
refresh
#
# octave
disp('-----------------------------------')
disp('  ')
disp(' The END ')
# %%%%%%%%%%%%%%%

