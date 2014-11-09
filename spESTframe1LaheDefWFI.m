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
disp(' The number of unknowns are 12*number_of_elements + support_reactions_count.')
disp(' The matrices are assembled into compressed column sparse matrices.')
disp('')


# Basic information
# =================
# Number of nodes, elements, support reactions.
# It is used to find possible errors in data later.
number_of_nodes = 8;
number_of_elements = 7;
support_reactions_count = 8;
number_of_unknowns = 12 * number_of_elements + support_reactions_count;

# Element properties
EIp = 20000;
EAp = 4.6 * 10^15;
EAr = 6.8 * 10^15;
GAp = 0.4 * EAp;
GAr = 0.4 * EAr;
divisions = 4;


# Load variants
# =============
# For shorter notation we express variants as a table, one per row.
# [ l    h    Is  p1   p2   p3   F1   F2   F3]
load_variants = [
    6    3    2   12    0    0   50    0    0
    8    4    2    0   14    0    0   50    0
   10    4    2    0    0    8    0    0   50
    6    4    2   14    0    0   40    0    0
    8    5    2    0    8    0    0   40    0
   10    5    3    0    0   12    0    0   40
    6    3    3    8    0    0   30    0    0
    8    4    3    0   16    0    0   30    0
   10    4    3    0    0   14    0    0   30
   12    5    3   10    0    0   20    0    0];

disp('=======================================')
disp(' Current loads according to the variant')
disp('---------------------------------------')
load_variant = 1
# Extract the values to their names.
[l h Isuhe p1 p2 p3 F1 F2 F3] = num2cell(load_variants(load_variant, :)){:}
disp('')

# Precalculations to fill in node coordinates and force matrices later.
EIr = Isuhe * EIp;
scale = EIp / h;
xi = 0.6;
h08 = 0.8 * h;
deltah = (1.0 * 0.2 * h) / l;
L1 = h08;
L2 = sqrt(deltah^2 + 1.0 ^ 2);
L3 = sqrt((0.2 * h)^2 + l ^ 2);
L4 = h;
L5 = sqrt((0.2 * h)^2 + l ^ 2);
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
q = max(p1, p2);

qz3v = q * l / L3; # load / length of the element
qz3 = qz3v * cosA3; # projection onto z - axis
qx3 = - qz3v * sinA3; # projection onto x - axis
qz2v = q * 1.0 / L2; # load / length of the element
qz2 = - qz2v * cosA2; # projection onto z - axis
qx2 = qz2v * sinA2; # projection onto x - axis

qz2 = - qz3;
qx2 = - qx3;

Fz4 = - F1;
aF4 = 0.4 * h;
Fz5 = F3 * cosA5; % Fz2 projection onto z-axis
Fx5 = F3 * sinA5; % Fx2 projection onto x-axis
aF5 = 0.4 * L5;
Fz6 = F2;
aF6 = 0.3 * h;


# Node coordinates
# ================
# [x z]
# Each row stands for one node, numbering starts from 1.
deltah = 0.2 * h / l;
node_coordinates = [#  X         Z
                      0.0       0.0;          # nid = 1
                      0.0       -h08;         # nid = 2
                      -1.0      -h08+deltah;  # nid = 3 etc
                      l         -h;
                      l         0.0;
                      l*2       -h08;
                      l*2       0.0;
                      l*2+1.0   -h08+deltah];
node_count = size(node_coordinates, 1);


# Restrictions on support displacements
# =====================================
# [nid x y fi]
# nid - Node id.
# x - Support fixation in x-direction.
# z - Support fixation in z-direction.
# fi- Support fixation for revolving around y-axis.
# Supports: 1 - fixed, 0 - open/free.
#              [nid x y fi]
support_nodes = [ 1 1 1 1;
                  5 1 1 1;
                  7 1 1 0];


# Element properties, topology and hinges
# =======================================
# [EI EA GA end start Nend Qend Mend Nstart Qstart Mstart]
# EI, EA, GA - element stiffess
# end, start - node numbers of element
# N, Q, M at end and start - hinges for axial, shear and moment forces
# Hinges: 1 - fixed, 0 - open/free
element_topology_hinges = [
EIp EAp GAp 2 1 0 0 1 0 0 0; % eid = 1
EIr EAr GAr 3 2 0 0 0 0 0 0; % eid = 2
EIr EAr GAr 4 2 0 0 0 0 0 0; % eid = 3 etc
EIp EAp GAp 5 4 0 0 0 0 0 0;
EIr EAr GAr 6 4 0 0 0 0 0 0;
EIp EAp GAp 7 6 0 0 1 0 0 0;
EIr EAr GAr 8 6 0 0 0 0 0 0];

element_count = size(element_topology_hinges, 1);


# Element loads
# =============
# Uniformly distributed element loads in local coordinate z and x direction.
#
# eid - Element id.
# lid - Load id, one element may carry multiple loads.
# qz - Load in local z direction, kN/m.
# qx - Load in local x direction, kN/m.
#  a - Load start position from the start of the element, m.
#  l - Load length in m (make sure to not exceed the length of element).
#               eid lid        [ qz  qx  a  l ]
distributed_loads = zeros(number_of_elements, 1, 4);
distributed_loads(1, 1, 1:4) = [q1  0.0 0.0 L1];
distributed_loads(2, 1, 1:4) = [qz2 qx2 0.0 L2];
distributed_loads(3, 1, 1:4) = [qz3 qx3 0.0 L3];


# Element forces
# ==============
# Elements forces F in local coordinates z and x direction.
#
# eid - Element id.
# fid - Force id, one element may have multiple forces.
# Fz - Force in local z direction, kN.
# Fx - Force in local x direction, kN.
#  a - Force position from the start of the element, m. Please do not exceed
#          the length of the element.
#          eid fid        [ Fz  Fx  a ]
point_forces = zeros(number_of_elements, 1, 3);
point_forces(1, 1, 1:3) = [0.0 0.0  L1];
point_forces(2, 1, 1:3) = [0.0 0.0  L2];
point_forces(3, 1, 1:3) = [0.0 0.0  L3];
point_forces(4, 1, 1:3) = [Fz4 0.0 aF4];
point_forces(5, 1, 1:3) = [Fz5 Fx5 aF5];
point_forces(6, 1, 1:3) = [Fz6 0.0 aF6];
point_forces(7, 1, 1:3) = [0.0 0.0  L7];


# Node forces
# ===========
# Node forces in global coordinates X and Z.
# nid - Node id.
# fid - Force id, one node can have multiple forces.
# Fx - Force in global X direction, kN.
# Fz - Force in global Z direction, kN.
# My - Moment around global Y axis, kNm.
#         nid fid        [ Fx  Fz  My]
node_forces = zeros(number_of_nodes, 1, 3);
node_forces(1, 1, 1:3) = [0.0 0.0 0.0];


# Support shifts
# ==============
# Support shift multiplied by scaling multiplier
# nid - Support node id.
# sid - Shift id.
# dx - Shift in global coordinate X direction, multiplied by scale.
# dz - Shift in global coordinate Z direction, multiplied by scale.
# fi - Revolving around global Y axis.
# Hinges: 1 - fixed, 0 - open/free
#           nid  sid       [ dx  dz  fi]
support_shift = zeros(number_of_nodes, 1, 3);
support_shift(1, 1, 1:3) = [0.0 0.0 0.0];


disp('')
disp('=============================')
disp('     Node coordinates        ')
disp('  nid       X        Z       ')
disp('-----------------------------')
for nid = 1 : node_count
    disp(sprintf('  %2i     %7.4f  %7.4f', nid, node_coordinates(nid, 1:2)))
endfor
disp('-----------------------------')

if number_of_nodes != node_count
    display(sprintf('Number of nodes must be %d, but is %d', number_of_nodes, node_count))
    error(' Faulty node data.')
endif


disp('===============================================')
disp('        Topology and  hinges                   ')
disp('  End of the element, beginning of the element,')
disp('  axial-, shear-, moment hinge; 1 - fixed hinge')
disp('-----------------------------------------------')
topology_and_hinges = element_topology_hinges(:, 4:11)


# Displacements and force numbers at the end and the beginning of element.
for eid = 1 : element_count
    displacements_forces_numbers(eid, :) = [eid*12 - 11 : eid*12];
endfor

# Combine all element information into one matrix.
#  1  2  3  4  5  6  7  8  9 10 11 12 EIp EAp GAp 2 1; # element 1
# 13 14 15 16 17 18 19 20 21 22 23 24 EIr EAr GAr 4 2; # element 2
# 25 26 27 28 29 30 31 32 33 34 35 36 EIp EAp GAp 4 3; # element 3
element_properties = [displacements_forces_numbers element_topology_hinges];


disp('===========================================================================')
disp(' Topology of elements                                          Moment hinge')
disp('   |    At the end    | At the beginning |                       | nid | ^ ')
disp('eid| u  w fi  N  Q  M | u  w fi  N  Q  M |   EI     EA      GA   | s  e|s e')
disp('---------------------------------------------------------------------------')
cols = [1:17 20 23];
sformat = '%2i  %2i %2i %2i %2i %2i %2i  %2i %2i %2i %2i %2i %2i  %7.1e %7.1e %7.1e %2i %2i %1i %1i';
for eid = 1:element_count
    disp(sprintf(sformat, eid, element_properties(eid, cols)))
endfor
disp('---------------------------------------------------------------------------')
disp('')


# Checking data
if number_of_elements != element_count
    disp(sprintf('Number of elements must  %d, but is %d', number_of_elements, element_count))
    error('Faulty elements data')
endif

disp('=========================================')
disp(' Restrictions on support displacements   ')
disp('  sid  nid   u  w fi                     ')
disp('-----------------------------------------')
for sid = 1 : size(support_nodes, 1)
    disp(sprintf('  %2i   %2i   %2i %2i %2i %2i %2i', sid, support_nodes(sid, 1:4)))
endfor
disp('-----------------------------------------')

# Checking data
support_reactions_found = sum(sum(support_nodes(:, 2:4)))
if support_reactions_count != support_reactions_found
    disp(sprintf('Number of support reactions must be %d, but is %d', support_reactions_count, support_reactions_found))
    error('Faulty support reactions data')
endif


# Do the actual calculation.
AlgPar = LaheFrameDFIm(scale, support_reactions_count, distributed_loads, point_forces, node_forces, support_nodes, support_shift, node_coordinates, element_properties);


disp('')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' Element displacements and forces determined by transfer matrix             ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
load_variant
disp('')

# Calculate lengths for all elements.
element_lengths = elementLengths(node_coordinates, element_properties);

# Names for rows types.
row_names = [ 'displace u'
              'displace w'
              'rotat   fi'
              'normal f N'
              'shear fo Q'
              'moment f M'];

# The element displacements and forces
for eid = 1 : element_count
    [EI, EA, GAr] = num2cell(element_properties(eid, 13:15)){:};
    Li = element_lengths(eid);
    x_step = Li / divisions; % element is divided into parts
    force = point_forces(eid, :, :);
    dist_load = distributed_loads(eid, :, :);
    AP = AlgPar(eid, :)';

    x = 0;
    for j = 1 : divisions + 1 # 5 - displacements and forces at x = 0.0
        vvF = ylfhlin(1.0, x, EA, GAr, EI);
        vB = ESTFrKrmus(1.0, x, Li, force, dist_load, EA, EI);
        Fvv(:, j) = vvF * AP + vB;
        x += x_step;
    endfor

    # Output element displacements and forces
    disp(sprintf('Element %i (l=%5.3f m)', eid, Li))
    for i = 1:3
        disp(sprintf(' %10s   %9.1e   %9.1e   %9.1e   %9.1e  %9.1e', row_names(i, :), Fvv(i, 1:5)))
    endfor
    for i = 4:6
        disp(sprintf(' %10s %9.3f   %9.3f   %9.3f   %9.3f  %9.3f', row_names(i, :), Fvv(i, 1:5)))
    endfor
    disp('----------------------------------------------------------------------')
endfor


disp('  ')
disp(' Testing for static equilibrium:  ')
disp('  ')
disp(' Support_reactions = [X(85) X(86) X(87) X(88) X(89) X(90) X(91) X(92)] ')
Support_reactions = [X(85) X(86) X(87) X(88) X(89) X(90) X(91) X(92)]
disp('The_lengths=[l1 l2 h08 h1  h] ')
The_lengths = [l l h08 h1 h]
disp(' q3=p1 or q3=p2  ')
q3 = q
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
disp('')
disp('  Sum of the moments acting  about point 5:  ')
disp(' sumM5=-p3*h08*h08/2+q3*(l1+1.0)*((l1+1.0)/2+l2)-F1*0.6*h... ')
disp('       +F3*(0.6*l2)+F2*0.5*h+X(87)+X(90)+l1*X(89)+(l1+l2)*X(86) ')
sumM5 = - p3 * h08 * h08 / 2 + q3 * (l + 1.0) * ((l + 1.0) / 2 + l) ...
  - F1 * 0.6 * h + F3 * (0.6 * l) + F2 * 0.5 * h + X(87) + X(90) + l * X(89) + (l + l) * X(86)
disp('')
disp(' Calculations verified the static equilibrium of the frame   ')



# Drawing figures
# ===============

# First we need to extract some data.
for eid = 1 : element_count
    siireVardaA(eid, 1:3) = element_properties(eid, 7:9);
    siireVardaL(eid, 1:3) = element_properties(eid, 1:3);
    joudVardaA(eid, 1:3) = element_properties(eid, 10:12);
    joudVardaL(eid, 1:3) = element_properties(eid, 4:6);

    LkoordN = element_properties(eid, 16);
    AkoordN = element_properties(eid, 17);
    DeltaX(eid) = node_coordinates(LkoordN, 1) - node_coordinates(AkoordN, 1);
    DeltaZ(eid) = node_coordinates(LkoordN, 2) - node_coordinates(AkoordN, 2);
    VGRx(eid, 1) = node_coordinates(AkoordN, 1);
    VGRx(eid, 2) = node_coordinates(LkoordN, 1);
    VGRz(eid, 1) = node_coordinates(AkoordN, 2);
    VGRz(eid, 2) = node_coordinates(LkoordN, 2);
endfor

# Draw the first figure, picture of frame.
figure(1)
hold on
ax = gca();
set(ax, 'linewidth', 1);
grid off
d = 4.0;
axis([-d Lp+d -H-8.0 0.5*d], 'ij')

# Draw elements as lines
for eid = 1:element_count
    plot(VGRx(eid, :), VGRz(eid, :), '3', 'linewidth', 5)
endfor

divs = 2;
title('spESTframe1DefWFI', 'fontsize', 12)
ax = gca();
xlabel('x', 'fontsize', 12)
ylabel('z', 'fontsize', 12)
xticks(1, :) = - d:divs:Lp + d;
set (ax, 'xtick', xticks)

# Numbering nodes
for nid = 1:node_count;
    text(node_coordinates(nid, 1), node_coordinates(nid, 2) + 0.25, sprintf('%i', nid))
endfor

# Numbering elements
for eid = 1:element_count;
    text(mean(VGRx(eid, 1:2)) + 0.1, mean(VGRz(eid, 1:2)) - 0.2, sprintf('%i', eid))
endfor

# Writing title and displacement and force numbers for elements
x1 = -1;
x2 = 7;
y1 = -H - 6;
y2 = -H - 6;
step = 0.5;

text(x1, y1 - 1.5*step, 'Numeration of displacements and forces', 'fontweight', 'bold')
text(x1, y1, '  u   w  fi   N   Q   M', 'fontsize', 9)
text(x2, y2, '  u   w  fi   N   Q   M', 'fontsize', 9)
text(x1-0.4, y1+3.5, 'at the beginning', 'rotation', 90)
text(x2-0.4, y2+2.8, 'at the end', 'rotation', 90)

for eid = 1:element_count;
    y1 += step;
    y2 += step;
    str = sprintf('%3i %3i %3i %3i %3i %3i', siireVardaA(eid, 1:3), joudVardaA(eid, 1:3));
    text(x1, y1, str, 'fontsize', 9)
    str = sprintf('%3i %3i %3i %3i %3i %3i', siireVardaL(eid, 1:3), joudVardaL(eid, 1:3));
    text(x2, y2, str, 'fontsize', 9)
endfor
y1 += 2*step;
text(x1, y1, 'Support reactions: 85 86 87 88 89 90 91 92 ')

# http://www.gnu.org/software/octave/doc/interpreter/Printing-Plots.html
print('-dpng', 'spESTframe1DefWFI.png');
print('-dfig', '-landscape', '-mono', '-solid', 'spESTframe1DefWFI.fig');
refresh

figure(2)
spy(spA)
title('spy(spA) - The non-zero elements of sparse matrix spA() [3%]')
refresh
print('-dpng', 'spESTframe1DefWFI_sparse_matrix.png');
print('-dfig', '-landscape', 'spESTframe1DefWFI_sparse_matrix.fig');
refresh

disp('')
disp(' The END ')

