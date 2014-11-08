##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-07-16
## Last modified: Andres Lahe, 2013-07-21
##                Mattias Põldaru, 2014-11-07
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

%======================================================================
%> @file LaheFramedfim.m
%> @brief Function to compile all conditions for the equation matrix.
%======================================================================
%> @brief Function to compile all conditions for the equation matrix.
%>
%> The compiler of a basic-, compatibility-, joint equilibrium-,
%> side conditions- and support restrtictions sparse equations.
%> Solving system of sparse equations. Find the initial parameter vectors
%> for elements displacements, forces and the support reactions.
%>
%> @param scale Scaling multiplier for the displacements.
%> @param number_of_support_reactions The number of the support reactions.
%> @param distributed_loads The uniformly distributed load in local
%>            coordinate z direction. [loadvalue, 3, element_count]
%> @param point_forces The point load in local coordinate z direction.
%>            [loadvalue, 2, element_count]
%> @param node_forces The node forces in global coordinates.
%>            [2, 1, node_count]
%> @param support_nodes The restrtictions on the support displacements.
%> @param support_shift The support shift. [2, 1, count_of_nodes]
%> @param coordinates The nodal coordinates.
%> @param element_properties The topology of elements. [elem_id, 1:16]
%>
%> @retval The initial parameter vector for elements.
%======================================================================
function AlgPar = LaheFrameDFIm(scale, support_reactions_count, distributed_loads, point_forces, node_forces, support_nodes, support_shift, coordinates, element_properties)
if nargin != 9
    error('Function LaheBeamDFIm() has wrong number of input arguments!')
end

global spA
global X
global siireVardaA

element_count = size(element_properties, 1);
node_count = size(coordinates, 1);
NNK = 12 * element_count + support_reactions_count;
support_nodes_count = size(support_nodes, 1);

number_of_support_reactions = 0;
for id = 1:support_nodes_count
    for j = 2:3
        if support_nodes(id, j) == 1
            number_of_support_reactions += 1;
        endif
    endfor
endfor

lvarras = VardaPikkus(node_count, element_count, coordinates, element_properties);

# Nr; daf''s at the end: u, w, fi, N, Q, M; Node number; axial -, shear -, moment hinge.
# disp(' Nr;  daf''s at the end:  u, w, fi, N, Q, M; Node number; axial-, shear-, moment hinge ')

for eid = 1:element_count
    # Beginning
    element_free_node(eid, 1) = eid;
    element_free_node(eid, 2:7) = element_properties(eid, 1:6);
    element_free_node(eid, 8) = element_properties(eid, 16);
    element_free_node(eid, 9:11) = element_properties(eid, 18:20);
    # End
    element_free_node(element_count+eid, 1) = eid;
    element_free_node(element_count+eid, 2:7) = element_properties(eid, 7:12);
    element_free_node(element_count+eid, 8) = element_properties(eid, 17);
    element_free_node(element_count+eid, 9:11) = element_properties(eid, 21:23);
endfor

disp('================================================================================')
disp(' Element number;     DaF numbers;    Node number;   axial, shear, moment hinge  ')
disp('--------------------------------------------------------------------------------')
disp(' Element number   u w fi N Q M        Node  N    Q    M  hinge -  true=1        ')
disp('--------------------------------------------------------------------------------')
AB = sortrows(element_free_node(:, 8));
[!, idx] = sortrows(element_free_node(:, 8));
ABB = element_free_node(idx, :)

# Define a sparse matrix spA for the left side of the equation.
spA = sparse(NNK, NNK);

# Define an empty vector B, the right-hand side of the equation.
B = zeros(NNK, 1);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('----- Writing basic equations of frame with the transfer matrix ---- ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
# http://digi.lib.ttu.ee/opik_eme/Ehitusmehaanika.pdf#page=397

for eid = 1:element_count
    EI = element_properties(eid, 13);
    EA = element_properties(eid, 14);
    GAr = element_properties(eid, 15);
    Li = lvarras(eid, 1);
    xx = Li;

    # Force matrix [Fx Fx a]
    forces = point_forces(eid, :, 1:3);
    
    # Load matrix [qz qx qA qL]
    loads = distributed_loads(eid, :, 1:4);

    # The transfer matrix equation
    spvF = ysplvfmhvI(scale, xx, Li, EA, GAr, EI);
    vB = ESTFrKrmus(scale, Li, Li, forces, loads, EA, EI);
    position_row = eid * 6 - 5;
    position_col = eid * 12 - 11;
    spA = spInsertBtoA(spA, position_row, position_col, spvF);
    B = addBtoA(B, vB, position_row, 1);
endfor


equation_count = size(spA);
non_zero_elements_in_the_basic_equations = nnz(spA);

disp('-------- Basic equations are inserted into spA --------')
disp(sprintf('Rows of basic equations: %d', equation_count(1)))
disp(sprintf('Cols of basic equations: %d', equation_count(2)))
disp(sprintf('Nonzero elements in basic equations: %d', nnz(spA)))
disp('')

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Compatibility equations of displacements at nodes ')
row = equation_count(1) + 1;
From_row = sprintf('Compatibility equations begin from row: %d', row)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

for nid = 1:node_count
    Node = nid;
    node_elems = VardadSolmes(node_count, element_count, nid, AB, ABB);
    node_elemsdim = size(node_elems, 1);

    if node_elemsdim > 1
        for eid = 2:node_elemsdim
            element_first = node_elems(1, 1);
            element_n = node_elems(eid, 1);

            node_elemssum = node_elems(1, 11) + node_elems(eid, 11);

            if node_elemssum == 0
                transform = spTransformation(3, element_first, coordinates, element_properties);
                spA = spInsertBtoA(spA, row, node_elems(1, 2), transform);
                transform = spTransformation(3, element_n, coordinates, element_properties);
                spA = spInsertBtoA(spA, row, node_elems(eid, 2), -transform);
                row += 3;
            else
                transform = spTransformation(2, element_first, coordinates, element_properties);
                spA = spInsertBtoA(spA, row, node_elems(1, 2), transform);
                transform = spTransformation(2, element_n, coordinates, element_properties);
                spA = spInsertBtoA(spA, row, node_elems(eid, 2), -transform);
                row += 2;
            endif
        endfor
    endif
endfor

spA_rowsD = size(spA, 1);
nnzD = nnz(spA);
disp('----- Compatibility equations of displacements are inserted into spA  ----  ')
disp('')
compatibility_equations_rows = spA_rowsD - equation_count(2)
non_zero_elements_in_compatibility_equations = nnzD - non_zero_elements_in_the_basic_equations
disp('')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Joint equilibrium equations at nodes ')
From_rows = sprintf('Joint equilibrium equations begin from row: %d', row)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

support_nodes_extension = zeros(node_count, 4);
for nid = 1:node_count
    support_nodes_extension(nid, 1) = nid;
    for sid = 1:support_nodes_count;
        if support_nodes_extension(nid, 1) == support_nodes(sid, 1)
            support_nodes_extension(nid, 2:4) = support_nodes(sid, 2:4);
        endif
    endfor
endfor

toemuutuja = equation_count(2) + 1;

for nid = 1:node_count
    node_elems = VardadSolmes(node_count, element_count, nid, AB, ABB);
    uvfiLaiend(nid, 1) = sum(support_nodes_extension(nid, 2:4));

    row = size(spA, 1) + 1;
    uvfi = 0;

    node_elemsdim = size(node_elems, 1);
    if node_elemsdim == 1
        if uvfiLaiend(nid, 1) == 0
            element_n = node_elems(1, 1);
            transform = spTransformation(3, element_n, coordinates, element_properties);
            spA = spInsertBtoA(spA, row, node_elems(1, 5), transform);
            Arv123 = 3;
            B(row : row+2) = node_forces(nid, 1, 1:3);
            row += 3;
        endif

        for k = 1 : size(support_nodes, 1)
            if support_nodes(k, 1) == nid
                uvfi = sum(support_nodes(k, 2:4));
                row = size(spA, 1) + 1;

                element_n = node_elems(1, 1);
                if uvfi == 3
                    transform = spTransformation(3, element_n, coordinates, element_properties);
                    spA = spInsertBtoA(spA, row, node_elems(1, 5), transform);
                    spA = spInsertBtoA(spA, row, toemuutuja, -speye(uvfi));
                    row += uvfi;
                    toemuutuja += uvfi;
                elseif uvfi == 2
                    transform = spTransformation(2, element_n, coordinates, element_properties);
                    spA = spInsertBtoA(spA, row, node_elems(1, 5), transform);
                    spA = spInsertBtoA(spA, row, toemuutuja, -speye(uvfi));
                    row += uvfi;
                    toemuutuja += uvfi;
                elseif uvfi == 1
                    TugiX = support_nodes(k, 2);
                    TugiZ = support_nodes(k, 3);
                    if TugiX == 1
                        reaction_x_vector = SpToeReaktsioonXvektor(node_count, element_count, element_n, coordinates, element_properties);
                        spA = spInsertBtoA(spA, row, node_elems(1, 5), reaction_x_vector);
                        row += uvfi
                        toemuutuja += uvfi;
                    elseif TugiZ == 1
                        reaction_z_vector = SpToeReaktsioonZvektor(node_count, element_count, element_n, coordinates, element_properties);
                        spA = spInsertBtoA(spA, row, node_elems(1, 5), reaction_z_vector);
                        spA = spSisestaArv(spA, row, toemuutuja, 1);
                        row += uvfi
                        toemuutuja += uvfi;
                    endif
                endif
            endif
        endfor
    endif

    if node_elemsdim > 1
        for eid = 1:node_elemsdim
            element_n = node_elems(eid, 1);
            node_elemstas = node_elems(eid, 11);
            b1 = row
            b2 = row + 2;
            Arv123 = 3;

            if node_elemstas == 0
                transform = spTransformation(3, element_n, coordinates, element_properties);
            else
                transform = spTransformation(2, element_n, coordinates, element_properties);
                b2 -= 1;
                Arv123 -= 1;
            endif
            spA = spInsertBtoA(spA, row, node_elems(eid, 5), transform);

            has_supports = 0;
            is_support_node = sum(support_nodes_extension(Node, 2:4));
            if is_support_node
                if uvfiLaiend(Node, 1) == 3
                    has_supports = 3;
                endif
                TSMLD1 = sum(support_nodes_extension(Node, 2:3));
                if TSMLD1 ==2
                    has_supports = 2;
                elseif TSMLD1 == 1
                    if support_nodes_extension(Node, 3) == 1
                        has_supports = 11;
                    elseif support_nodes_extension(Node, 2) == 1
                        has_supports = 1;
                    endif
                endif
            endif
        endfor

        if has_supports == 11
            transmatrix = -speye(1);
            spA = spInsertBtoA(spA, row, toemuutuja, transmatrix);
            row += 1;
            toemuutuja += 1;
        elseif has_supports != 0
            transmatrix = -speye(has_supports);
            spA = spInsertBtoA(spA, row, toemuutuja, transmatrix);
            row += has_supports;
            toemuutuja += has_supports;
        endif

    B(b1:b2) = node_forces(Node, 1, 1:Arv123);
 
    row = size(spA, 1) + 1;
    endif
endfor

spA_rowsE = size(spA, 1);
nnzE = nnz(spA);
disp('')
disp('-----Equilibrium equations are inserted into spA  ----  ')
disp('')
equilibrium_equations_rows = spA_rowsE - spA_rowsD
non_zero_elements_in_equilibrium_equations = nnzE - nnzD

row = size(spA, 1) + 1;
ABBdim = size(ABB);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Side conditions (hinges) ')
From_rows = sprintf('Side conditions begin from row: %d', row)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

for i = 1:ABBdim(1)
    if ABB(i, 10) == 1
        spA = spSisestaArv(spA, row, ABB(i, 6), 1);
        row += 1;
    endif
    if ABB(i, 11) == 1
        spA = spSisestaArv(spA, row, ABB(i, 7), 1);
        row += 1;
    endif
endfor

nnzN = nnz(spA);
disp('')
disp('----- Side conditions are inserted into spA  ----  ')
disp('')
side_condition_rows = size(spA, 1) - spA_rowsE
non_zero_elements_in_side_condition = nnzN - nnzE

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Restrictions on support displacements ')
support_reactions_start = size(spA, 1) + 1;
From_rows = sprintf('Restrictions on support displacements begin from row: %d', support_reactions_start)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

support_nodes_count = size(support_nodes, 1);

number_of_support_reactions = sum(sum(support_nodes(:, 2:4)))
support_reactions = zeros(number_of_support_reactions, 6);

in_support_node = zeros(support_nodes_count, 8);
in_support_node(:, 1:4) = support_nodes(:, 1:4);

for i = 1:support_nodes_count
    support_node_number = support_nodes(i, 1);
    for j = 1:ABBdim(1)
        if support_node_number == ABB(j, 8)
            in_support_node(i, 5:8) = ABB(j, 1:4);
            break;
        endif
    endfor
endfor

# Support node
ToeSiireRealNR = support_reactions_start - 1;
nmr = 1;
for j = 1:support_nodes_count
    for k = 2:4
        if support_nodes(j, k) == 1
            ToeSiireRealNR += 1;
            support_reactions(nmr, 1) = ToeSiireRealNR;
            support_reactions(nmr, 2) = in_support_node(j, 5);
            support_reactions(nmr, 3) = k - 1;
            support_reactions(nmr, 4:6) = in_support_node(j, 6:8);
            support_reactions(nmr, 7) = in_support_node(j, 1);
            nmr += 1;
        endif
    endfor
endfor

# Support reactions
for i = 1:number_of_support_reactions
    n = support_reactions(i, 1);
    eid = support_reactions(i, 2);
    UWFi = support_reactions(i, 3);
    NodeA = support_reactions(i, 4);
    NodeB = support_reactions(i, 7);
    shift = support_shift(NodeB, 1, UWFi);

    switch (UWFi);
        case{1}
            sup_displ_vector_U = SpToeSiirdeUvektor(node_count, element_count, eid, coordinates, element_properties);
            spA = spInsertBtoA(spA, n, NodeA, sup_displ_vector_U);
            B(n) = shift;
            d = 'x';
        case{2}
            sup_displ_vector_W = SpToeSiirdeWvektor(node_count, element_count, eid, coordinates, element_properties);
            spA = spInsertBtoA(spA, n, NodeA, sup_displ_vector_W);
            d = 'z';
        case{3}
            sup_displ_vector = sparse(1, 3, 1);
            spA = spInsertBtoA(spA, n, NodeA, sup_displ_vector);
            d = 'y';
    endswitch
    disp(sprintf('support_shift(%i, 1, %i) = %i in %s direction at the node %i', UWFi, NodeB, shift, d, NodeB))
endfor

nnzR = nnz(spA);
disp('')
disp('----- Restriction equations are inserted into spA  ----  ')
disp('')
non_zero_elements_in_restrtictions_equations = nnzR - nnzN

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

B_size = size(B, 1);

disp(' No     B   ')
disp('------------')
for i = 1:B_size
    b = B(i);
    if b != 0
        disp(sprintf('%3i %8.3f', i, b))
    endif
endfor
disp('Note: only nonzero b values are shown.')
disp('')

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Solving system of sparse equations. ')
disp(' X = spA \ B; ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

% Solving of the main matrix. Many interesting. Much wow.
% spA * X = B, solve for X
X = spA \ B;

disp('   No       X     ')
disp('------------------')
for i = 1:B_size
    disp(sprintf('  %3i   %10.3e', i, X(i, 1)))
endfor
disp('')

display(sprintf('Support reactions begin from row X: %d', support_reactions_start))
disp('')

for eid = 1:element_count
    siireVardaA(eid, 1:3) = element_properties(eid, 7:9);
    siireVardaL(eid, 1:3) = element_properties(eid, 1:3);
    joudVardaA(eid, 1:3) = element_properties(eid, 10:12);
    joudVardaL(eid, 1:3) = element_properties(eid, 4:6);
endfor


disp('=======================================================')
disp('     Support displacements/(Cx, Cz, CMy) at nodes      ')
disp('    X_No    eid     u/w/fi_No    u_No   w_No   fi_No   ')
disp('-------------------------------------------------------')
for rid = 1:number_of_support_reactions
    disp(sprintf('     %2i     %2i         %2i        %3i    %3i     %3i  ', support_reactions(rid, 1:6)))
endfor
disp('-------------------------------------------------------')
disp('')
disp('=================================================================')
disp(' Displacements and force numbers at the beginning of the element ')
disp('  eid     u    w   fi      N    Q    M                           ')
disp('-----------------------------------------------------------------')
for eid = 1:element_count
    disp(sprintf('  %2i     %2i   %2i   %2i     %2i   %2i   %2i  ', eid, siireVardaA(eid, 1:3), joudVardaA(eid, 1:3)))
endfor
disp('-----------------------------------------------------------')
disp('')
disp('===========================================================')
disp(' Displacements and force numbers at the end of the element ')
disp('  eid     u    w   fi      N    Q    M                     ')
disp('-----------------------------------------------------------')
for eid = 1:element_count
    disp(sprintf('  %2i     %2i   %2i   %2i     %2i   %2i   %2i ', eid, siireVardaL(eid, 1:3), joudVardaL(eid, 1:3)))
endfor
disp('-----------------------------------------------------------')
disp('')

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Initial parameter vectors for elements displacements and forces  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp('')
disp('-------- Descaling multiplier for the displacements = 1/scale --------  ')
disp('')

AlgPar = zeros(element_count, 6);
for i = 1:element_count
    displacement = X(siireVardaA(i, 1:3))' / scale;
    force = X(joudVardaA(i, 1:3), 1)';
    SiireJoud(2*i-1, :) = [displacement force];
    AlgPar(i, :) = [displacement force];
    displacement = X(siireVardaL(i, 1:3))' / scale;
    force = X(joudVardaL(i, 1:3), 1)';
    SiireJoud(2*i, :) = [displacement force];
endfor

disp('============================================================================')
disp(' Unscaled initial parameter vector ')
disp('  eid      u           w          fi              N          Q          M   ')
disp('----------------------------------------------------------------------------')
for eid = 1 : element_count
    disp(sprintf('  %2i   %9.2e   %9.2e   %9.2e    %9.2f  %9.2f  %9.2f', eid, AlgPar(eid, 1:6)))
endfor
disp('----------------------------------------------------------------------------')
disp('============================================================================')
disp('      Unscaled daf-s at beginning/end of elements                           ')
disp('  eid      u           w          fi              N          Q          M   ')
disp('----------------------------------------------------------------------------')
for i = 1 : element_count*2
    disp(sprintf('  %2i   %9.2e   %9.2e   %9.2e    %9.2f  %9.2f  %9.2f', ceil(i/2), SiireJoud(i, :)))
endfor
disp('----------------------------------------------------------------------------')
disp('')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('Support reactions begin from row X: %d', support_reactions_start))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('')
disp('======================================================================')
disp(' No         X           nid        (Cx - 1 / Cz - 2 / Cy(moment) - 3) ')
disp('----------------------------------------------------------------------')
disp('')
for i = support_reactions_start:B_size
    toenr = i - support_reactions_start + 1;
    disp(sprintf('%3i   %+12.6e    %3i            %3i ', i, X(i, 1), support_reactions(toenr, 7), support_reactions(toenr, 3)))
endfor
disp('----------------------------------------------------------------------')
disp('')

AlgPar;
endfunction

