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

##disp(' The compiler of a basic-, compatibility-, joint equilibrium-, ')
##disp(' side conditions- and support restrtictions sparse equations.  ')
##disp(' Solving system of sparse equations. Find the initial parameter vectors ')
##disp(' for elements displacements, forces and the support reactions.  ')
##disp(' OUTPUT: AlgParm -- the initial parameter vector for elements. ')
##
## baasi0 - scaling multiplier for the displacements,
## Ntoerkt - the number of the support reactions,
## esQkoormus - the uniformly distributed load in local coordinate z direction
##              esQkoormu(LoadsqONelement,3,ElementideArv),
## esFjoud - the point load in local coordinate z direction
##           esFjoud(LoadsF_on_Element,2,ElementideArv),
## sNodeF - the node forces in global  coordinates
##          sNodeF(2,1,NodeedeArv),
## support_nodes - the restrtictions on the support displacements
## tSiire - the support shift
##          tSiire(2,1,NodeedeArv)
## coordinates - the nodal coordinates
## element_properties - the topology of elements
##         element_properties(ElementideArv,1:16)

function AlgPar = LaheFrameDFIm(baasi0, support_reactions_count, esQkoormus, esFjoud, sNodeF, support_nodes, tSiire, coordinates, element_properties)
if nargin != 9
    error('Function LaheBeamDFIm() has wrong number of input arguments!')
end

global spA
global X
global siireVardaA

element_count = size(element_properties)(1);
node_count = size(coordinates)(1);
NNK = 12 * element_count + support_reactions_count;
support_nodes_count = size(support_nodes)(1);

number_of_support_reactions = 0;
for i = 1:support_nodes_count
    for j = 2:3
        if support_nodes(i, j) == 1
            number_of_support_reactions += 1;
        endif
    endfor
endfor

lvarras = VardaPikkus(node_count, element_count, coordinates, element_properties);

# Nr; daf''s at the end: u, w, fi, N, Q, M; Node number; axial -, shear -, moment hinge.
# disp(' Nr;  daf''s at the end:  u, w, fi, N, Q, M; Node number; axial-, shear-, moment hinge ')

for i = 1:element_count
    # Beginning
    elemvabNode(i, 1) = i;
    elemvabNode(i, 2:7) = element_properties(i, 1:6);
    elemvabNode(i, 8) = element_properties(i, 16);
    elemvabNode(i, 9:11) = element_properties(i, 18:20);
    # End
    elemvabNode(element_count+i, 1) = i;
    elemvabNode(element_count+i, 2:7) = element_properties(i, 7:12);
    elemvabNode(element_count+i, 8) = element_properties(i, 17);
    elemvabNode(element_count+i, 9:11) = element_properties(i, 21:23);
endfor

disp('================================================================================')
disp('  Element number;     DaF numbers;    Node number;   axial, shear, moment hinge ')
disp('--------------------------------------------------------------------------------')
disp(' Element number   u w fi N Q M        Node  N    Q    M  hinge -  true=1        ')
disp('--------------------------------------------------------------------------------')
AB = sortrows(elemvabNode(:, 8));
[!, idx] = sortrows(elemvabNode(:, 8));
ABB = elemvabNode(idx, :)

# Define a sparse matrix spA for the left side of the equation.
spA = sparse(NNK, NNK);

# Define an empty vector B, the right-hand side of the equation.
B = zeros(NNK, 1);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('----- Writing basic equations of frame with the transfer matrix ---- ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
# http://digi.lib.ttu.ee/opik_eme/Ehitusmehaanika.pdf#page=397

for i = 1:element_count
    EI = element_properties(i, 13);
    EA = element_properties(i, 14);
    GAr = element_properties(i, 15);
    Li = lvarras(i, 1);
    xx = Li;

    Fjoud = esFjoud(:, 1:3, i);
    %       Fz    Fx     a
    %Fjoud=[0.0    0.0   0.0;
    qkoormus = esQkoormus(:, 1:4, i);
    %       qz     qx     qA      qL
    %qkoormus=[0.0   0.0   0.0  0.0   0.0;

    # --------- The transfer matrix equation --------
    spvF = ysplvfmhvI(baasi0, xx, Li, EA, GAr, EI);
    vB = ESTFrKrmus(baasi0, Li, Li, Fjoud, qkoormus, EA, EI);
    IIv = i * 6 - 5;
    IJv = i * 12 - 11;
    spA = spInsertBtoA(spA, IIv, IJv, spvF);
    B = addBtoA(B, vB, IIv, 1);
endfor

equation_count = size(spA);
non_zero_elements_in_the_basic_equations = nnz(spA);

disp('-------- Basic equations are inserted into spA --------')
disp(sprintf('Rows of basic equations: %d', equation_count(1)))
disp(sprintf('Cols of basic equations: %d', equation_count(2)))
disp(sprintf('Nonzero elements in basic equations: %d', nnz(spA)))
disp('  ')
Nr1 = equation_count(1) + 1;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Compatibility equations of displacements at nodes ')
From_row = sprintf('Compatibility equations begin from row: %d', Nr1)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

for i = 1:node_count
    Node = i;
    AVS = VardadSolmes(node_count, element_count, i, AB, ABB);
    AVSdim = size(AVS);

    if AVSdim(1) > 1
        for j = 2:AVSdim(1)
            Varras1 = AVS(1, 1);
            VarrasN = AVS(j, 1);

            AVSsum = AVS(1, 11) + AVS(j, 11);

            if AVSsum == 0
                unit = SpTeisendusMaatriks(node_count, element_count, Varras1, coordinates, element_properties);
                spA = spInsertBtoA(spA, Nr1, AVS(1, 2), unit);
                unit = SpTeisendusMaatriks(node_count, element_count, VarrasN, coordinates, element_properties);
                spA = spInsertBtoA(spA, Nr1, AVS(j, 2), -unit);
                Nr1 += 3;
            else
                unit = SpTeisendusMaatriks2x2(node_count, element_count, Varras1, coordinates, element_properties);
                spA = spInsertBtoA(spA, Nr1, AVS(1, 2), unit);
                unit = SpTeisendusMaatriks2x2(node_count, element_count, VarrasN, coordinates, element_properties);
                spA = spInsertBtoA(spA, Nr1, AVS(j, 2), -unit);
                Nr1 += 2;
            endif
        endfor
    endif
endfor

disp('-----')
spA_rowsD = size(spA)(1)
nnzD = nnz(spA);
disp('  ')
disp('----- Compatibility equations of displacements are inserted into spA  ----  ')
disp('  ')
compatibility_equations_rows = spA_rowsD - equation_count(2)
non_zero_elements_in_compatibility_equations = nnzD - non_zero_elements_in_the_basic_equations
disp('  ')
# %
toemuutuja = equation_count(2) + 1;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Joint equilibrium equations at nodes ')
From_rows = sprintf('Joint equilibrium equations begin from row: %d', Nr1)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

support_nodes_extension = zeros(node_count, 4);
for i = 1:node_count
    support_nodes_extension(i, 1) = i;
    for j = 1:support_nodes_count;
        if support_nodes_extension(i, 1) == support_nodes(j, 1)
            support_nodes_extension(i, 2:4) = support_nodes(j, 2:4);
        endif
    endfor
endfor

for i = 1:node_count
    Node = i;
    AVS = VardadSolmes(node_count, element_count, Node, AB, ABB);
    uvfiLaiend(Node, 1) = sum(support_nodes_extension(Node, 2:4));

    Nr1 = size(spA)(1) + 1;
    uvfi = 0;

    AVSdim = size(AVS);
    if AVSdim(1) == 1
        if uvfiLaiend(Node, 1) == 0
            VarrasN = AVS(1, 1);
            unit = SpTeisendusMaatriks(node_count, element_count, VarrasN, coordinates, element_properties);
            spA = spInsertBtoA(spA, Nr1, AVS(1, 5), unit);
            Arv123 = 3;
            B(Nr1 : Nr1+2, 1) = sNodeF(1:3, 1, Node);
            Nr1 += 3;
        endif

        for k = 1 : size(support_nodes)(1)
            if support_nodes(k, 1) == Node
                uvfi = sum(support_nodes(k, 2:4));
                Nr1 = size(spA)(1) + 1;

                VarrasN = AVS(1, 1);
                if uvfi == 3
                    unit = SpTeisendusMaatriks(node_count, element_count, VarrasN, coordinates, element_properties);
                    spA = spInsertBtoA(spA, Nr1, AVS(1, 5), unit);
                    spA = spInsertBtoA(spA, Nr1, toemuutuja, -sparse(eye(uvfi)));
                    Nr1 += uvfi;
                    toemuutuja += uvfi;
                elseif uvfi == 2
                    unit = SpTeisendusMaatriks2x2(node_count, element_count, VarrasN, coordinates, element_properties);
                    spA = spInsertBtoA(spA, Nr1, AVS(1, 5), unit);
                    spA = spInsertBtoA(spA, Nr1, toemuutuja, -sparse(eye(uvfi)));
                    Nr1 += uvfi;
                    toemuutuja += uvfi;
                elseif uvfi == 1
                    TugiX = support_nodes(k, 2);
                    TugiZ = support_nodes(k, 3);
                    if TugiX == 1
                        unit = SpToeReaktsioonXvektor(node_count, element_count, VarrasN, coordinates, element_properties);
                        spA = spInsertBtoA(spA, Nr1, AVS(1, 5), unit);
                        Nr1 += uvfi
                        toemuutuja += uvfi;
                    endif

                    if TugiZ == 1
                        unit = SpToeReaktsioonZvektor(node_count, element_count, VarrasN, coordinates, element_properties);
                        spA = spInsertBtoA(spA, Nr1, AVS(1, 5), unit);
                        spA = spSisestaArv(spA, Nr1, toemuutuja, 1);
                        Nr1 += uvfi
                        toemuutuja += uvfi;
                    endif
                endif
            endif
        endfor
    endif

    if AVSdim(1) > 1
        for j = 1:AVSdim(1)
            VarrasN = AVS(j, 1);
            AVStas = AVS(j, 11);
            b1 = Nr1
            b2 = Nr1 + 2;
            Arv123 = 3;

            if AVStas == 0
                unit = SpTeisendusMaatriks(node_count, element_count, VarrasN, coordinates, element_properties);
            else
                unit = SpTeisendusMaatriks2x2(node_count, element_count, VarrasN, coordinates, element_properties);
                b2 -= 1;
                Arv123 -= 1;
            endif
            spA = spInsertBtoA(spA, Nr1, AVS(j, 5), unit);

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
            transmatrix = -SpTeisendusUhikMaatriks0x1v(1);
            spA = spInsertBtoA(spA, Nr1, toemuutuja, transmatrix);
            Nr1 += 1;
            toemuutuja += 1;
        elseif has_supports != 0
            transmatrix = -sparse(eye(has_supports));
            spA = spInsertBtoA(spA, Nr1, toemuutuja, transmatrix);
            Nr1 += has_supports;
            toemuutuja += has_supports;
        endif

    B(b1:b2, 1) = sNodeF(1:Arv123, 1, Node);
 
    Nr1 = size(spA)(1) + 1;
    endif
endfor

spA_rowsE = size(spA)(1);
nnzE = nnz(spA);
disp('  ')
disp('-----Equilibrium equations are inserted into spA  ----  ')
disp('  ')
equilibrium_equations_rows = spA_rowsE - spA_rowsD
non_zero_elements_in_equilibrium_equations = nnzE - nnzD

Nr1 = size(spA)(1) + 1;
ABBdim = size(ABB);
#
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Side conditions (hinges) ')
From_rows = sprintf('Side conditions begin from row: %d', Nr1)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

for i = 1:ABBdim(1)
    if ABB(i, 10) == 1
        spA = spSisestaArv(spA, Nr1, ABB(i, 6), 1);
        Nr1 += 1;
    endif
    if ABB(i, 11) == 1
        spA = spSisestaArv(spA, Nr1, ABB(i, 7), 1);
        Nr1 += 1;
    endif
endfor

spA_rows = size(spA)(1);

nnzN = nnz(spA);
disp('  ')
disp('----- Side conditions are inserted into spA  ----  ')
disp('  ')
side_condition_rows = spA_rows - spA_rowsE
non_zero_elements_in_side_condition = nnzN - nnzE

Nr1 = size(spA)(1) + 1;

toereaktsioonid_algavad = Nr1;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Restrictions on support displacements ')
From_rows = sprintf('Restrictions on support displacements begin from row: %d', Nr1)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

support_nodes_count = size(support_nodes)(1);

number_of_support_reactions = sum(sum(support_nodes(:, 2:4)))
support_reactions = zeros(number_of_support_reactions, 6);

ToeNodees = zeros(support_nodes_count, 8);
ToeNodees(:, 1:4) = support_nodes(:, 1:4);

for i = 1:support_nodes_count
    support_node_number = support_nodes(i, 1);
    for j = 1:ABBdim(1)
        if support_node_number == ABB(j, 8)
            ToeNodees(i, 5:8) = ABB(j, 1:4);
            break;
        endif
    endfor
endfor

# Support node
ToeSiireRealNR = toereaktsioonid_algavad - 1;
nmr = 1;
for j = 1:support_nodes_count
    for k = 2:4
        if support_nodes(j, k) == 1
            ToeSiireRealNR += 1;
            support_reactions(nmr, 1) = ToeSiireRealNR;
            support_reactions(nmr, 2) = ToeNodees(j, 5);
            support_reactions(nmr, 3) = k - 1;
            support_reactions(nmr, 4:6) = ToeNodees(j, 6:8);
            support_reactions(nmr, 7) = ToeNodees(j, 1);
            nmr += 1;
        endif
    endfor
endfor

# Support reactions
for i = 1:number_of_support_reactions
    n = support_reactions(i, 1);
    VarrasS = support_reactions(i, 2);
    UWFi = support_reactions(i, 3);
    NodeA = support_reactions(i, 4);
    NodeB = support_reactions(i, 7);
    shift = tSiire(UWFi, 1, NodeB);

    switch (UWFi);
        case{1}
            SpTv = SpToeSiirdeUvektor(node_count, element_count, VarrasS, coordinates, element_properties);
            B(n, 1) = shift;
            spA = spInsertBtoA(spA, n, NodeA, SpTv);
            d = 'x';
        case{2}
            SpTv = SpToeSiirdeWvektor(node_count, element_count, VarrasS, coordinates, element_properties);
            spA = spInsertBtoA(spA, n, NodeA, SpTv);
            d = 'z';
        case{3}
            SpTV = sparse(1, 3, 1);
            spA = spInsertBtoA(spA, n, NodeA, SpTV);
            d = 'y';
    endswitch
    disp(sprintf('Support shift tSiire(%i, 1, %i) = %i in %s direction at the node %i', UWFi, NodeB, shift, d, NodeB))
endfor

nnzR = nnz(spA);
disp('  ')
disp('----- Restriction equations are inserted into spA  ----  ')
disp('  ')
non_zero_elements_in_restrtictions_equations = nnzR - nnzN

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

Bsuurus = size(B)(1);

disp(' No    B ')
disp('  ')
for i = 1:Bsuurus
    disp(sprintf('%2i %7.3f', i, B(i, 1)))
endfor
disp('  ')

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Solving system of sparse equations. ')
disp(' X = spA \ B; ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

% Solving of the main matrix. Many interesting. Much wow.
% spA * X = B, solve for X
X = spA \ B;

disp(' No     X ')
disp('  ')
for i = 1:Bsuurus
    disp(sprintf('%2i   %9.3e', i, X(i, 1)))
endfor
disp('  ')

display(sprintf('Support reactions begin from row X: %d', toereaktsioonid_algavad))
disp('  ')

for i = 1:element_count
    elemendiN(i, 1) = i;
    siireVardaA(i, 1:3) = element_properties(i, 7:9);
    siireVardaL(i, 1:3) = element_properties(i, 1:3);
    joudVardaA(i, 1:3) = element_properties(i, 10:12);
    joudVardaL(i, 1:3) = element_properties(i, 4:6);
endfor


disp('=======================================================')
disp('     Support displacements/(Cx, Cz, CMy) at nodes  ')
disp('   X_No,  Element_No,  u/w/fi_No, u_No,  w_No,  fi_No ')
disp('-------------------------------------------------------')
for i = 1:number_of_support_reactions
    disp(sprintf('     %2i     %2i         %2i        %3i    %3i     %3i  ', support_reactions(i, 1:6)))
endfor
disp('==================================================================')
disp(' Displacements and force numbers at beginning of the element  ')
disp('   No,    u,   w,   fi     N,   Q,   M  ')
disp('------------------------------------------------------------------')
for i = 1:element_count
    disp(sprintf('  %2i     %2i   %2i   %2i     %2i   %2i   %2i  ', i, siireVardaA(i, 1:3), joudVardaA(i, 1:3)))
endfor

disp('-------------------------------------------------------')
disp(' ')
disp('==================================================================')
disp(' Displacements and forces numbers at end of the element  ')
disp('   No,    u,   w,   fi     N,   Q,   M  ')
disp('------------------------------------------------------------------')
for i = 1:element_count
    disp(sprintf('  %2i     %2i   %2i   %2i     %2i   %2i   %2i ', i, siireVardaL(i, 1:3), joudVardaL(i, 1:3)))
endfor
#
disp('-------------------------------------------------------')
disp(' ')
#
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Initial parameter vectors for elements displacements and forces  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' ')
disp('-------- Descaling multiplier for the displacements = 1/baasi0 --------  ')
disp(' ')

AlgPar = zeros(element_count, 6)
for i = 1:element_count
    displacement = X(siireVardaA(i, 1:3))' / baasi0;
    force = X(joudVardaA(i, 1:3), 1)';
    SiireJoud(2*i-1, :) = [displacement force];
    AlgPar(i, :) = [displacement force];
    displacement = X(siireVardaL(i, 1:3))' / baasi0;
    force = X(joudVardaL(i, 1:3), 1)';
    SiireJoud(2*i, :) = [displacement force];
endfor

disp('============================================================================')
disp(' Unscaled initial parameter vector ')
disp('Element No    u          w          fi             N          Q          M ')
disp('----------------------------------------------------------------------------')
for i = 1:element_count
    disp(sprintf('  %2i   %9.2e   %9.2e   %9.2e    %9.2f  %9.2f  %9.2f', i, AlgPar(i, 1:6)))
endfor
disp('----------------------------------------------------------------------------')

for i = 1:element_count
    JrN2(2*i - 1) = i;
    JrN2(2*i) = i;
endfor

#SiirdeJoud
disp('============================================================================')
disp('      Unscaled daf-s at beginning/end of elements ')
disp('Element No    u          w          fi             N          Q          M  ')
disp('----------------------------------------------------------------------------')
for i = 1 : element_count*2
    disp(sprintf('  %2i   %9.2e   %9.2e   %9.2e    %9.2f  %9.2f  %9.2f', JrN2(i), SiireJoud(i, :)))
endfor
disp('----------------------------------------------------------------------------')
disp('')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
Support_reactions = sprintf('Support reactions begin from row X: %d', toereaktsioonid_algavad)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp('  ')
disp('======================================================================')
disp(' No         X          Node   (Cx <=> 1 /Cz <=> x / Cy(moment) <=> 3) ')
disp('----------------------------------------------------------------------')
disp('')
for i = toereaktsioonid_algavad:Bsuurus
    toenr = i - toereaktsioonid_algavad + 1
    disp(sprintf('%3i   %+12.6e    %3i            %3i ', i, X(i, 1), support_reactions(toenr, 7), support_reactions(toenr, 3)))
endfor
disp('--------------------------------------------------------------')
disp('  ')

AlgPar;
endfunction

