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

##disp('==================================================================')
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
## sSolmF - the node forces in global  coordinates
##          sSolmF(2,1,SolmedeArv),
## support_nodes - the restrtictions on the support displacements
## tSiire - the support shift
##          tSiire(2,1,SolmedeArv)
## coordinates - the nodal coordinates
## element_properties - the topology of elements
##         element_properties(ElementideArv,1:16)
##
function AlgParm = LaheFrameDFIm(baasi0, support_reactions_count, esQkoormus, esFjoud, sSolmF, support_nodes, tSiire, coordinates, element_properties)
if nargin != 9
    error('Function LaheBeamDFI has wrong number of input arguments!')
end

global spA
global X
global siireVardaA
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Nr; daf''s at the end: u, w, fi, N, Q, M; Node number; axial -, shear -, moment hinge.
# disp(' Nr;  daf''s at the end:  u, w, fi, N, Q, M; Node number; axial-, shear-, moment hinge ')

for i = 1:element_count
    # Beginning
    elemvabsolm(i, 1) = i;
    elemvabsolm(i, 2:7) = element_properties(i, 1:6);
    elemvabsolm(i, 8) = [element_properties(i, 16)];
    elemvabsolm(i, 9:11) = [element_properties(i, 18:20)];
    # End
    elemvabsolm(element_count+i, 1) = i;
    elemvabsolm(element_count+i, 2:7) = element_properties(i, 7:12);
    elemvabsolm(element_count+i, 8) = element_properties(i, 17);
    elemvabsolm(element_count+i, 9:11) = element_properties(i, 21:23);
endfor

disp('================================================================================')
disp('  Element number;     DaF numbers;    Node number;   axial, shear, moment hinge ')
disp('--------------------------------------------------------------------------------')
disp(' Element number   u w fi N Q M        Node  N    Q    M  hinge -  true=1        ')
disp('--------------------------------------------------------------------------------')
AB = sortrows(elemvabsolm(:, 8));
[!, idx] = sortrows(elemvabsolm(:, 8));
ABB = elemvabsolm(idx, :)

# --------- Define a matrix spA to be sparse -----------
disp('-------- Sparse matrix instantiation --------  ')
spA = sparse(NNK, NNK); % nullistan võrrandisüsteemi kordajate maatriksi
# --------- Define a vector B ------------
disp('-------- Right-hand side of the equations (RHS). --------  ')
B = zeros(NNK, 1); % nullistan võrrandisüsteemi vabaliikmed
%%%%%%%%%%%%%%%%%%%%%%%%%

# ------------ The basic equations of frame ----------------
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp('----- Writing basic equations of frame with the transfer matrix ----  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
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
rows = sprintf('rows_of_basic_equations: %d', equation_count(1))
col = sprintf('cols_of_basic_equations: %d', equation_count(2))
spA_nnz = sprintf('non_zero_elements_in_the_basic_equations: %d', nnz(spA))
nnzB = spA_nnz;
disp('  ')
Nr1 = equation_count(1) + 1;
#
#disp('  ')
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VarrasteArv = element_count;
SolmedeArv = node_count;
# %%%%
Vord = '=';
Kom = ',';
Semi = ';';
#Tyhk = '  ';
Arv1 = num2str(1);
Arv2 = num2str(2);
Arv3 = num2str(3);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Compatibility equations of displacements at nodes ')
From_row = sprintf('Compatibility equations begin from row: %d', Nr1)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')

for i = 1:node_count
    Solm = i;
    AVS = VardadSolmes(node_count, element_count, i, AB, ABB);
    AVSdim = size(AVS);

    if AVSdim(1) > 1
        for j = 2:AVSdim(1)
            Varras1 = AVS(1, 1);
            VarrasN = AVS(j, 1);

            AVSsum = AVS(1, 11) + AVS(j, 11);

            if AVSsum == 0
                SpTM3x3 = SpTeisendusMaatriks(node_count, element_count, Varras1, coordinates, element_properties);
                spA = spInsertBtoA(spA, Nr1, AVS(1, 2), SpTM3x3);

                SpTM3x3 = SpTeisendusMaatriks(node_count, element_count, VarrasN, coordinates, element_properties);
                spA = spInsertBtoA(spA, Nr1, AVS(j, 2), -SpTM3x3);
                Nr1 += 3;
            else
                SpTM2x2 = SpTeisendusMaatriks2x2(node_count, element_count, Varras1, coordinates, element_properties);
                spA = spInsertBtoA(spA, Nr1, AVS(1, 2), SpTM2x2);
                SpTM2x2 = SpTeisendusMaatriks2x2(node_count, element_count, VarrasN, coordinates, element_properties);
                spA = spInsertBtoA(spA, Nr1, AVS(j, 2), -SpTM2x2);
                Nr1 += 2;
            endif
        endfor
    endif
endfor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----- ')
vorrandeidD = size(spA);
#
disp(' ')
spA_rows = vorrandeidD(1, 1)
spA_rowsD = spA_rows;
spA_cols = vorrandeidD(1, 2)
#
nnzD = nnz(spA);
spA_nnz = sprintf('non_zero_elements_in_spA: %d', nnz(spA))
disp('  ')
disp('----- Compatibility equations of displacements are inserted into spA  ----  ')
disp('  ')
compatibility_equations_rows = spA_rowsD - equation_count(2)
non_zero_elements_in_compatibility_equations = nnzD - non_zero_elements_in_the_basic_equations
#nnzB = nnz(spA)
disp('  ')
# %
# %pohivorrandiXe
toemuutuja = equation_count(2) + 1;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Joint equilibrium equations at nodes ')
From_rows = sprintf('Joint equilibrium equations begin from row: %d', Nr1)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
#
str1 = num2str(Nr1);
strSolm = num2str(Solm);
koolonV = ':';

support_nodesLaiend = zeros(node_count, 4);
for i = 1:node_count
    support_nodesLaiend(i, 1) = i;
    for j = 1:support_nodes_count;
        if (support_nodesLaiend(i, 1) - support_nodes(j, 1)) == 0
            support_nodesLaiend(i, 2:4) = support_nodes(j, 2:4);
        endif
    endfor
endfor
# support_nodes
# support_nodesLaiend
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:node_count % sõlmede arv
    # uvfiLaiend(i) = support_nodesLaiend(i, 2) + support_nodesLaiend(i, 3) + support_nodesLaiend(i, 4);
    disp('---------------------------------------------------------')
    Solm = i;
    Node = Solm
    AVS = VardadSolmes(node_count, element_count, Solm, AB, ABB);
    uvfiLaiend(Node, 1) = support_nodesLaiend(Node, 2) + support_nodesLaiend(Node, 3) + support_nodesLaiend(Node, 4);
    #
    disp('---------------------------------------------------------')
    vorrandeid = size(spA);
    spRidu = vorrandeid(1, 1);
    Nr1 = spRidu + 1;
    uvfi = 0;
    # uvfi = support_nodes(k, 2) + support_nodes(k, 3) + support_nodes(k, 4); # uvfi - the number of reactions
    AVSdim = size(AVS);
    if AVSdim(1) == 1
        if uvfiLaiend(Node, 1) == 0
            VarrasN = AVS(1, 1);
            SpTM3x3 = SpTeisendusMaatriks(node_count, element_count, VarrasN, coordinates, element_properties);
            spA = spInsertBtoA(spA, Nr1, AVS(1, 5), SpTM3x3);
            Arv123 = 3;
            B(Nr1 : Nr1+2, 1) = sSolmF(1:3, 1, Solm);
            Nr1 = Nr1 + 3;
        endif
        support_nodes_count = size(support_nodes)(1);
        for k = 1:support_nodes_count
            if (support_nodes(k, 1) - Solm) == 0
                uvfi = support_nodes(k, 2) + support_nodes(k, 3) + support_nodes(k, 4);
                Number_of_reaction = sprintf('The number of reaction at node (node_nr, reactions):  %d, %d', Node, uvfi)
                vorrandeid = size(spA);
                spRidu = vorrandeid(1, 1);
                Nr1 = spRidu + 1;
                str1 = num2str(Nr1);
                strSolm = num2str(Solm);
                koolonV = ':';
                AVStas = AVS(1, 11);
                if (uvfi - 3) == 0
                    VarrasN = AVS(1, 1);
                    SpTM3x3 = SpTeisendusMaatriks(node_count, element_count, VarrasN, coordinates, element_properties);
                    spA = spInsertBtoA(spA, Nr1, AVS(1, 5), SpTM3x3);
                    spA = spInsertBtoA(spA, Nr1, toemuutuja, -sparse(eye(3)));
                    Nr1 += 3;
                    toemuutuja += 3;
                endif
             
                if (uvfi - 2) == 0
                    VarrasN = AVS(1, 1);
                    SpTM2x2 = SpTeisendusMaatriks2x2(node_count, element_count, VarrasN, coordinates, element_properties);
                    spA = spInsertBtoA(spA, Nr1, AVS(1, 5), SpTM2x2);
                    spA = spInsertBtoA(spA, Nr1, toemuutuja, -sparse(eye(2)));
                    Nr1 += 2;
                    toemuutuja += 2;
                endif
                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (uvfi - 1) == 0
                    VarrasN = AVS(1, 1);
                    TugiX = support_nodes(k, 2);
                    TugiZ = support_nodes(k, 3);
                    if (TugiX - 1) == 0
                        SpTM1xX = SpToeReaktsioonXvektor(node_count, element_count, VarrasN, coordinates, element_properties);
                        spA = spInsertBtoA(spA, Nr1, AVS(1, 5), SpTM1xX);
                        Nr1 += 1
                        toemuutuja += 1;
                    endif

                    if (TugiZ - 1) == 0
                        SpTM1xZ = SpToeReaktsioonZvektor(node_count, element_count, VarrasN, coordinates, element_properties);
                        spA = spInsertBtoA(spA, Nr1, AVS(1, 5), SpTM1xZ);
                        spA = spSisestaArv(spA, Nr1, toemuutuja, 1);
                        Nr1 += 1
                        toemuutuja += 1;
                    endif
                endif
            endif
        endfor
    endif

    toemuutujaS = toemuutuja - 1;
    toemuutujaid = toemuutuja - equation_count(2) - 1;

    if AVSdim(1) > 1
        k = 1;
        taskVorreid(1, k) = 0;
     
        for j = 1:AVSdim(1)
            VarrasN = AVS(j, 1);
            str1 = num2str(Nr1);
            AVS(1, 11);
            AVS(j, 11);
            AVStas = AVS(j, 11);

            if AVStas == 0
                SpTM3x3 = SpTeisendusMaatriks(node_count, element_count, VarrasN, coordinates, element_properties);
                spA = spInsertBtoA(spA, Nr1, AVS(j, 5), SpTM3x3);
                k += 1;
                taskVorreid(1, k) = 3;
             
                b1 = Nr1
                b2 = Nr1 + 2;
                Arv123 = 3;
                ontugesid = 0;
                Kas_on_toe_solm = support_nodesLaiend(Node, 2) + support_nodesLaiend(Node, 3) + support_nodesLaiend(Node, 4);
                if (Kas_on_toe_solm != 0)
                    if ((uvfiLaiend(Node, 1)) - 3) == 0
                        ontugesid = 3;
                    endif
                    TSMLD1 = support_nodesLaiend(Node, 2) + support_nodesLaiend(Node, 3);
                    if (TSMLD1 - 2) == 0
                        ontugesid = 2;
                    endif
                    if (TSMLD1 - 1) == 0
                        if ((support_nodesLaiend(Node, 2)) - 1) == 0
                            ontugesid = 1;
                        endif
                        if ((support_nodesLaiend(Node, 3)) - 1) == 0
                            ontugesid = 11;
                        endif
                    endif
                endif
            else
                SpTM2x2 = SpTeisendusMaatriks2x2(node_count, element_count, VarrasN, coordinates, element_properties);
                spA = spInsertBtoA(spA, Nr1, AVS(j, 5), SpTM2x2);
             
                k += 1;
                taskVorreid(1, k) = 2;
             
                b1 = Nr1;
                b2 = Nr1 + 1;
                Arv123 = 2;
                ontugesid = 0;
                Kas_on_toe_solm = support_nodesLaiend(Node, 2) + support_nodesLaiend(Node, 3) + support_nodesLaiend(Node, 4);
                Kas_on = Kas_on_toe_solm;
                if (Kas_on_toe_solm != 0)
                    if (uvfiLaiend(Node) - 3) == 0
                        ontugesid = 3;
                    endif
                    TSMLD1 = support_nodesLaiend(Node, 2) + support_nodesLaiend(Node, 3);
                    if (TSMLD1 - 2) == 0
                        ontugesid = 2;
                    endif
                    if (TSMLD1 - 1) == 0
                        if ((support_nodesLaiend(Node, 2)) - 1) == 0
                            ontugesid = 1;
                        endif
                        if ((support_nodesLaiend(Node, 3)) - 1) == 0
                            ontugesid = 11;
                        endif
                    endif
                endif
            endif
        endfor

        if ontugesid == 11
            transmatrix = -SpTeisendusUhikMaatriks0x1v(1);
            spA = spInsertBtoA(spA, Nr1, toemuutuja, transmatrix);
            Nr1 += 1;
            toemuutuja += 1;
        elseif ontugesid != 0
            transmatrix = -sparse(eye(ontugesid));
            spA = spInsertBtoA(spA, Nr1, toemuutuja, transmatrix);
            Nr1 += ontugesid;
            toemuutuja += ontugesid;
        endif

    B(b1:b2, 1) = sSolmF(1:Arv123, 1, Solm);
 
    vorrandeid = size(spA);
    spRidu = vorrandeid(1);
    Nr1 = spRidu + 1;
    endif
endfor
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vorrandeid = size(spA);
spRidu = vorrandeid(1);
vorrandeidE = size(spA);
spA_rows = vorrandeidE(1, 1)
spA_rowsE = spA_rows;
spA_cols = vorrandeidE(1, 2)
nnzE = nnz(spA);
spA_nnz = sprintf('non_zero_elements_in_spA: %d', nnz(spA))
disp('  ')
disp('-----Equilibrium equations are inserted into spA  ----  ')
disp('  ')
equilibrium_equations_rows = spA_rowsE - spA_rowsD
non_zero_elements_in_equilibrium_equations = nnzE - nnzD
%nnzE=nnz(spA)
disp('  ')
vorrandeid = size(spA);
spRidu = vorrandeid(1, 1);
Nr1 = spRidu + 1;
ABB = ABB;
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
Nr1 -= 1;
#
vorrandeidN = size(spA);
spA_rows = vorrandeidN(1);
spA_cols = vorrandeidN(2);

nnzN = nnz(spA);
spA_nnz = sprintf('non_zero_elements_in_spA: %d', nnz(spA))
disp('  ')
disp('----- Side conditions are inserted into spA  ----  ')
disp('  ')
side_condition_rows = spA_rows - spA_rowsE
non_zero_elements_in_side_condition = nnzN - nnzE
#nnzN = nnz(spA)
disp('  ')
# %
#disp('----- ')
#
vorrandeid = size(spA);
spRidu = vorrandeid(1, 1);
Nr1 = spRidu + 1;
#
ABB;
support_nodes = support_nodes;
TSARV = size(support_nodes);
support_nodes_count = TSARV(1, 1);
toereaktsioonid_algavad = Nr1;
#
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Restrictions on support displacements ')
From_rows = sprintf('Restrictions on support displacements begin from row: %d', Nr1)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
#
toesiirdedN = 1;
TSARV = size(support_nodes);
support_nodes_count = TSARV(1, 1);
ABBMoot = ABBdim(1, 1);
number_of_support_reactions = 0;
ToeSiireRealNR = toereaktsioonid_algavad - 1;

for i = 1:support_nodes_count
    tsolNr = support_nodes(i, 1);
    for j = 2:4
        if (support_nodes(i, j) == 1)
            number_of_support_reactions = number_of_support_reactions + 1;
        endif
    endfor
endfor
number_of_support_reactions;
disp('   ')
ToeSidemed = zeros(number_of_support_reactions, 6);
ToeSolmes = zeros(support_nodes_count, 8);
ToeSolmes(:, 1:4) = support_nodes(:, 1:4);

for i = 1:support_nodes_count
    tsolNr = support_nodes(i, 1);
    LeidsinVarda = 0;
    for j = 1:ABBMoot
        if (LeidsinVarda == 0)
            if (tsolNr - ABB(j, 8)) == 0
                ToeSolmes(i, 5:8) = ABB(j, 1:4);
                LeidsinVarda = 1;
            endif
        endif
    endfor
endfor
# %%%%%%%%%
#ToeSolmes
# %%%%%%%%
nmr = 0;
for j = 1:support_nodes_count
    for k = 2:4
        KusOnYks = support_nodes(j, k);
        if (KusOnYks - 1) == 0
            nmr = nmr + 1;
            ToeSiireRealNR = ToeSiireRealNR + 1;
            ToeSidemed(nmr, 1) = ToeSiireRealNR;
            ToeSidemed(nmr, 2) = ToeSolmes(j, 5);
            ToeSidemed(nmr, 3) = k - 1;
            ToeSidemed(nmr, 4:6) = ToeSolmes(j, 6:8);
            ToeSidemed(nmr, 7) = ToeSolmes(j, 1);
        endif
    endfor
endfor
# %%%%%%%%%%%%%
#ToeSidemed
# %%%%%%%%%%%%%
for i = 1:number_of_support_reactions
    n = ToeSidemed(i, 1);
    UWFi = ToeSidemed(i, 3);

    VarrasS = ToeSidemed(i, 2);
    SolmA = ToeSidemed(i, 4);
    SolmB = ToeSidemed(i, 7);
    shift = tSiire(UWFi, 1, SolmB);
    
    switch (UWFi);
        case{1}
            SpTv = SpToeSiirdeUvektor(node_count, element_count, VarrasS, coordinates, element_properties);
            B(n, 1) = shift;
            spA = spInsertBtoA(spA, n, SolmA, SpTv);
            d = 'x';
        case{2}
            SpTv = SpToeSiirdeWvektor(node_count, element_count, VarrasS, coordinates, element_properties);
            spA = spInsertBtoA(spA, n, SolmA, SpTv);
            d = 'z';
        case{3}
            SpTV = sparse(1, 3, 1);
            spA = spInsertBtoA(spA, n, SolmA, SpTV);
            d = 'y';
    endswitch
    disp(sprintf('Support shift tSiire(%i, 1, %i) = %i in %s direction at the node %i', UWFi, SolmB, shift, d, SolmB))
endfor

spRidu = size(spA)(1);
vorrandeidR = size(spA);
spA_rows = vorrandeidR(1, 1)
spA_rowsR = spA_rows;
spA_cols = vorrandeidR(1, 2)
#
nnzR = nnz(spA);
spA_nnz = sprintf('non_zero_elements_in_spA: %d', nnz(spA))
disp('  ')
disp('----- Restriction equations are inserted into spA  ----  ')
disp('  ')
restrtictions_equations_rows = spA_rowsR - spA_rows
non_zero_elements_in_restrtictions_equations = nnzR - nnzN
#
disp('  ')
#
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
#
spA_rank = sprank(spA);
spA;
#
B;
BsuurusS = size(B);
Bsuurus = BsuurusS(1, 1);

#
disp(' No    B ')
disp('  ')
for i = 1:Bsuurus
    # %9.3e %9.3f
    vorrandiNr = i;
    disp(sprintf('%3i %7.3f', vorrandiNr, B(i, 1)))
endfor
disp('  ')
#
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Solving system of sparse equations. ')
disp(' X=spA\B; ')
#Support_reactions = sprintf('are the last rows: %d', restrtictions_equations_rows)
Support_reactions = sprintf('Support reactions begin from row X: %d', toereaktsioonid_algavad)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
#
%%%%%%%%%%%%%%%%%%%%
X = spA \ B; % lahend  võrrandisõteemist spA*X=B
%%%%%%%%%%%%%%%%%%%%%
#
disp(' No     X ')
disp('  ')
for i = 1:Bsuurus
    # %9.3e   %9.3f
    vorrandiNr = i;
    disp(sprintf('%3i   %9.3e', vorrandiNr, X(i, 1)))
endfor
#
disp('  ')
Support_reactions = sprintf('Support reactions begin from row X: %d', toereaktsioonid_algavad)
#
disp('  ')
# %disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
for i = 1:element_count
    elemendiN(i, 1) = i;
    siireVardaA(i, 1:3) = element_properties(i, 7:9);
    siireVardaL(i, 1:3) = element_properties(i, 1:3);
    joudVardaA(i, 1:3) = element_properties(i, 10:12);
    joudVardaL(i, 1:3) = element_properties(i, 4:6);
endfor
#
siireVardaA;
siireVardaL;
joudVardaA;
joudVardaL;
# %%%%%%%%%%%%
for i = 1:element_count
    JrN(i, 1) = i;
endfor
#
JrNT = JrN';
#
disp('=======================================================')
disp('     Support displacements/(Cx, Cz, CMy) at nodes  ')
disp('   X_No,  Element_No,  u/w/fi_No, u_No,  w_No,  fi_No ')
disp('-------------------------------------------------------')
for i = 1:number_of_support_reactions
    disp(sprintf('     %2i     %2i         %2i        %3i    %3i     %3i  ', ToeSidemed(i, 1), ToeSidemed(i, 2), ToeSidemed(i, 3), ToeSidemed(i, 4), ToeSidemed(i, 5), ToeSidemed(i, 6)))
endfor
disp('-------------------------------------------------------')
#
#The_displacement_DaF - s_numbers_at_beginning_of_the_element = sprintf('elementNr, displacement: %d', JrN(i, 1) siireVardaA)
disp('==================================================================')
disp(' Displacements and force numbers at beginning of the element  ')
disp('   No,    u,   w,   fi     N,   Q,   M  ')
disp('------------------------------------------------------------------')
for i = 1:element_count
    disp(sprintf('  %2i     %2i   %2i   %2i     %2i   %2i   %2i  ', JrN(i, 1), siireVardaA(i, 1), siireVardaA(i, 2), siireVardaA(i, 3), joudVardaA(i, 1), joudVardaA(i, 2), joudVardaA(i, 3)))
endfor
%
disp('-------------------------------------------------------')
disp(' ')
disp('==================================================================')
disp(' Displacements and forces numbers at end of the element  ')
disp('   No,    u,   w,   fi     N,   Q,   M  ')
disp('------------------------------------------------------------------')
for i = 1:element_count
    disp(sprintf('  %2i     %2i   %2i   %2i     %2i   %2i   %2i ', JrN(i, 1), siireVardaL(i, 1), siireVardaL(i, 2), siireVardaL(i, 3), joudVardaL(i, 1), joudVardaL(i, 2), joudVardaL(i, 3)))
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
#
for i = 1:element_count
    # %siireVardaA(i,1:3)
    AlgPar(i, 1:6) = [X(siireVardaA(i, 1), 1) / baasi0 X(siireVardaA(i, 2), 1) / baasi0 X(siireVardaA(i, 3), 1) / baasi0 ...
    X(joudVardaA(i, 1), 1) X(joudVardaA(i, 2), 1) X(joudVardaA(i, 3), 1)];
endfor
#
disp('============================================================================')
disp(' Unscaled initial parameter vector ')
disp('Element No    u          w          fi             N          Q          M ')
disp('----------------------------------------------------------------------------')
for i = 1:element_count
    disp(sprintf('  %2i   %9.2e   %9.2e   %9.2e    %9.2f  %9.2f  %9.2f', JrN(i, 1), AlgPar(i, 1), AlgPar(i, 2), AlgPar(i, 3), AlgPar(i, 4), AlgPar(i, 5), AlgPar(i, 6)))
    #AlgPar
endfor
disp('----------------------------------------------------------------------------')
#
JJ = 0;
for i = 1:element_count
    #siireVardaL(i, 1:3)
    JJ = JJ + 1;
    SiireJoud(JJ, 1:6) = [X(siireVardaA(i, 1), 1) / baasi0 X(siireVardaA(i, 2), 1) / baasi0 X(siireVardaA(i, 3), 1) / baasi0 ...
      X(joudVardaA(i, 1), 1) X(joudVardaA(i, 2), 1) X(joudVardaA(i, 3), 1)];
    JJ = JJ + 1;
    SiireJoud(JJ, 1:6) = [X(siireVardaL(i, 1), 1) / baasi0 X(siireVardaL(i, 2), 1) / baasi0 X(siireVardaL(i, 3), 1) / baasi0 ...
      X(joudVardaL(i, 1), 1) X(joudVardaL(i, 2), 1) X(joudVardaL(i, 3), 1)];
endfor
#
element_count2 = element_count * 2;
Nn = 0;
Nn2 = 0;
for i = 1:element_count
    Nn = Nn + 1;
    Nn2 = Nn2 + 1;
    JrN2(Nn, 1) = Nn2;
    Nn = Nn + 1;
    JrN2(Nn, 1) = Nn2;
endfor
#
#JrN2
#SiireJoud
#
disp('============================================================================')
disp('      Unscaled daf-s at beginning/end of elements ')
disp('Element No    u          w          fi             N          Q          M ')
disp('----------------------------------------------------------------------------')
for i = 1:element_count2
    disp(sprintf('  %2i   %9.2e   %9.2e   %9.2e    %9.2f  %9.2f  %9.2f', JrN2(i, 1), SiireJoud(i, 1), SiireJoud(i, 2), SiireJoud(i, 3), SiireJoud(i, 4), SiireJoud(i, 5), SiireJoud(i, 6)))
endfor
disp('----------------------------------------------------------------------------')

#disp('============================================================================')
#
disp('  ')
#disp('============================================================================')
#
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
##disp('The support reactions, toereaktsioonid_algavad  ')
Support_reactions = sprintf('Support reactions begin from row X: %d', toereaktsioonid_algavad)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp('  ')
#ToeSidemed 6 ja 3
#ToeSidemed
disp('======================================================================')
disp(' No         X          Node   (Cx <=> 1 /Cz <=> x / Cy(moment) <=> 3) ')
disp('----------------------------------------------------------------------')
disp('  ')
for i = toereaktsioonid_algavad:Bsuurus
    # %9.3e   %9.3f
    vorrandiNr = i; %%9.3f
    disp(sprintf('%3i   %+12.6e    %3i            %3i ', vorrandiNr, X(i, 1), ToeSidemed(i - toereaktsioonid_algavad + 1, 7), ToeSidemed(i - toereaktsioonid_algavad + 1, 3)))
endfor
#
disp('--------------------------------------------------------------')
disp('  ')

AlgParm = AlgPar;
endfunction

