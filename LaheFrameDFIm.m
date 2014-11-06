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
##disp(' Solving system of sparse equations. Find the initial parmeter vectors ')
##disp(' for elements displacements, forces and the support reactions.  ')
##disp(' OUTPUT: AlgParm -- the initial parmeter vector for elements. ')
##
## baasi0 - scaling multiplier for the displacements,
## Ntoerkt - the number of the support reactions,
## esQkoormus - the uniformly distributed load in local coordinate z direction
##              esQkoormu(LoadsqONelement,3,ElementideArv),
## esFjoud - the point load in local coordinate z direction
##           esFjoud(LoadsF_on_Element,2,ElementideArv),
## sSolmF - the node forces in global  coordinates
##          sSolmF(2,1,SolmedeArv),
## tsolm - the restrtictions on the support displacements
## tSiire - the support shift
##          tSiire(2,1,SolmedeArv)
## krdn - the nodal coordinates
## selem - the topology of elements
##         selem(ElementideArv,1:16)
##
function AlgParm = LaheFrameDFIm(baasi0, Ntoerkts, esQkoormus, esFjoud, sSolmF, tsolm, tSiire, krdn, selem)
#if nargin != 9
if ~ (nargin == 9)
    error(' function LaheBeamDFI has wrong number of input arguments!')
end
#
global spA
global X
global siireVardaA
#disp('---  ')
# %selem=[selemjl(:,1:23)];
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EARV = size(selem);
NEARV = EARV(1, 1);
SARV = size(krdn);
NSARV = SARV(1, 1);
NNK = 12 * NEARV + Ntoerkts;
# ==========
TSARV = size(tsolm);
NTSARV = TSARV(1, 1);
## %%%%%%%%%%%
ToeSidemeteArv = 0;
for i = 1:NTSARV
    sJrN(i) = i;
    for j = 2:3EARV = size(selem);
        NEARV = EARV(1, 1);
        SARV = size(krdn);
        NSARV = SARV(1, 1);
        NNK = 12 * NEARV + Ntoerkts;
        # ==========
        TSARV = size(tsolm);
        NTSARV = TSARV(1, 1);
        ## %%%%%%%%%%%
        ToeSidemeteArv = 0;
        for i = 1:NTSARV
            sJrN(i) = i;
            for j = 2:3
                if (tsolm(i, j) == 1)
                    ToeSidemeteArv = ToeSidemeteArv + 1;
                endif
            endfor
        endfor
        #
        if (tsolm(i, j) == 1)
            ToeSidemeteArv = ToeSidemeteArv + 1;
        endif
    endfor
endfor
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selemjl = selem;
#1
#
lvarras = VardaPikkus(NSARV, NEARV, krdn, selem);
#lvarras
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('   ')
NEARV;
NNK = 12 * NEARV + Ntoerkts;

#lvarras;
#
#selemjl;
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Nr; daf''s at the end: u, w, fi, N, Q, M; Node number; axial -, shear -, moment hinge.
#disp(' Nr;  daf''s at the end:  u, w, fi, N, Q, M; Node number; axial-, shear-, moment hinge ')
#
for i = 1:NEARV
    selemlopp(i, 1) = i;
    selemlopp(i, 2:7) = selemjl(i, 1:6);
    selemlopp(i, 8) = [selemjl(i, 16)];
    selemlopp(i, 9:11) = [selemjl(i, 18:20)];
endfor
#
selemlopp;
# %%%%%%%%%%%
# Nr; daf''s at the beginning: u, w, fi, N, Q, M; Node number; axial -, shear -, moment hinge.
#disp(' Nr;  daf''s at the beginning:  u, w, fi, N, Q, M; Node number; axial-, shear-, moment hinge ')
#
for i = 1:NEARV
    selemaluses(i, 1) = i;
    selemaluses(i, 2:7) = selemjl(i, 7:12);
    selemaluses(i, 8) = [selemjl(i, 17)];
    selemaluses(i, 9:11) = [selemjl(i, 21:23)];
endfor
#
selemaluses;
# %%%%%%%%%%%
#elemvabsolm
for i = 1:NEARV
    elemvabsolm(i, 1:11) = selemlopp(i, 1:11);
endfor
#
for i = 1:NEARV
    elemvabsolm(NEARV + i, 1:11) = selemaluses(i, 1:11);
endfor
#
elemvabsolm;
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %spA=sparse(60,60);
AB = sortrows(elemvabsolm(:, 8));
#AB;
#
##[~, idx] = sortrows(elemvabsolm(:, 8));
[tmp, idx] = sortrows(elemvabsolm(:, 8));
ABB = elemvabsolm(idx, :);
disp('================================================================================= ')
disp('  Element number;     DaF numbers;    Node number;   axial, shear, moment hinge')
disp('-------------------------------------------------------------------------------- ')
disp(' Element number   u w fi N Q M        Node  N    Q    M  hinge -  true=1         ')
disp('-------------------------------------------------------------------------------- ')
ABB = ABB

# --------- Define a matrix spA to be sparse -----------
disp('-------- Sparse matrix instantiation --------  ')
disp('  spA=sparse(NNK,NNK)  ')
spA = sparse(NNK, NNK) % nullistan võrrandisüsteemi kordajate maatriksi
# --------- Define a vector B ------------
disp('-------- Right-hand side of the equations (RHS). --------  ')
disp('  B=zeros(NNK,1);  ')
B = zeros(NNK, 1); % nullistan võrrandisüsteemi vabaliikmed
disp('   ')
%%%%%%%%%%%%%%%%%%%%%%%%%

# ------------ The basic equations of frame ----------------
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp('----- Writing basic equations of frame with the transfer matrix ----  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
# http: / / digi.lib.ttu.ee / opik_eme / Ehitusmehaanika.pdf#page = 397 page = 693
#
IIv = 0;
IJv = 0;
#
for i = 1:NEARV
    qkoormus = zeros(1, 3);
    Fjoud = zeros(1, 2);
    spvF = zeros(4, 8);
    vB = zeros(4, 1);
    vFz = zeros(4, 1);
    #Sj = zeros(4, 1);
 
    krda = i;
    EI = selem(i, 13);
    EA = selem(i, 14);
    GAr = selem(i, 15);
    Li = lvarras(i, 1);
    xx = Li;
    #SII = SI(i, 1);
    # SII = 0;
    # qx = qxZ(i, 1);
    # qz = qzZ(i, 1);
    # aLx = aLXx(i, 1);
    # Fz = FZz(i, 1);
    # Fx = FZx(i, 1);
 
    #Fjoud = esFjoud(:, 1:2, i);
    Fjoud = esFjoud(:, 1:3, i);
    %       Fz    Fx     a
    %Fjoud=[0.0    0.0   0.0;
    #qkoormus = esQkoormus(:, 1:3, i);
    qkoormus = esQkoormus(:, 1:4, i);
    %       qz     qx     qA      qL
    %qkoormus=[0.0   0.0   0.0  0.0   0.0;
 
    # --------- The transfer matrix equation --------
    spvF = ysplvfmhvI(baasi0, xx, Li, EA, GAr, EI);
    # vB = yzhqz(baasi0, Li, qx, qz, EA, EI);
    # vB = yzhqzm(baasi0, xx, Li, qx, qz, EA, EI);
    # vFz = yzfzv(baasi0, Li, aLx, Fx, Fz, EA, EI);
    # Sj = ESTSKrmus(xx, Li, Fjoud, qkoormus) % vB4=yzShqz(xx,qx4,qz4); % vFz5=yzSfzv(xx,aLx5,Fx5,Fz5)
    Sj = ESTFrKrmus(baasi0, Li, Li, Fjoud, qkoormus, EA, EI);
    # Sj = ESTFrKrmus(baasi0, xx, Li, Fjoud, qkoormus, EA, EI)
    # vB = vB + vFz;
    vB = Sj;
    IIv = krda * 6 - 5;
    IJv = krda * 12 - 11;
    spA = spInsertBtoA(spA, IIv, IJv, spvF);
    B = InsertBtoA(B, NNK, 1, IIv, 1, vB, 6, 1);
    #
endfor
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vorrandeidB = size(spA);
rows_of_basic_equations = vorrandeidB(1, 1);
cols_of_basic_equations = vorrandeidB(1, 2);
non_zero_elements_in_the_basic_equations = nnz(spA);

spRidu = vorrandeidB(1, 1);
pohivorrandiXe = vorrandeidB(1, 2);
disp('  ')
disp('--------Basic equations are inserted into spA  --------  ')
rows = sprintf('rows_of_basic_equations: %d', spRidu)
col = sprintf('cols_of_basic_equations: %d', pohivorrandiXe)
spA_nnz = sprintf('non_zero_elements_in_the_basic_equations: %d', nnz(spA))
nnzB = spA_nnz;
disp('  ')
Nr1 = spRidu + 1;
#
#disp('  ')
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VarrasteArv = NEARV;
SolmedeArv = NSARV;
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
#
for i = 1:NSARV
    #disp('---------------------------------------------------------')
    Solm = i;
    Node = Solm
    AVS = VardadSolmes(NSARV, NEARV, Solm, AB, ABB);
    #AVS
    #
    #disp('---------------------------------------------------------')
    #
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NEARV2 = NEARV * 2;
    #
    SpTM3x3 = sparse(3, 3);
    SpTM3x3m = sparse(3, 3);
    SpTM2x2 = sparse(2, 2);
    SpTM2x2m = sparse(2, 2);
    SpTM3x3xz = sparse(3, 3);
    SpTM2x2xz = sparse(2, 2);
    #
    AVSdim = size(AVS);
 
    if AVSdim(1, 1) > 1
     
        for j = 2:AVSdim(1, 1)
         
            Varras1 = AVS(1, 1);
            VarrasN = AVS(j, 1);
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            #
            str1 = num2str(Nr1);
            AVS(1, 11);
            AVS(j, 11);
            AVSsum = AVS(1, 11) + AVS(j, 11);
            #
            if AVSsum == 0
                SpTM3x3 = SpTeisendusMaatriks(NSARV, NEARV, Varras1, krdn, selem);
                TAVS2 = eval(sprintf('%3i', AVS(1, 2)));
                TVS2 = int2str(TAVS2);
                str2 = TVS2;
                str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM3x3', ')', Semi];
                cmd = sprintf(str)
                eval(cmd)
                SpTM3x3 = SpTeisendusMaatriks(NSARV, NEARV, VarrasN, krdn, selem);
                SpTM3x3m = - SpTM3x3;
                TAVS2 = eval(sprintf('%3i', AVS(j, 2)));
                TVS2 = int2str(TAVS2);
                str2 = TVS2;
                str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM3x3m', ')', Semi];
                cmd = sprintf(str)
                eval(cmd);
                Nr1 = Nr1 + 3;
            else
                #
                SpTM2x2 = SpTeisendusMaatriks2x2(NSARV, NEARV, Varras1, krdn, selem);
                # SpTM2x2;
                str1 = num2str(Nr1);
                TAVS2 = eval(sprintf('%3i', AVS(1, 2)));
                TVS2 = int2str(TAVS2);
                str2 = TVS2;
                str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM2x2', ')', Semi];
                cmd = sprintf(str)
                eval(cmd);
                SpTM2x2 = SpTeisendusMaatriks2x2(NSARV, NEARV, VarrasN, krdn, selem);
                SpTM2x2m = - SpTM2x2;
                TAVS2 = eval(sprintf('%3i', AVS(j, 2)));
                TVS2 = int2str(TAVS2);
                str2 = TVS2;
                str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM2x2m', ')', Semi];
                cmd = sprintf(str)
                eval(cmd);
                Nr1 = Nr1 + 2;
            endif
            #
        endfor # % for j=2:AVSdim(1,1)
    endif # %if AVSdim(1,1) > 1
endfor # % for i=1:NSARV
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
compatibility_equations_rows = spA_rowsD - rows_of_basic_equations
non_zero_elements_in_compatibility_equations = nnzD - non_zero_elements_in_the_basic_equations
#nnzB = nnz(spA)
disp('  ')
# %
# %pohivorrandiXe
toemuutuja = pohivorrandiXe + 1;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Joint equilibrium equations at nodes ')
From_rows = sprintf('Joint equilibrium equations begin from row: %d', Nr1)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
#
str1 = num2str(Nr1);
strSolm = num2str(Solm);
koolonV = ':';
#strSolmB
tsolmLaiend = zeros(NSARV, 4);
for i = 1:NSARV
    tsolmLaiend(i, 1) = i;
    for j = 1:NTSARV;
        if (tsolmLaiend(i, 1) - tsolm(j, 1)) == 0
            tsolmLaiend(i, 2:4) = tsolm(j, 2:4);
        endif
    endfor
endfor
# tsolm
# tsolmLaiend
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:NSARV % sõlmede arv
    # uvfiLaiend(i) = tsolmLaiend(i, 2) + tsolmLaiend(i, 3) + tsolmLaiend(i, 4);
    disp('---------------------------------------------------------')
    Solm = i;
    Node = Solm
    AVS = VardadSolmes(NSARV, NEARV, Solm, AB, ABB);
    uvfiLaiend(Node, 1) = tsolmLaiend(Node, 2) + tsolmLaiend(Node, 3) + tsolmLaiend(Node, 4);
    #
    disp('---------------------------------------------------------')
    vorrandeid = size(spA);
    spRidu = vorrandeid(1, 1);
    Nr1 = spRidu + 1;
    # %%%%%%%%%%%%%%%%%%
    NEARV2 = NEARV * 2;
    SpTM3x3 = sparse(3, 3);
    SpTM3x3m = sparse(3, 3);
    SpTM2x2 = sparse(2, 2);
    SpTM2x2m = sparse(2, 2);
    SpTM3x3xz = sparse(3, 3);
    SpTM2x2xz = sparse(2, 2);
    uvfi = 0;
    # uvfi = tsolm(k, 2) + tsolm(k, 3) + tsolm(k, 4); # uvfi - the number of reactions
    #
    AVSdim = size(AVS);
    if AVSdim(1, 1) == 1
     
        # %%%%%%%%%%%%%%%%
        if uvfiLaiend(Node, 1) == 0
            VarrasN = AVS(1, 1);
            %      AVS
            %%      AVS =2   19   20   21   22   23   24    3    0    0    0
         
            SpTM3x3 = SpTeisendusMaatriks(NSARV, NEARV, VarrasN, krdn, selem);
            %        TAVS2=eval(sprintf('%3i',AVS(1,2)));
            TAVS2 = eval(sprintf('%3i', AVS(1, 5)));
            TVS2 = int2str(TAVS2);
            str2 = TVS2;
            str1 = num2str(Nr1);
            str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM3x3', ')', Semi];
            cmd = sprintf(str)
            eval(cmd)
            Arv123 = Arv3;
            sF(:, 1) = sSolmF(:, 1, Solm);
            % siia tehtud parandus
            strB1 = num2str(Nr1);
            strB2 = num2str(Nr1 + 2);
            %  strSolmB=Arv1;
            strB = ['B(', strB1, koolonV, strB2, Kom, Arv1, ')', '=', 'sF(', Arv1, koolonV, Arv123, Kom, Arv1, ')', Semi];
            cmd = sprintf(strB)
            eval(cmd)
         
            Nr1 = Nr1 + 3;
        endif
        # %%%%%%%%%%%%%%%%
     
        TSARV = size(tsolm);
        NTSARV = TSARV(1, 1);
        for k = 1:NTSARV
            %olentoesolmes=tsolm(k,1)
            tsolm(k, 1);
            tsolSol = (tsolm(k, 1) - Solm);
            if tsolSol == 0
                uvfi = tsolm(k, 2) + tsolm(k, 3) + tsolm(k, 4);
                Number_of_reaction = sprintf('The number of reaction at node (node_nr, reactions):  %d, %d', Node, uvfi)
                vorrandeid = size(spA);
                spRidu = vorrandeid(1, 1);
                Nr1 = spRidu + 1;
                str1 = num2str(Nr1);
                strSolm = num2str(Solm);
                koolonV = ':';
                AVStas = AVS(1, 11);
                if (uvfi - 3) == 0
                    # olentoesolmes = tsolm(k, 1); # ADDED
                    VarrasN = AVS(1, 1);
                    SpTM3x3 = SpTeisendusMaatriks(NSARV, NEARV, VarrasN, krdn, selem);
                    TAVS2 = eval(sprintf('%3i', AVS(1, 5)));
                    TVS2 = int2str(TAVS2);
                    str2 = TVS2;
                    str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM3x3', ')', Semi];
                    cmd = sprintf(str)
                    eval(cmd)
                    SpTM3x3xz = SpTeisendusUhikMaatriks(1);
                    SpTM3x3xz = - SpTM3x3xz;
                    TAVS2 = toemuutuja;
                    TVS2 = int2str(TAVS2);
                    str2 = TVS2;
                    str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM3x3xz', ')', Semi];
                    cmd = sprintf(str)
                    eval(cmd)
                 
                    Nr1 = Nr1 + 3;
                    toemuutuja = toemuutuja + 3;
                 
                endif
             
                if (uvfi - 2) == 0
                    #
                    olentoesolmes = tsolm(k, 1); # ADDED
                    VarrasN = AVS(1, 1);
                    SpTM2x2 = SpTeisendusMaatriks2x2(NSARV, NEARV, VarrasN, krdn, selem);
                    #
                    str1 = num2str(Nr1);
                    TAVS2 = eval(sprintf('%3i', AVS(1, 5)));
                    TVS2 = int2str(TAVS2);
                    str2 = TVS2;
                    str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM2x2', ')', Semi];
                    cmd = sprintf(str)
                    eval(cmd)
                 
                    SpTM2x2xz = SpTeisendusUhikMaatriks2x2(1);
                    SpTM2x2xz = - SpTM2x2xz;
                    #
                    TAVS2 = toemuutuja;
                    TVS2 = int2str(TAVS2);
                    str2 = TVS2;
                    str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM2x2xz', ')', Semi];
                    cmd = sprintf(str)
                    eval(cmd)
                 
                    Nr1 = Nr1 + 2;
                    # taskVorreid(1, k) = 2;
                    toemuutuja = toemuutuja + 2;
                endif
                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (uvfi - 1) == 0
                 
                    olentoesolmes = tsolm(k, 1); # ADDED
                    VarrasN = AVS(1, 1);
                    #AVS
                    TugiX = tsolm(k, 2);
                    TugiZ = tsolm(k, 3);
                    mArv1 = num2str(1);
                    if (TugiX - 1) == 0
                        TAVS2 = eval(sprintf('%3i', AVS(1, 5)));
                        TVS2 = int2str(TAVS2);
                        str2 = TVS2
                        SpTM1xX = SpToeReaktsioonXvektor(NSARV, NEARV, VarrasN, krdn, selem);
                        str1 = num2str(Nr1)
                        str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM1xX', ')', Semi];
                        cmd = sprintf(str)
                        eval(cmd)
                        #
                        TAVS2 = toemuutuja;
                        TVS2 = int2str(TAVS2);
                        str2 = TVS2;
                     
                        Nr1 = Nr1 + 1
                        toemuutuja = toemuutuja + 1;
                    endif
                    # SIIN
                    if (TugiZ - 1) == 0
                        TAVS2 = eval(sprintf('%3i', AVS(1, 5)));
                        TVS2 = int2str(TAVS2);
                        str2 = TVS2;
                        #
                        SpTM1xZ = SpToeReaktsioonZvektor(NSARV, NEARV, VarrasN, krdn, selem);
                        str1 = num2str(Nr1);
                        str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM1xZ', ')', Semi];
                        cmd = sprintf(str)
                        eval(cmd)
                        #
                        TAVS2 = toemuutuja;
                        TVS2 = int2str(TAVS2);
                        str2 = TVS2
                        #
                        str = ['spA=spSisestaArv(spA', Kom, str1, Kom, str2, Kom, mArv1, ')', Semi];
                        cmd = sprintf(str)
                        eval(cmd)
                     
                        Nr1 = Nr1 + 1
                        toemuutuja = toemuutuja + 1;
                    endif
                endif # uvfi - 1
                # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            endif
        endfor # % for i=1:NTSARV
        # endif
    endif
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toemuutujaS = toemuutuja - 1;
    toemuutujaid = toemuutuja - pohivorrandiXe - 1;
    #
    if AVSdim(1, 1) > 1
        k = 1;
        taskVorreid(1, k) = 0;
     
        for j = 1:AVSdim(1, 1)
            VarrasN = AVS(j, 1);
            str1 = num2str(Nr1);
            AVS(1, 11);
            AVS(j, 11);
            AVStas = AVS(j, 11);
            #
            if AVStas == 0
                #
                SpTM3x3 = SpTeisendusMaatriks(NSARV, NEARV, VarrasN, krdn, selem);
                TAVS2 = eval(sprintf('%3i', AVS(j, 5)));
                TVS2 = int2str(TAVS2);
                str2 = TVS2;
                str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM3x3', ')', Semi];
                cmd = sprintf(str)
                eval(cmd)
                k = k + 1;
                taskVorreid(1, k) = 3;
             
                strB1 = str1;
                NrB2 = Nr1 + 2;
                strB2 = num2str(NrB2);
                Arv123 = Arv3;
                ontugesid = 0;
                Kas_on_toe_solm = tsolmLaiend(Node, 2) + tsolmLaiend(Node, 3) + tsolmLaiend(Node, 4);
                if (Kas_on_toe_solm != 0)
                    if ((uvfiLaiend(Node, 1)) - 3) == 0
                        ontugesid = 3;
                    endif
                    TSMLD1 = tsolmLaiend(Node, 2) + tsolmLaiend(Node, 3);
                    if (TSMLD1 - 2) == 0
                        ontugesid = 2;
                    endif
                    if (TSMLD1 - 1) == 0
                        if ((tsolmLaiend(Node, 2)) - 1) == 0
                            ontugesid = 1;
                        endif
                        if ((tsolmLaiend(Node, 3)) - 1) == 0
                            ontugesid = 11;
                        endif
                    endif
                endif
                #
            else
                #
                SpTM2x2 = SpTeisendusMaatriks2x2(NSARV, NEARV, VarrasN, krdn, selem);
                #
                str1 = num2str(Nr1);
                TAVS2 = eval(sprintf('%3i', AVS(j, 5)));
                TVS2 = int2str(TAVS2);
                str2 = TVS2;
                str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM2x2', ')', Semi];
                cmd = sprintf(str)
                eval(cmd);
             
                k = k + 1;
                taskVorreid(1, k) = 2;
             
                strB1 = str1;
                NrB2 = Nr1 + 1;
                strB2 = num2str(NrB2);
                Arv123 = Arv2;
                ontugesid = 0;
                Kas_on_toe_solm = tsolmLaiend(Node, 2) + tsolmLaiend(Node, 3) + tsolmLaiend(Node, 4);
                Kas_on = Kas_on_toe_solm;
                if (Kas_on_toe_solm != 0)
                    if (uvfiLaiend(Node) - 3) == 0
                        ontugesid = 3;
                    endif
                    TSMLD1 = tsolmLaiend(Node, 2) + tsolmLaiend(Node, 3);
                    if (TSMLD1 - 2) == 0
                        ontugesid = 2;
                    endif
                    if (TSMLD1 - 1) == 0
                        if ((tsolmLaiend(Node, 2)) - 1) == 0
                            ontugesid = 1;
                        endif
                        if ((tsolmLaiend(Node, 3)) - 1) == 0
                            ontugesid = 11;
                        endif
                    endif
                endif
                #
            endif
            #
        endfor %% for j=2:AVSdim(1,1)
     
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%r835 59    1    2    7    8    9    1
        if (ontugesid != 0)
            switch (ontugesid);
                case{1}
                Mitu = ontugesid
                SpTM1x0xz = SpTeisendusUhikMaatriks1x0v(1);
                SpTM1x0xz = - SpTM1x0xz;
                TAVS2 = toemuutuja;
                TVS2 = int2str(TAVS2);
                str2 = TVS2;
                str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM1x0xz', ')', Semi];
                cmd = sprintf(str)
                eval(cmd)
                Nr1 = Nr1 + 1;
                # taskVorreid(1, k) = 2;
                toemuutuja = toemuutuja + 1;
             
                #
                # %%%%%%%%%%%
                case{2}
                SpTM2x2xz = SpTeisendusUhikMaatriks2x2(1);
                SpTM2x2xz = - SpTM2x2xz;
                TAVS2 = toemuutuja;
                TVS2 = int2str(TAVS2);
                str2 = TVS2;
                str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM2x2xz', ')', Semi];
                cmd = sprintf(str)
                eval(cmd)
                Nr1 = Nr1 + 2;
                # taskVorreid(1, k) = 2;
                toemuutuja = toemuutuja + 2;
             
                #
                # %%%%%%%%%%%
                case{3}
                SpTM3x3xz = SpTeisendusUhikMaatriks(1);
                SpTM3x3xz = - SpTM3x3xz;
                TAVS2 = toemuutuja;
                TVS2 = int2str(TAVS2);
                str2 = TVS2;
                str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM3x3xz', ')', Semi];
                cmd = sprintf(str)
                eval(cmd)
             
                Nr1 = Nr1 + 3;
                toemuutuja = toemuutuja + 3;
                #
                #
                # %%%%%%%%%%%
                case{4}
                #
                # %%%%%%%%%%%
                case{5}
                #
                # %%%%%%%%%%%
                case{6}
                #
                # %%%%%%%%%%%
                case{7}
                #
                # %%%%%%%%%%%
                case{8}
                #
                # %%%%%%%%%%%
                case{9}
                #
                # %%%%%%%%%%%
                case{10}
                #
                # %%%%%%%%%%%
                case{11}
                SpTM0x1xz = SpTeisendusUhikMaatriks0x1v(1);
                SpTM0x1xz = - SpTM0x1xz;
                TAVS2 = toemuutuja;
                TVS2 = int2str(TAVS2);
                str2 = TVS2;
                str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTM0x1xz', ')', Semi];
                cmd = sprintf(str)
                eval(cmd)
                Nr1 = Nr1 + 1;
                # taskVorreid(1, k) = 2;
                toemuutuja = toemuutuja + 1;
                # %%%%%
        endswitch
    endif
    #
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   strSolmB=num2str(Solm);
    %    disp('      Nodal forces at the node')
    %      strB=['B(',strB1, koolonV, strB2, Kom,Arv1,')','=','sSolmF(', Arv1, koolonV, Arv123, Kom, %strSolmB,')',Semi];
    %      cmd = sprintf(strB)
    %      eval(cmd)
    # %%%%%%%%%%%%%%%%
    sF(:, 1) = sSolmF(:, 1, Solm);
    strSolmB = Arv1;
    strB = ['B(', strB1, koolonV, strB2, Kom, Arv1, ')', '=', 'sF(', Arv1, koolonV, Arv123, Kom, Arv1, ')', Semi];
    cmd = sprintf(strB)
    eval(cmd)
 
    # %%%%%%%%%%%%%%%%
    maxtaskVorreid = max(taskVorreid(1, k));
    #
    vorrandeid = size(spA);
    spRidu = vorrandeid(1, 1);
    Nr1 = spRidu + 1;
 
endif # %if AVSdim(1,1) > 1
endfor # % for i=1:NSARV
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vorrandeid = size(spA);
spRidu = vorrandeid(1, 1);
disp('----- ')
vorrandeidE = size(spA);
disp(' ')
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
#
for i = 1:ABBdim(1, 1)
    #
    if ABB(i, 10) == 1
        str1 = num2str(Nr1);
        Arv1 = num2str(1);
        TAVS2 = eval(sprintf('%3i', ABB(i, 6)));
        TVS2 = int2str(TAVS2);
        str2 = TVS2;
        str = ['spA=spSisestaArv(spA', Kom, str1, Kom, str2, Kom, Arv1, ')', Semi];
        #
        cmd = sprintf(str)
        eval(cmd);
        Nr1 = Nr1 + 1;
    endif
    #
    if ABB(i, 11) == 1
        str1 = num2str(Nr1);
        Arv1 = num2str(1);
        TAVS2 = eval(sprintf('%3i', ABB(i, 7)));
        TVS2 = int2str(TAVS2);
        str2 = TVS2;
        str = ['spA=spSisestaArv(spA', Kom, str1, Kom, str2, Kom, Arv1, ')', Semi];
        #
        cmd = sprintf(str)
        eval(cmd);
        Nr1 = Nr1 + 1;
    endif
endfor
Nr1 = Nr1 - 1;
#
disp('----- ')
vorrandeidN = size(spA);
#
disp(' ')
spA_rows = vorrandeidN(1, 1)
spA_rowsN = spA_rows;
spA_cols = vorrandeidN(1, 2)
#
nnzN = nnz(spA);
spA_nnz = sprintf('non_zero_elements_in_spA: %d', nnz(spA))
disp('  ')
disp('----- Side conditions are inserted into spA  ----  ')
disp('  ')
side_condition_rows = spA_rowsN - spA_rowsE
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
tsolm = tsolm;
TSARV = size(tsolm);
NTSARV = TSARV(1, 1);
toereaktsioonid_algavad = Nr1;
#
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Restrtictions on support displacements ')
From_rows = sprintf('Restrtictions on support displacements begin from row: %d', Nr1)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
#
toesiirdedN = 1;
TSARV = size(tsolm);
NTSARV = TSARV(1, 1);
ABBMoot = ABBdim(1, 1);
ToeSidemeteArv = 0;
ToeSiireRealNR = toereaktsioonid_algavad - 1;

for i = 1:NTSARV
    tsolNr = tsolm(i, 1);
    for j = 2:4
        if (tsolm(i, j) == 1)
            ToeSidemeteArv = ToeSidemeteArv + 1;
        endif
    endfor
endfor
ToeSidemeteArv;
disp('   ')
ToeSidemed = zeros(ToeSidemeteArv, 6);
ToeSolmes = zeros(NTSARV, 8);
ToeSolmes(:, 1:4) = tsolm(:, 1:4);

for i = 1:NTSARV
    tsolNr = tsolm(i, 1);
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
for j = 1:NTSARV
    for k = 2:4
        KusOnYks = tsolm(j, k);
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
for i = 1:ToeSidemeteArv
    Nr11 = ToeSidemed(i, 1);
    str1 = num2str(Nr11);
    Arv1 = num2str(1);
    UWFi = ToeSidemed(i, 3);
    #
    switch (UWFi);
        case{1}
        VarrasS = ToeSidemed(i, 2);
        SpTUv = SpToeSiirdeUvektor(NSARV, NEARV, VarrasS, krdn, selem);
        TAVS = ToeSidemed(i, 4);
        TVS2 = int2str(TAVS);
        str2 = TVS2;
     
        str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTUv', ')', Semi];
        cmd = sprintf(str)
        eval(cmd)
        #
        strB1 = str1;
        SolmB = ToeSidemed(i, 7);
        strSolmB = num2str(SolmB);
        ToeSidemed(1, 1);
        Arv1x = num2str(1);
        shift1 = tSiire(1, 1, SolmB);
     
        strShift1 = num2str(shift1);
        %    disp('  ')
        strTxt1 = ['Support shift  ' 'tSiire(' Arv1x Kom Arv1 Kom strSolmB ')' ' = ' strShift1 '  in x direction at the node ' strSolmB];
        disp(strTxt1)
        strB = ['B(', strB1, Kom, Arv1, ')', '=', 'tSiire(', Arv1x, Kom, Arv1, Kom, strSolmB, ')', Semi];
        cmd = sprintf(strB)
        eval(cmd)
        disp('  ')
        case{2}
        VarrasS = ToeSidemed(i, 2);
        SpTWv = SpToeSiirdeWvektor(NSARV, NEARV, VarrasS, krdn, selem);
        %TAVS=ToeSidemed(i,5);
        TAVS = ToeSidemed(i, 4);
        TVS2 = int2str(TAVS);
        str2 = TVS2;
     
        str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTWv', ')', Semi];
        cmd = sprintf(str)
        eval(cmd)
        #
        strB1 = str1;
        SolmB = ToeSidemed(i, 7);
        strSolmB = num2str(SolmB);
        ToeSidemed(1, 1);
        Arv1z = num2str(2);
        shift2 = tSiire(2, 1, SolmB);
        strShift2 = num2str(shift2);
        %    disp('  ')
        strTxt2 = ['Support shift  ' 'tSiire(' Arv1z Kom Arv1 Kom strSolmB ')' ' = ' strShift2 '  in z direction at the node ' strSolmB];
        disp(strTxt2)
        strB = ['B(', strB1, Kom, Arv1, ')', '=', 'tSiire(', Arv1z, Kom, Arv1, Kom, strSolmB, ')', Semi];
        cmd = sprintf(strB)
        eval(cmd)
        disp('  ')
        case{3}
        VarrasS = ToeSidemed(i, 2);
        SpTFiV = SpToeSiirdeFiVektor(1);
        %TAVS=ToeSidemed(i,6);
        TAVS = ToeSidemed(i, 4);
        TVS2 = int2str(TAVS);
        str2 = TVS2;
     
        str = ['spA=spInsertBtoA(spA', Kom, str1, Kom, str2, Kom, 'SpTFiV', ')', Semi];
        cmd = sprintf(str)
        eval(cmd)
        #
        strB1 = str1;
        SolmB = ToeSidemed(i, 7);
        strSolmB = num2str(SolmB);
        ToeSidemed(1, 1);
        Arv1y = num2str(3);
        shift3 = tSiire(3, 1, SolmB);
        strShift3 = num2str(shift3);
        %    disp('  ')
        strTxt3 = ['Support shift  ' 'tSiire(' Arv1y Kom Arv1 Kom strSolmB ')' ' = ' strShift3 '  in y direction at the node ' strSolmB];
        disp(strTxt3)
        strB = ['B(', strB1, Kom, Arv1, ')', '=', 'tSiire(', Arv1y, Kom, Arv1, Kom, strSolmB, ')', Semi];
        cmd = sprintf(strB)
        eval(cmd)
        disp('  ')
        # %%%%%
endswitch
#
endfor # for ToeSidemeteArv

# %%%%%%%%%%%%%%%
#
vorrandeid = size(spA);
spRidu = vorrandeid(1, 1);
#
disp('----- ')
vorrandeidR = size(spA);
#
disp(' ')
spA_rows = vorrandeidR(1, 1)
spA_rowsR = spA_rows;
spA_cols = vorrandeidR(1, 2)
#
nnzR = nnz(spA);
spA_nnz = sprintf('non_zero_elements_in_spA: %d', nnz(spA))
disp('  ')
disp('----- Restrtictions equations are inserted into spA  ----  ')
disp('  ')
restrtictions_equations_rows = spA_rowsR - spA_rowsN
non_zero_elements_in_restrtictions_equations = nnzR - nnzN
#
disp('  ')
#
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
#
spA_rank = sprank(spA)
spA
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
for i = 1:NEARV
    elemendiN(i, 1) = i;
    siireVardaA(i, 1:3) = selemjl(i, 7:9);
    siireVardaL(i, 1:3) = selemjl(i, 1:3);
    joudVardaA(i, 1:3) = selemjl(i, 10:12);
    joudVardaL(i, 1:3) = selemjl(i, 4:6);
endfor
#
siireVardaA;
siireVardaL;
joudVardaA;
joudVardaL;
# %%%%%%%%%%%%
for i = 1:NEARV
    JrN(i, 1) = i;
endfor
#
JrNT = JrN';
#
disp('=======================================================')
disp('     Support displacements/(Cx, Cz, CMy) at nodes  ')
disp('   X_No,  Element_No,  u/w/fi_No, u_No,  w_No,  fi_No ')
disp('-------------------------------------------------------')
for i = 1:ToeSidemeteArv
    disp(sprintf('     %2i     %2i         %2i        %3i    %3i     %3i  ', ToeSidemed(i, 1), ToeSidemed(i, 2), ToeSidemed(i, 3), ToeSidemed(i, 4), ToeSidemed(i, 5), ToeSidemed(i, 6)))
endfor
disp('-------------------------------------------------------')
#
#The_displacement_DaF - s_numbers_at_beginning_of_the_element = sprintf('elementNr, displacement: %d', JrN(i, 1) siireVardaA)
disp('==================================================================')
disp(' Displacements and force numbers at beginning of the element  ')
disp('   No,    u,   w,   fi     N,   Q,   M  ')
disp('------------------------------------------------------------------')
for i = 1:NEARV
    disp(sprintf('  %2i     %2i   %2i   %2i     %2i   %2i   %2i  ', JrN(i, 1), siireVardaA(i, 1), siireVardaA(i, 2), siireVardaA(i, 3), joudVardaA(i, 1), joudVardaA(i, 2), joudVardaA(i, 3)))
endfor
%
disp('-------------------------------------------------------')
disp(' ')
disp('==================================================================')
disp(' Displacements and forces numbers at end of the element  ')
disp('   No,    u,   w,   fi     N,   Q,   M  ')
disp('------------------------------------------------------------------')
for i = 1:NEARV
    disp(sprintf('  %2i     %2i   %2i   %2i     %2i   %2i   %2i ', JrN(i, 1), siireVardaL(i, 1), siireVardaL(i, 2), siireVardaL(i, 3), joudVardaL(i, 1), joudVardaL(i, 2), joudVardaL(i, 3)))
endfor
#
disp('-------------------------------------------------------')
disp(' ')
#
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' Initial parmeter vectors for elements displacements and forces  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
disp(' ')
disp('-------- Descaling multiplier for the displacements = 1/baasi0 --------  ')
disp(' ')
#
for i = 1:NEARV
    # %siireVardaA(i,1:3)
    AlgPar(i, 1:6) = [X(siireVardaA(i, 1), 1) / baasi0 X(siireVardaA(i, 2), 1) / baasi0 X(siireVardaA(i, 3), 1) / baasi0 ...
    X(joudVardaA(i, 1), 1) X(joudVardaA(i, 2), 1) X(joudVardaA(i, 3), 1)];
endfor
#
disp('============================================================================')
disp(' Unscaled initial parmeter vector ')
disp('Element No    u          w          fi             N          Q          M ')
disp('----------------------------------------------------------------------------')
for i = 1:NEARV
    disp(sprintf('  %2i   %9.3e   %9.3e   %9.3e    %9.3f  %9.3f  %9.3f', JrN(i, 1), AlgPar(i, 1), AlgPar(i, 2), AlgPar(i, 3), AlgPar(i, 4), AlgPar(i, 5), AlgPar(i, 6)))
    #AlgPar
endfor
disp('----------------------------------------------------------------------------')
#
JJ = 0;
for i = 1:NEARV
    #siireVardaL(i, 1:3)
    JJ = JJ + 1;
    SiireJoud(JJ, 1:6) = [X(siireVardaA(i, 1), 1) / baasi0 X(siireVardaA(i, 2), 1) / baasi0 X(siireVardaA(i, 3), 1) / baasi0 ...
      X(joudVardaA(i, 1), 1) X(joudVardaA(i, 2), 1) X(joudVardaA(i, 3), 1)];
    JJ = JJ + 1;
    SiireJoud(JJ, 1:6) = [X(siireVardaL(i, 1), 1) / baasi0 X(siireVardaL(i, 2), 1) / baasi0 X(siireVardaL(i, 3), 1) / baasi0 ...
      X(joudVardaL(i, 1), 1) X(joudVardaL(i, 2), 1) X(joudVardaL(i, 3), 1)];
endfor
#
NEARV2 = NEARV * 2;
Nn = 0;
Nn2 = 0;
for i = 1:NEARV
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
for i = 1:NEARV2
    disp(sprintf('  %2i   %9.3e   %9.3e   %9.3e    %9.3f  %9.3f  %9.3f', JrN2(i, 1), SiireJoud(i, 1), SiireJoud(i, 2), SiireJoud(i, 3), SiireJoud(i, 4), SiireJoud(i, 5), SiireJoud(i, 6)))
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

