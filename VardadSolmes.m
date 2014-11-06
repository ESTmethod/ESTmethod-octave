# VardadSolmes.m 
function Sv=VardadSolmes(NSARV,NEARV,Solm,AB,ABB)
#
# Sort the rows of matrix ABB according to the order of the (column=8) node 
# number. Returned matrix with the order of the row reversed.
# 
##=========================================================================
## This Program is writtwn by Andres Lahe,   2013-06-26
##                    e-mail: andres.lahe@ttu.ee
## LAST MODIFIED: Andres Lahe,   2013-06-27
## Copyright (c)  2013 by Tallinn University of Technology 
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
#%%%%%%%%%%%%%%%%%
#if nargin != 5
if ~(nargin==5)
error(' function VardadSolmes have wrong number of input arguments!') 
end
#%%%%%%%%%%%%%%%%%
##
## NSARV - the number of frame_nodes
## NEARV - the number of elements
## Solm - the node number
## AB -  the elements connected at the nodes
## ABB - matrix[Element number, u w fi N Q M, Node, hinges N    Q    M ]  
##
#
NEARV2=NEARV*2;
NSARV;
k=1;
for i=1:NSARV
   for j=1:NEARV2
       ks=AB(j,1);
      if ks == i
        ABBU(k,1:11,i)=ABB(j,1:11);
        k=k+1;
      endif
   endfor
endfor

ABBU;
jl=1;
for i=1:NSARV
   ABBUXX=ABBU(:,:,i);
   [~,idx] = sortrows(ABBUXX(:,8));
   ABBYY=ABBUXX(idx,:);
   ABBYYY(:,:,jl)=flipud(ABBYY);
   jl=jl+1;
endfor
ABBYYY;
# 

#%%http://www.network-theory.co.uk/docs/octave3/octave_87.html
#if (AF1-eps1) > 0
switch (Solm);
case{1}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
ABBYYY1(:,:) = ABBYYY(tmp~=s(2),:,i);
ABBYYY1;
#ascending
[~,idx] = sortrows(ABBYYY1(:,11));
#%
ABBYYY1=ABBYYY1(idx,:);
ABBYYY0=ABBYYY1;

case{2}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY2(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY2(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY2;
#ascending
[~,idx] = sortrows(ABBYYY2(:,11));
#%
ABBYYY2=ABBYYY2(idx,:);
ABBYYY0=ABBYYY2;

case{3}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY3(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY3(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY3;
#ascending
[~,idx] = sortrows(ABBYYY3(:,11));
#%
ABBYYY3=ABBYYY3(idx,:);
ABBYYY0=ABBYYY3;

case{4}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY4(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY4(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY4; 
#ascending
[~,idx] = sortrows(ABBYYY4(:,11));
#%
ABBYYY4=ABBYYY4(idx,:);
ABBYYY0=ABBYYY4;


case{5}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY5(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY5(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY5;
#ascending
[~,idx] = sortrows(ABBYYY5(:,11));
#%
ABBYYY5=ABBYYY5(idx,:);
ABBYYY0=ABBYYY5;

case{6}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY6(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY6(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY6;
#ascending
[~,idx] = sortrows(ABBYYY6(:,11));
#%
ABBYYY6=ABBYYY6(idx,:);
ABBYYY0=ABBYYY6;


case{7}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY7(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY7(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY7;
#ascending
[~,idx] = sortrows(ABBYYY7(:,11));
#%
ABBYYY7=ABBYYY7(idx,:);
ABBYYY0=ABBYYY7;


case{8}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY8(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY8(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY8;
#ascending
[~,idx] = sortrows(ABBYYY8(:,11));
#%
ABBYYY8=ABBYYY8(idx,:);
ABBYYY0=ABBYYY8;


case{9}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY9(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY9(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY9;
#ascending
[~,idx] = sortrows(ABBYYY9(:,11));
#%
ABBYYY9=ABBYYY9(idx,:);
ABBYYY0=ABBYYY9;



case{10}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY10(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY10(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY10;
#ascending
[~,idx] = sortrows(ABBYYY10(:,11));
#%
ABBYYY10=ABBYYY10(idx,:);
ABBYYY0=ABBYYY10;


case{11}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY11(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY11(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY11;
#ascending
[~,idx] = sortrows(ABBYYY11(:,11));
#%
ABBYYY11=ABBYYY11(idx,:);
ABBYYY0=ABBYYY11;


case{12}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY12(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY12(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY12;
#ascending
[~,idx] = sortrows(ABBYYY12(:,11));
#%
ABBYYY12=ABBYYY12(idx,:);
ABBYYY0=ABBYYY12;


case{13}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY13(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY13(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY13;
#ascending
[~,idx] = sortrows(ABBYYY13(:,11));
#%
ABBYYY13=ABBYYY13(idx,:);
ABBYYY0=ABBYYY13;



case{14}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY14(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY14(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY14;
#ascending
[~,idx] = sortrows(ABBYYY14(:,11));
#%
ABBYYY14=ABBYYY14(idx,:);
ABBYYY0=ABBYYY14;


case{15}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY15(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY15(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY15;
#ascending
[~,idx] = sortrows(ABBYYY15(:,11));
#%
ABBYYY15=ABBYYY15(idx,:);
ABBYYY0=ABBYYY15;


case{16}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY16(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY16(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY16;
#ascending
[~,idx] = sortrows(ABBYYY16(:,11));
#%
ABBYYY16=ABBYYY16(idx,:);
ABBYYY0=ABBYYY16;


case{17}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY17(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY17(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY17;
#ascending
[~,idx] = sortrows(ABBYYY17(:,11));
#%
ABBYYY17=ABBYYY17(idx,:);
ABBYYY0=ABBYYY17;


case{18}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY18(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY18(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY18;
#ascending
[~,idx] = sortrows(ABBYYY18(:,11));
#%
ABBYYY18=ABBYYY18(idx,:);
ABBYYY0=ABBYYY18;


case{19}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY19(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY19(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY19;
#ascending
[~,idx] = sortrows(ABBYYY19(:,11));
#%
ABBYYY19=ABBYYY19(idx,:);
ABBYYY0=ABBYYY19;


case{20}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY20(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY20(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY20;
#ascending
[~,idx] = sortrows(ABBYYY20(:,11));
#%
ABBYYY20=ABBYYY20(idx,:);
ABBYYY0=ABBYYY20;  


case{21}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
%ABBYYY21(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY21(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY21;
#ascending
[~,idx] = sortrows(ABBYYY21(:,11));
#%
ABBYYY21=ABBYYY21(idx,:);
ABBYYY0=ABBYYY21;


case{22}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY22(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY22(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY22;
#ascending
[~,idx] = sortrows(ABBYYY22(:,11));
##
ABBYYY22=ABBYYY22(idx,:);
ABBYYY0=ABBYYY22;


case{23}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY23(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY23(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY23;
#ascending
[~,idx] = sortrows(ABBYYY23(:,11));
##
ABBYYY23=ABBYYY23(idx,:);
ABBYYY0=ABBYYY23;


case{24}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY24(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY24(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY24;
#ascending
[~,idx] = sortrows(ABBYYY24(:,11));
##
ABBYYY24=ABBYYY24(idx,:);
ABBYYY0=ABBYYY24;


case{25}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY25(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY25(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY25;
#ascending
[~,idx] = sortrows(ABBYYY25(:,11));
##
ABBYYY25=ABBYYY25(idx,:);
ABBYYY0=ABBYYY25;


case{26}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY26(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY26(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY26;
#ascending
[~,idx] = sortrows(ABBYYY26(:,11));
##
ABBYYY26=ABBYYY26(idx,:);
ABBYYY0=ABBYYY26;


case{27}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY27(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY27(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY27;
#ascending
[~,idx] = sortrows(ABBYYY27(:,11));
##
ABBYYY27=ABBYYY27(idx,:);
ABBYYY0=ABBYYY27;


case{28}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY28(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY28(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY28;
#ascending
[~,idx] = sortrows(ABBYYY28(:,11));
##
ABBYYY28=ABBYYY28(idx,:);
ABBYYY0=ABBYYY28;


case{29}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY29(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY29(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY29;
#ascending
[~,idx] = sortrows(ABBYYY29(:,11));
##
ABBYYY29=ABBYYY29(idx,:);
ABBYYY0=ABBYYY29;


case{30}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY30(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY30(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY30;
#ascending
[~,idx] = sortrows(ABBYYY30(:,11));
##
ABBYYY30=ABBYYY30(idx,:);
ABBYYY0=ABBYYY30;  


case{31}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY31(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY31(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY31;
#ascending
[~,idx] = sortrows(ABBYYY31(:,11));
##
ABBYYY31=ABBYYY31(idx,:);
ABBYYY0=ABBYYY31;


case{32}    
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY32(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY32(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY32;
#ascending
[~,idx] = sortrows(ABBYYY32(:,11));
##
ABBYYY32=ABBYYY32(idx,:);
ABBYYY0=ABBYYY32;


case{33}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY33(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY33(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY33;
#ascending
[~,idx] = sortrows(ABBYYY33(:,11));
##
ABBYYY33=ABBYYY33(idx,:);
ABBYYY0=ABBYYY33;


case{34}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY34(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY34(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY34;
#ascending
[~,idx] = sortrows(ABBYYY34(:,11));
##
ABBYYY34=ABBYYY34(idx,:);
ABBYYY0=ABBYYY34;


case{35}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY35(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY35(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY35;
#ascending
[~,idx] = sortrows(ABBYYY35(:,11));
##
ABBYYY35=ABBYYY35(idx,:);
ABBYYY0=ABBYYY35;


case{36}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY36(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY36(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY36;
#ascending
[~,idx] = sortrows(ABBYYY36(:,11));
##
ABBYYY36=ABBYYY36(idx,:);
ABBYYY0=ABBYYY36;


case{37}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY37(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY37(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY37;
#ascending
[~,idx] = sortrows(ABBYYY37(:,11));
##
ABBYYY37=ABBYYY37(idx,:);
ABBYYY0=ABBYYY37;


case{38}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY38(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY38(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY38;
#ascending
[~,idx] = sortrows(ABBYYY38(:,11));
##
ABBYYY38=ABBYYY38(idx,:);
ABBYYY0=ABBYYY38;


case{39}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY39(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY39(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY39;
#ascending
[~,idx] = sortrows(ABBYYY39(:,11));
##
ABBYYY39=ABBYYY39(idx,:);
ABBYYY0=ABBYYY39;


case{40}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY40(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY40(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY40;
#ascending
[~,idx] = sortrows(ABBYYY40(:,11));
##
ABBYYY40=ABBYYY40(idx,:);
ABBYYY0=ABBYYY40;


case{41}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY41(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY41(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY41;
#ascending
[~,idx] = sortrows(ABBYYY41(:,11));
##
ABBYYY41=ABBYYY41(idx,:);
ABBYYY0=ABBYYY41;


case{42}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY42(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY42(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY42;
#ascending
[~,idx] = sortrows(ABBYYY42(:,11));
##
ABBYYY42=ABBYYY42(idx,:);
ABBYYY0=ABBYYY42;


case{43}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY43(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY43(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY43;
#ascending
[~,idx] = sortrows(ABBYYY43(:,11));
##
ABBYYY43=ABBYYY43(idx,:);
ABBYYY0=ABBYYY43;


case{44}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY44(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY44(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY44;
#ascending
[~,idx] = sortrows(ABBYYY44(:,11));
##
ABBYYY44=ABBYYY44(idx,:);
ABBYYY0=ABBYYY44;


case{45}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY45(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY45(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY45;
#ascending
[~,idx] = sortrows(ABBYYY45(:,11));
##
ABBYYY45=ABBYYY45(idx,:);
ABBYYY0=ABBYYY45;


case{46}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY46(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY46(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY46;
#ascending
[~,idx] = sortrows(ABBYYY46(:,11));
##
ABBYYY46=ABBYYY46(idx,:);
ABBYYY0=ABBYYY46;


case{47}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY47(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY47(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY47;
#ascending
[~,idx] = sortrows(ABBYYY47(:,11));
##
ABBYYY47=ABBYYY47(idx,:);
ABBYYY0=ABBYYY47;


case{48}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY48(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY48(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY48;
#ascending
[~,idx] = sortrows(ABBYYY48(:,11));
##
ABBYYY48=ABBYYY48(idx,:);
ABBYYY0=ABBYYY48;


case{49}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY49(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY49(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY49;
#ascending
[~,idx] = sortrows(ABBYYY49(:,11));
##
ABBYYY49=ABBYYY49(idx,:);
ABBYYY0=ABBYYY49;


case{50}
i=Solm;
#
s = size(ABBYYY(:,:,i));
inds = (ABBYYY(:,:,i)==0);
tmp = sum(inds');
#ABBYYY50(:,:,i) = ABBYYY(tmp~=s(2),:,i);
ABBYYY50(:,:) = ABBYYY(tmp~=s(2),:,i);
#
ABBYYY50;
#ascending
[~,idx] = sortrows(ABBYYY50(:,11));
##
ABBYYY50=ABBYYY50(idx,:);
ABBYYY0=ABBYYY50;


case{51}
disp(' Pole veel 51-est sõlme ')

case{52}
case{53}
case{54}
case{55}
case{56}
case{57}
case{58}
case{59}
case{60}


#%?otherwise
endswitch
#else
# disp(' Sõlm puudub  '); 
%aaaaaaaaaa=0;
#endif

Sv=ABBYYY0;
#
endfunction