## InsertBtoA.m
function A=InsertBtoA(A,I,J,IM,JN,B,M,N)
#
# Inserted matrix B (dimension M, N)  into matrix A (dimension I, J)
# starting at row index IM and column index JN.
#
##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2002-02-15
## Last modified: Andres Lahe, 2009-01-03
## Copyright (c)  2002 by Tallinn University of Technology
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
#%%%%%%%%%%%%%%%%%
#if nargin != 8
if ~(nargin==8)
error(' function InsertBtoA have wrong number of input arguments!')
end
#%%%%%%%%%%%%%%%%%
#
IMB=IM+M-1;
JNB=JN+N-1;
n1=0;
m1=0;
for i=IM:IMB
     n1=0;
     m1=m1+1;
  for j=JN:JNB
        n1=n1+1;
        CC=B(m1,n1);
       A(i,j)=A(i,j)+CC;
  endfor
endfor
A;
#
endfunction
