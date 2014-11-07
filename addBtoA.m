##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-07-16
## Last modified: Mattias Põldaru, 2014-11-06
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
%> @file addBtoA.m
%> @brief Function to add one matrix to another at a specific position.
%======================================================================
%> @brief addBtoA adds matrix B to matrix A by element.
%>
%> Add matrix B elements to matrix A elements starting at row i and column j.
%> Please note that matrix A has to be large enough for matrix B at given
%> position or one will get "row/element index out of bounds" error.
%>
%> @param A The original matrix.
%> @param B The matrix to be inserted at row i and column j.
%> @param i Starting row of to be inserted matrix B or position matrix [i j].
%> @param j Starting column of to be inserted matrix B (optional).
%>
%> @retval A The new combined matrix.
%======================================================================
function A = addBtoA(A, B, i, j)
if nargin == 4
elseif nargin == 3
    j = i(2);
    i = i(1);
else
    error('Function addBtoA needs 3 or 4 input arguments!')
end

endpos = size(B) + [i j] - 1;
m = 0;
for k = i : endpos(1)
    m += 1;
    n = 0;
    for l = j : endpos(2)
        n += 1;
        A(k, l) = A(k, l) + B(m, n);
    endfor
endfor
endfunction

