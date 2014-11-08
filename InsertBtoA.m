##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2002-02-15
## Last modified: Andres Lahe, 2009-01-03
##                Mattias Põldaru 2014-11-07
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

%======================================================================
%> @file InsertBtoA.m
%> @brief Function to add one matrix to another at a given position
%> (deprecated, use addBtoA()).
%======================================================================
%> @brief InsertBtoA adds one matrix B into matrix A (deprecated).
%>
%> Add matrix B (dimension M, N) elements to matrix A (dimension I, J)
%> elements starting at row index IM and column index JN.
%> Deprecated in favour of addBtoA().
%>
%> @param A The original matrix.
%> @param I Height of matrix A (not used).
%> @param J Width of matrix A (not used).
%> @param IM Starting row of to be inserted matrix B.
%> @param JN Starting column of to be inserted matrix B.
%> @param B The matrix to be inserted at row IM and column JN.
%> @param M Height of matrix B (not used).
%> @param N Width of matrix B (not used).
%>
%> @retval The new combined matrix.
%======================================================================
function A = InsertBtoA(A, I, J, IM, JN, B, M, N)
if nargin != 8
    error('Function InsertBtoA has wrong number of input arguments!')
end
A = addBtoA(A, B, IM, JN)
warning('InsertBtoA() is deprecated, use addBtoA() instead.')
endfunction

