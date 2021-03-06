##=========================================================================
## This program is written by Andres Lahe <andres.lahe@ttu.ee>, 2010-07-16
## Last modified: Andres Lahe, 2010-08-03
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

function SF = SisejoudPunktism(VardaNr, X, AlgPar, lvarras, selem, esQkoormus, esFjoud, suurused)
## The forces of frame element, ElementNr, at x = X

if nargin != 8
    error(' function SisejoudPunktis has wrong number of input arguments!')
endif

i = VardaNr;
xx = X;
baasi0 = 1.0;

EI = selem(i, 13); %  from topology
EA = selem(i, 14);
GAr = selem(i, 15);
Li = lvarras(i, 1);

Fjoud = esFjoud(:, 1:3, i);
qkoormus = esQkoormus(:, 1:4, i);

AP = AlgPar(i, :)';
# --------- The transfer matrix equation --------
vvF = ylfhlin(1.0, xx, EA, GAr, EI);
Sj = ESTFrKrmus(baasi0, xx, Li, Fjoud, qkoormus, EA, EI);
vB = Sj;
Fvv(:, 1) = vvF * AP + vB;

disp(sprintf('%15s %2i %9s  %4.7f   ', 'Forces of element', VardaNr, ' at x =', X))

for i = 1:3
    disp(sprintf('%14s %9.5e   ', suurused(i, :), Fvv(i, 1)))
endfor

for i = 4:6
    disp(sprintf('%14s %9.5f   ', suurused(i, :), Fvv(i, 1)))
endfor
SF = Fvv(:, 1);
endfunction

