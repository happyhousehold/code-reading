FPC: Fixed Point Continuation
Copyright (C) 2007-2008 Elaine Hale, Wotao Yin and Yin Zhang
-------------------------------------------------------------

1.  Introduction
================

Fixed-point continuation (FPC) is an algorithm for solving large-scale l-1 
regularized least-squares problems

min_x ||x||_1 + (mu/2)*||Ax - b||_M^2

using fixed-point iterations paired with a continuation strategy for mu.  

2.  Getting Started
===================

The files in the drivers folder generate example compressed sensing 
problems and solve them with fpc.m or fpc_bb.m.  To run them, set your 
current directory to .../fpc/drivers, and call bb_run or basic_run from the 
command line.  The files automatically add the appropriate folders to your 
Matlab path (temporarily).  Note that some of the example problem types 
(determined by Ameth) require the signal processing toolbox.

3.  Optional External Packages
==============================

FPC is compatible with the SPARCO toolbox, which contains a number of 
problems for benchmarking sparse reconstruction algorithms, along with a
library of linear operators that can be used to generate additional 
problems.  

    SPARCO          http://www.cs.ubc.ca/labs/scl/sparco/

4.  Contact Information
=======================

FPC is available at http://www.caam.rice.edu/~optimization/L1/.  Please
feel free to e-mail the authors with any comments or suggestions:

    Elaine Hale     <ehale@rice.edu>
    Wotao Yin       <wotao.yin@rice.edu>
    Yin Zhang       <zhang@caam.rice.edu>

5.  Copyright Notice
====================

FPC is free software; you can redistribute it and/or modify it under the 
terms of the GNU General Public License as published by the Free Software 
Foundation; either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.