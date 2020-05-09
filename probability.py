#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  probability.py
#  
#  Copyright 2020 Jordan <jordan@HPEnvy>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

from math import log

def probabilities (X) -> dict:
	""" This function maps the set of outcomes found in the sequence of events, 'X', to their respective probabilty of occuring in 'X'.
	The return value is a python dictionary where the keys are the set of outcomes and the values are their associated probabilities."""
	# The set of outcomes, denoted as 'C', and the total events, denoted as 'T'.
	C, T = set(X), len(X)
	return {c: X.count(c) / T for c in C}
	
def shannon (P, b : int = None) -> float:
	""" Shannon or Informational Entropy 
	Calculates the shannon entropy, denoted as 'H', of a discrete random variable, 'X'.
	Parameters:
	1.) P:  An iterable of probabilities associated with the set of outcomes of 'X'.
	2.) b:  The logarithmic base used in calculating 'H'.  The default value is len(P). """
	# Check params.
	if b == None: b = len(P)
	if b <= 0: raise ValueError('Argument Value Error:\tThe logarithmic base must be a positive integer.')
	# Init return value.
	res = 0
	for p in P:
		res += p*log(p)/log(b)
	return -res
	
	

def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
