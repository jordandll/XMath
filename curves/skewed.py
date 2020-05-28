""" Skewed Parabola Module
A skewed parabola is given by a function of the form:
	f(x | a, b, c, m, s) := ax^2 + b*sqrt(mx + s) + c, where a != 0;  b != 0; m in {-1,1}; and a*b < 0."""

from math import sqrt

def make_f (a, b, c, m, s):
	""" Makes a function that gives a skewed parabola using the coefficients passed to this function.
	The generalized and parametric form of said functon is defined in the documentation of this module."""
	check_params(a, b, m)
	return lambda x: a*x**2 + b*sqrt(m*x + s) + c
	
	
def dir_x (m):
	""" Direction of the curve along the x-axis.  Can be either positive or negative. """
	return m
	
def dir_y (a):
	""" Direction of the curve along the y-axis.  Can be either positive or negative. """
	return -a
	
def starting_point (m, s):
	""" The x-coordinate of the starting point of the curve on the x-axis. """
	return -m*s
	
def extremum (a, b, c, m, s):
	""" The x, y coordinates of the extremum.  This is the point where the slope of the curve is zero."""
	check_params(a, b, m)
	x = (m*b**2)/(4*a**2) - s/m
	f = make_f(a, b, c, m, s)
	return(x, f(x))
	
def roots (a, b, c, m, s) -> set:
	""" Returns the roots of the function that gives the skewed parabola, denoted as ‘f’.  These are the set of values that satisfy the equation, f(r)=0, and are the y-intercept(s) of the curve. """
	check_params(a, b, m)
	q = (2*a*c - m*b**2) / (2*a**2)
	d = q**2 + (s-c**2)/(a**2)
	delta = sqrt(d)
	return {delta - q, -delta - q}
	
def check_params (a, b, m):
	""" Check the the parameters against the conditions found in the clause of the generalized function defition that gives a skewed parabola.
	Said definition can be found in the documentation of this module. """
	if a == 0: raise ValueError()
	elif b == 0: raise ValueError()
	elif m not in {-1, 1}: raise ValueError()
	elif a*b > 0: raise ValueError()
	
	
