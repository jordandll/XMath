# polynomial.py

from math import sqrt
from collections.abc import Sequence
from numbers import Number

def deg2 (a, b, c):
	""" Generate a Degree 2 Polynomial 
	Returns a functon, denoted as 'f(x | a, b, c)=ax^2+bx+c', where 'a', 'b', and 'c' are equal to the arguments passed to their respective parameters in this function."""
	return lambda x: a*x**2 + b*x + c
	
def degN (*args):
	""" Generate a Degree N Polynomial
	Returns a functon, denoted as 'f(x | a_0, a_1, ... , a_i, ... , a_N)= a_N*x^N + ... +a_i*x^i + ... + a_1*x + a_0, where N=len(args)-1.
	The elements in 'args' equal the coefficients of their corresponding term in the function, 'f';  And the index of each element in 'args' is equal to the 
	exponent of the variable, 'x', in it's corresponding term in 'f'.
	
	Example:  An argument list of [5, 1, 2] will result in the function, f(x) = 2x^2 + x + 5, being returned.
	"""
	return lambda x: sum(a*x**i for i, a in enumerate(args))
	
def extremum (a, b, c) -> tuple:
	""" Returns the (x, y) coordinates of the extremum of the curve given by the polynomial, f(x | a, b, c).
	The extremum can refer to either a maximum or minimum value.  When 'a' is negative, the max or top of the curve is returned.  Otherwise, the min is returned.
	The value of the x-coordinate can be thought of as the midpoint of the curve."""
	# Check params.
	if a == 0:
		raise ValueError('Argument Value Error:\tThe parameter, \'a\', in a polynomial, f(x) = ax^2 + bx + c, cannot be equal to zero.')
	x = -b/(2*a)
	return (x, a*x**2 + b*x + c)
	
def roots (a, b, c) -> set:
	""" Find the Roots of a Quadratic Polynomial
	Returns the set of values that satisfies the equation, ax^2 + bx + c = 0. """
	# Check for argument errors.
	if a == 0: raise ValueError('Argument Value Error:\tThe parameter, \'a\', in a polynomial, f(x) = ax^2 + bx + c, cannot be equal to zero.')
	t = b**2 - 4*a*c
	if b == 0 and c == 0: res = {0}
	elif t < 0: res = set()
	else:
		mp = -b/(2*a)
		delta = sqrt(t)/(2*a)
		res = {mp+delta, mp-delta}
	return res
	
def add (p1, p2):
	""" Adds the two Polynomials, 'p1' and 'p2'
	The arguments can be a sequence of coefficients or an instance of the Polynomial class. """
	res = [x[0] + x[1] for x in zip(p1, p2)]
	n = len(res)
	res.extend(max(p1, p2, key=len)[n:])
	return res
	
	
class Polynomial:
	""" The base class for all Polynomials """
	
	def __init__ (self, *args):
		__doc__ = degN.__doc__
		
		self.coefficients = tuple(args)
		self.f = degN(*args)
		
	def __call__ (self, x = None):
		""" This method behaves differently depending on the argument type of 'x'.
		1.)  When 'x' is a Number, it returns the f(x), where 'f' is the underyling function of this Polynomial.
		2.)  When 'x' is a Polynomial or a Sequence of coefficients it returns the resulting Polynomial of performing function composition on this Polynomial and 'x'.
		3.)  When 'x' is None, it returns this Polynomial, which is useful when performing function composition.  For example, 'P1(P2())', will perform function composition on 
		the Polynomials, P1 and P2, and return the resulting Polynomial."""
		pass
		if isinstance(x, Number): return self.f(x)
		elif isinstance(x, Polynomial):
			# Perform function composition and return a new Polynomial.
			# TODO:  Implement this part of the method.
			raise NotImplementedError()
		elif isinstance(x, Sequence):
			# Perform function composition and return a new Polynomial.
			# TODO:  Implement this part of the method.
			raise NotImplementedError()
		elif x == None: return self
		else: raise TypeError(f'Argument Type Error:\t{__name__} only takes instances of Numbers, Polynomials, and Sequences as arguments.')

	def degree(self) -> int: return len(self.coefficients)-1
	
	def __len__(self) -> int:
		return len(self.coefficients) 
	
	def __getitem__ (self, i: int):
		""" Get the coefficient at index 'i'. Note that 'i' equals the exponent of the variable, 'x', in the corresponding term of the polynomial."""
		return self.coefficients[i]
		
	def __add__ (self, rhs):
		__doc__ = add.__doc__
		if isinstance(rhs, Polynomial): res = Polynomial(*add(self.coefficients, rhs.coefficients))
		elif isinstance(rhs, Sequence): res = Polynomial(*add(self.coefficients, rhs))
		else: raise TypeError('Argument Type Error:\tOnly a Polynomial or a sequence of coefficients can be added to a Polynomial.')
		return res
		
	def __mul__ (self, rhs):
		"""  Multiply two Polynomials """
		if isinstance(rhs, Polynomial):
			# Add the degrees to find the degree, 'n', of the resulting polynomial.  Then initialize a list, denoted as 'res', of coefficients for said polynomial.
			n = self.degree() + rhs.degree()
			res = [0 for i in range(n+1)]
			# Distribute the right hand side, denoted as 'rhs', polynomial to the left hand side, 'lhs', polynomial -- or 'self' -- and then perform polynomial addition.
			for idx_lhs, val_lhs in enumerate(self.coefficients):
				for idx_rhs, val_rhs in enumerate(rhs.coefficients):
					res[idx_lhs + idx_rhs] += val_lhs * val_rhs
		else:
			raise TypeError('Argument Type Error:\tOnly a Polynomial or a sequence of coefficients can be added to a Polynomial.')
			
		return Polynomial(*res)
				
			
		
		
		
		
	
		
	
