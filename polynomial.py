# polynomial.py

from math import sqrt, factorial
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

def sub (p1, p2):
	""" Subtracts the two Polynomials, 'p1' and 'p2'
	The arguments can be a sequence of coefficients or an instance of the Polynomial class. """
	res = [x[0] - x[1] for x in zip(p1, p2)]
	n = len(res)
	T = max(p1, p2, key=len)[n:]
	res.extend((-t for t in T) if len(p2) > len(p1) else T)
	return res
	
class Polynomial:
	""" The base class for all Polynomials """
	
	def __init__ (self, *args):
		__doc__ = degN.__doc__
		
		if args[-1] == 0: raise ValueError('Arugment Value Error:\tThe leading coefficient of a polynomial cannot be equal to zero.')

		self.coefficients = tuple(args)
		self.f = degN(*self.coefficients)
	
	@staticmethod
	def from_critical_points(*CP):
		""" Initialize a polynomial from a set of critical points, denoted as 'CP', instead of a collection of coefficients(the default method of initialization).
		Notes:
			1.)  Each argument should be a pair of (x, y) coordinates for each critical point.
			2.)  Said pair can be in the form of a builtin python 2-tuple or some other indexable or subscriptable collection such that for all arguments, 'a',
			in 'CP', a[0] returns the x-coordinate of the critical point and a[1] returns the y-coordinate.
			3.)  'CP' should be ordered under magnitude in ascending order with respect to the x-coordinates -- i.e. CP[0].x < CP[1].x < ... < CP[-1].x."""
		pass
		
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

	def __sub__ (self, rhs):
		__doc__ = add.__doc__
		if isinstance(rhs, Polynomial): res = Polynomial(*sub(self.coefficients, rhs.coefficients))
		elif isinstance(rhs, Sequence): res = Polynomial(*sub(self.coefficients, rhs))
		else: raise TypeError('Argument Type Error:\tOnly a Polynomial or a sequence of coefficients can be subtracted from a Polynomial.')
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

	def __pow__(self, n: int):
		""" Raise this polynomial to the n'th power, where 'n' is a non-negative whole number. """
		from itertools import product as cartesian_product
		# Check 'n'.
		if n < 0:
			raise ValueError('Argument Value Error:\tRaising a polynomial to a negative power is prohibited (atleast for now).')
		elif n > 1:
			deg1 = self.degree() * n
			C = [0 for i in range(deg1 + 1)]
			for T in cartesian_product(range(len(self.coefficients)), repeat=n):
				p = sum(T)
				c = 1
				for t in T: c *= self.coefficients[t]
				C[p] += c
			return Polynomial(*C)
		elif n == 0: return 1
		else: return Polynomial(*self.coefficients)

	def derivative (self, n: int = 1):
		""" Perform the derivative function on this polynomial 'n' times.
		Returns d^n(f)/(dx)^2, where 0 < n <= degree(f) and f(x) is this polynomial."""
		
		if n < 1: raise ValueError('Argument Value Error:\tThe power of the derivative function must be a positive integer.')
		elif n > self.degree(): raise ValueError('Argument Value Error:\tThe power of the derivative function cannot be greater than the degree of it\'s polynomial argument.')
		
		return Polynomial(*tuple(c * factorial(i) // factorial(i-n) for i, c in enumerate(self.coefficients[n:], n)))
		
	def anti_derivative(self, n: int =1):
		"""Perform the anti-derivative, sometimes known as the primitive integral, function on this polynomial 'n' times, where n > 0."""
		# Check parameter.
		if n < 1: raise ValueError('Argument Value Error:\tThe power of the derivative function must be a positive integer.')
		
		C = tuple(c * factorial(i) / factorial(i+n) for i, c in enumerate(self.coefficients))
		return tuple(0 for i in range(n)) + C

class Linear (Polynomial):
	""" Linear Polynomial
	A linear polynomial is of the form:
	f(x | a, b) := ax + b, where a != 0."""

	def __init__(self, a, b):
		""" Initialize a linear polynomial of the form:
		f(x | a, b) := ax + b, where a != 0. """
		
		if a == 0: raise ValueError('Argument Value Error:\tThe leading coefficent of a polynomial cannot be equal to zero.')
		
		Polynomial.__init__(self, b, a)
		
		b, a = self.coefficients
		self._root = -b / a

	@staticmethod
	def from_kw_args(a, b):
		"""Initialize a linear polynomial from keyword arguments.  Said polynomial has the form:
		f(x | a, b) := ax + b, where a != 0"""
		return Linear(a, b)
		
	@staticmethod
	def from_pos_args(*args):
		""" Initialize a linear polynomial of the form:
		f(x | a, b) := ax + b, where a != 0. 
		Note that 'b' = args[0] and 'a' = args[1]."""
		if len(args) != 2: raise IndexError('Argument Count Error:\tA linear polynomial must be initialized with exactly 2 coefficients.')
		return Linear(a=args[1], b=args[0])

	@property	
	def root(self):
		return self._root

	def _binom(self, n: int):
		""" Invoke binomial theorem to raise this polynomial to the n'th power, where 'n' is a non-negative whole number."""
		if n < 0: raise ValueError('Argument Value Error:\tRaising a polynomial to a negative power is prohibited (atleast for now).')
		from num_theory import binomial_coefficient as bc
		deg = n
		b, a = self.coefficients
		return Polynomial(*reversed(tuple(bc(n,m)*a**(n-m)*b**(m) for m in range(deg + 1))))

	__pow__ = _binom


class Quadratic (Polynomial):
	""" Quadratic Polynomial 
	A quadratic polynomial is of the form:  f(x | a, b, c) := ax^2 + bx + c, where a != 0."""
	
	def __init__ (self, a, b, c):
		""" Initialize a Quadtratic Polynomial of the form:
		f(x | a, b, c) := ax^2 + bx + c, where a != 0."""
		# Check parameters.
		if a == 0: raise ValueError('Argument Value Error:\tThe leading coefficent of a polynomial cannot be equal to zero.')
		Polynomial.__init__(self, c, b, a)
		self._roots = roots(a, b, c)
		self._extremum = extremum(a, b, c)

	@staticmethod
	def from_pos_args(*args):
		if len(args) != 3: raise IndexError('Argument Count Error:\tA quadratic or degree two polynomial must be initialized with exactly three coefficients.')
		return Quadratic(args[2], args[1], args[0])
		
	
	@staticmethod
	def from_props(ex_x, ex_y, w, y = 0):
		""" Create a quadratic polynomial from a set of properties instead of coefficients.
		The following properties of the parabola given by the polynomial to be created are used:
		1.)  (ex_x, ex_y):  The x, y coordinates of the extremum.  This is the point where the slope of the curve is zero.
		2.)  w(y):  The width of the parabola at any valid y-coordinate, ‘y’."""
		
		a = (4*y - 4*ex_y) / (w**2)
		b = -2*a*ex_x
		c = ex_y + a * ex_x**2
		
		return Quadratic(a, b, c)
	
	@property	
	def roots(self):
		__doc__ = roots.__doc__
		return self._roots
		
	@property
	def extremum (self):
		__doc__ = extremum.__doc__
		return self._extremum
		
	def width (self, y = 0):
		""" Caclulates the width of the parabola given by the quadratic polynomial at the y-coordinate, 'y'.
		When 'y' equals zero, the distance between the two roots (if they exist) is returned.  
		For example, width(0) = max(R) - min(R), where 'R' is the set of the roots of the polynomial.  
		Zero is returned if only one root exists, and an exception is raised if none exist."""
		if self.coefficients[2] < 0 and y > self.extremum[1]: raise ValueError('Argument Value Error:\tWhen the leading coefficient of a quadratic polynomial, \
			\'f\', is negative, \'y\' cannot be greater than the max(f).')
		elif self.coefficients[2] > 0 and y < self.extremum[1]: raise ValueError('Argument Value Error:\tWhen the leading coefficient of a quadratic polynomial, \
			\'f\', is positive, \'y\' cannot be less than the min(f).')
		elif y == self.extremum[1]: return 0
		
		a = self.coefficients[2]
		b = self.coefficients[1]
		c = self.coefficients[0]
		
		return abs(sqrt(b**2 + 4*a*(y-c)) / a)

class Cubic (Polynomial):
	"""A cubic polynomial is a third degree polynomial of the form:
	f(x | a, b, c, d) := ax^3 + bx^2 + cx + d, where a != 0."""

	def __init__(self, *args):
		if len(args) != 4: raise IndexError('Argument Count Error:\tA cubic or 3rd degree polynomal must be initialized with exactly four coefficients, even if they are \
			equal to zero.')
		Polynomial.__init__(self, *args)

		self._critical_points = self._find_critical_pnts()
		self._inflection_pnt = self._find_inflection_pnt()

	@staticmethod
	def from_props(cp0x, cp0y, cp1x, cp1y):
		"""Initialize a Cubic Polynomial from the x,y coordinates of a pair of critical points.  
		This is very useful in interpolation.  See /doc/CubicFromProps for more info.
		Parameters:
			'cp0x' is the x-coordinate of the zero'th critical point.
			'cp' stands for and denotes 'critical point';
			'0' equals the index of the critical point;
			And 'x' is the axis of the coordinate."""

		# Calculate the x-coordinate of inflection point, denoted as IPx.
		IPx = (cp1x + cp0x)/2

		# Calculate the leading coefficient, 'a', and then the remaining coefficients.
		a = (cp1y - cp0y) / (cp1x**3 - cp0x**3 - 3*IPx*cp1x**2 + 3*cp0x*cp1x**2 - 3*cp1x*cp0x**2 + 3*IPx*cp0x**2)
		b = -3*a*IPx
		c = 3*a*cp0x*cp1x
		d = cp0y - a*cp0x**3 - 3*a*cp1x*cp0x**2 + 3*a*IPx*cp0x**2

		return Cubic(d, c, b, a)

	@property
	def critical_points(self):
		return self._critical_points

	@property
	def inflection_point(self):
		return self._inflection_pnt

	def _find_critical_pnts(self):
		""" A 'private' method to find the critical points of the curve given by this polynomial.
		Returns an n-tuple of x,y coordinates, where 'n' is the number of roots in the derivative of this polynomial and 'n' in {0, 1, 2}."""
		# Find the roots of the derivative, denoted as 'p2'.  These are the x-coordinates of each critical point.
		df = self.derivative()
		p2 = Quadratic(df.coefficients[2], df.coefficients[1], df.coefficients[0])
		return tuple((r, self.f(r)) for r in p2.roots)

	def _find_inflection_pnt(self):
		d2f = self.derivative(2)
		p1 = Linear.from_pos_args(*d2f.coefficients)
		return (p1.root, self.f(p1.root))


def find_cubic_linearity(f, Ix: float):
	"""Let the argument, 'f', be a cubic polynomial and 'l' be defined as a linear polynomial of the form:
	l(x | m, s) := mx + s, where m != 0.

	Suppose that 'l' intesects 'f' at a point with an x-coordinate that is equal to the argument, 'Ix'.
	As long as the point of intersection, ‘I’, is not equal to the inflection point of ‘f’, ‘IP’, 
	it is possible to find  the coefficients of 'l' such that there is exactly one other point of incidence between ‘f’ and ‘l’ and ‘l’ is tangential to ‘f’ at that point.

	If 'I == IP' then the only point at which ‘l’ is tangential to ‘f’ is the inflection point;  In that case ‘f’ minus ‘l’ is a complete cube of the form, '(ax+b)^3',
	where a != 0.
	
	This function will find the coefficients of 'l' and return them in the form of an instance of the polynomial.Linear class.
	
	See doc/CubicLinearities.odt for more info."""
	pass

	# Check args.
	if isinstance(f, Polynomial) == False: raise TypeError(f'Argument Type Error:\tOnly instances of the Polynomial class may be passed to the \'f\' parameter of {__name__}.')
	if len(f) != 4: raise ValueError(f'Only degree 3 or cubic polynomials may be passed to \'f\'.')
	
	# Definitions
	d, c, b, a = f.coefficients
	p = Ix

	# 'C' is the coefficients of the quotient of dividing 'f-l' by 'x-p', exluding the leading coefficient which requires 'm' to be known.
	C = (a, a*p+b)

	# Solve for 'm' and 's'.
	m = a*p**2 + p*b + c - C[1]**2 / (4*C[0])
	s = f(p) - m*p

	return Linear(m, s)
		
		
		
		
		
	
		
	
