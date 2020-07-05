# Modular Exponentiation
def mod_exp(b: int, a: int, n: int) -> int:
    """Fast Modular Exponentiation:
    Raises 'm' to the 'a' power modulo 'n'.  
    This implementation has a computational complexity of O(log2(a))."""    
    res = 1
    for i in range(a.bit_length()):
        if bool((1<<i) & a): res = (res * b) % n
        b = (b**2) % n
    return res

# Multiplicative Inverse
def mul_inv(m: int, n: int) -> int:
	pass

def mul_group(n: int) -> set:
	"""" Multiplicative Group of Integers modulo 'n'.
	Returns the underlying set of the aforementioned group, which is simply the set integers less than 'n' that are coprime to 'n'."""
	return set(from_sieve(coprime_sieve(n, n)))

def mul_ord(a: int, n, cn = None) -> int:
	""" Order of 'a' in the multiplicative groups of integers modulo 'n'.
	Notes on parameters:
	1.)  'a' should be coprime to 'n'.
	2.)  If 'n' is too large to factor, then either 'n' or 'cn' should be a dictionary of prime factors and powers.
	3.)  'cn' is the Carmichael number of 'n'.  If left at 'None', it is computed by this function, but only if 'n' can be factored or already is."""
	if cn is None: cn = carmichael_num(n)
	if type(n) is not int:
		tmp = n
		n = 1
		for p, r in tmp.items(): n*= p**r
	# Store the prime factors of the Carmichael number, 'cn', in a set denoted as 'PF'.
	P = _factor(cn)
	
	if type(cn) is not int:
		cn = 1
		for p, r in P: cn *= p**r

	# Navigate the lattice of divisors of 'cn' ordered by divisibility, which is isomorphic to the lattice of subgroups of Zn/Z^X ordered by inclusion.
	
	# 'top' and 'bottom' refer to the current sublattice.
	# 'divisors' is a set of prime factors/divisors that are assigned to the dimensions of the current sublattice.
	# 'div_cur' is the set of divisors of the link(s) connecting the top of the current or parent sublattice to the top of the next or child sublattice.
	# When 'div_cur' is empty that indicates that the search is over and the current top is returned as the order.
	top, bottom = cn, 1
	divisors = set(P.keys())
	
	while top != bottom:
		# Check whether any element the top directly covers satisfies the congruence relation:  a**c % n = 1;  Where 'c' is the directly covered element.
		# If that relation is satisfied, add the divisor of the link connecting the top to 'c' to the set of current divisors, 'div_cur'.
		div_cur = set()
		for d in divisors:
			if mod_exp(a, top // d, n) == 1:
				div_cur.add(d)
			
		# Set the dimensions of the next sublattice.
		divisors = div_cur.copy()
		
		# Calculate the top and bottom of the next/child sublattice.
		top = top // from_set(div_cur)
		P = _factor(top); P[1] = 1
		bottom = from_dict({x: P[x] for x in set(P.keys()) - div_cur})
	
			
	return top

def from_set (X: set) -> int:
	""" Converts a set of factors into an integer.  This function merely returns the product of the numbers in 'x'."""
	r = 1
	for x in X: r *= x
	return r
	

def gcd (a, b , *args):
	""" Computes the greatest common divisor(GCD) of the arguments given."""
	if issubclass(type(a), dict):
		return gcd_pre_factored(a, b, *args)
	elif issubclass(type(a), int):
		return _gcd(a, b, *args)
	else: raise TypeError('Invalid Argument Type Error:\t\'gcd(...)\' only takes Integer and Mapping type arguments.')
	

def gcd_pre_factored (a: dict, b: dict, *args) -> dict:
	""" Computes the greatest common divisor(GCD) of the arguments given.
	This acts on dictionaries containing the prime factors and powers of an integer, a variation of the unique form that every integer has according to the 
	fundemental theorem of arithmetic."""
	if len(args) > 0:
		arg_list = [a, b]
		arg_list.extend(args)
		t = a
		for n in arg_list[1:]: t = _gcd_pre_factored(t, n)
	else: t = _gcd_pre_factored(a, b)
	return t
	
def _gcd_pre_factored (a: dict, b: dict) -> dict:
	# The GCD is given by the intersection of the two dictionaries of prime factors and powers, which is given ...
	# ... by performing elementwise logical conjuction on each prime factor to yield the new set of keys, coupled with performing elementwise fuzzy logical conjuction, given by 'min(...)', ... 
	# ... on each corresponding prime power, to yield the values associated with each key.
	return {p: min(a[p], b[p]) for p in set(a.keys()) & set(b.keys())}
	
	
def from_dict(x: dict) -> int:
	""" Converts a dictionary containing the prime factors and powers of an integer into it's integer form."""
	res = 1
	for p, r in x.items(): res *= p**r
	return res
	
def _gcd(a: int, b: int, *args) -> int:
	if len(args) > 0:
		arg_list = [a, b]
		arg_list.extend(args)
		t = a
		for n in arg_list[1:]: t = euclid_algo(t, n)
	else: t = euclid_algo(a, b)
	return t

def euclid_algo(a: int, b: int) -> int:
	"""Euclidean Algorithm:
	
	Computes the greatest common divisor(GCD) of 'a' and 'b' by 
	implementing the Euclidean algorithm.
	This function has a computational complexity of log10(max(a, b))*5."""	
	while b != 0:
		t = b
		b = a % b
		a = t
	return a
	
# Extended Euclidean Algorithm
def euclid_algo_ext(a: int, b: int) -> tuple:
	"""Extended Euclidean Algorithm:
	
	Computes the greatest common divisor of 'a' and 'b' as well as the
	Bezout's coefficients, 'x' and 'y', of the equation: 
	gcd(a, b) = ax + by
	
	Returns the previously mentioned values as a 3-tuple in the form of 
	(gcd(a, b), x, y) .
	
	NOTE:  If 'a' and 'b' are coprime and 'a' is less than 'b' then 'x' 
	is the multiplicative inverse of 'a' modulo 'b'.  And 'y' is the 
	multiplicative of 'b' modulo 'a' if 'b' is less than 'a'."""
	sa, sb = 1, 0
	ta, tb = 0, 1
	while b != 0:
		t = b
		q, b = divmod(a, b)
		a = t
		st = sb
		sb = sa - q*sb
		sa = st
		tt = tb
		tb = ta - q*tb
		ta = tt
	return (a, sa, ta)

def lcm(a: int, b: int, *args) -> int:
	"""Least Common Multiple:
	Computes the LCM of 'a' and 'b' and any other args passed to it.
	This function has a computational complexity of log10(max(a, b, *args))*5.  
	It simply divides the product of 'a' and 'b' by the gcd(a, b) and returns the result."""
	if len(args) > 0:
		arg_list = [a, b]
		arg_list.extend(args)
		t = 1
		for n in arg_list: t = _lcm(t, n)
	else: t = _lcm(a, b)
	return t
	
def _lcm(a: int, b: int) -> int:
	"""Least Common Multiple:
	Computes the LCM of 'a' and 'b'.
	This function has a computational complexity of log10(max(a, b))*5.  
	It simply divides the product of 'a' and 'b' by the gcd(a, b) and returns the result."""		
	_gcd = euclid_algo(a,b)
	return (a // _gcd) * b

def coprime_sieve(x, n: int = 100) -> list:
	""" Takes a number, 'x', and returns a sequence of boolean elements whose index equals a number coprime to 'x' if set to 'True'.
	The length of the list/sequence is capped at the argument,'n'.
	NOTE:  The argument, 'x', can be a dictionary of prime factors and powers, a set of prime factors, or an integer."""
	try:
		if type(x) is not set: x = set(_factor(x).keys())
	except TypeError as ex:
		# Check if 'x' is set-like.
		if len(x) != len(set(x)): raise TypeError('Argument Error:\tThe argument passed to \'x\' must be a mapping, integer, or set-like type.')
			
	l=[True for i in range(n)]
	l[0], l[1] = False, True
	for p in x:
		for i in range(p, n, p): l[i]=False
	return l
	
def from_sieve(l: list) -> list:
	""" Takes a sieve, a sequence of boolean elements whose index is some type of number -- e.g. a prime number -- if it is set to 'True', and returns that list of numbers. """
	p=list()
	for i, val in enumerate(l):
		if val: p.append(i)
	return p

def sieve(n: int = 100):
	""" Prime Sieve:
	Creates a sieve that can be used to generate prime numbers up to and including 'n'.  If 'n' is prime, it will indicated as such by this sieve."""
	n += 1
	# Create a sequence of boolean values whose index equals the number being checked for primality.  If the element at index 'i' is true, 'i' is prime.
	l=[True for i in range(n)]
	l[0]=False; l[1]=False
	# Enumerate the elements in 'l'.  'p' is the index and equals the number being checked for primality.  'val' is the booleans value that indicates that 'p' is prime if it is set to 'True'.
	for p, val in enumerate(l):
		if val:
		# 'p' is prime.
			for i in range(p*2, n, p): l[i]=False
	return l

def factor(n: int) -> dict:
	""" Finds the prime factors and powers of 'n' and returns them as a dictionary with the keys beings the set of prime factors and their associated values being the prime powers.
	This implementation has a computational complexity of O(sqrt(n)).  Note that there are faster algorithms, which should be implemented later."""
	# Generate a list of primes, 'p', up to the square root of 'n' rounded up.
	from math import ceil, floor
	P = from_sieve(sieve(floor(n**(1/2))))
	pf = dict()
	for p in P:
		if p > floor(n**(1/2)): break
		n, r = divmod(n, p)
		while r == 0:
			if p in pf.keys(): pf[p] += 1
			else: pf[p] = 1
			n, r = divmod(n, p)
		n = n*p + r
	if n not in pf.keys() and n != 1: pf[n] = 1
	return pf
	
def _factor(x) -> dict:
	"""This function takes either a dictionary of prime factors and powers or an integer and returns a dictionary of prime factors as key and prime powers as values.
	NOTE:  This funtion is for implementors of this module (as opposed to consumers of this module). """
	pf = dict()
	if type(x) == dict: pf = x
	elif type(x) == int: pf = factor(x)
	else: raise TypeError(f'Implementation Error:\t  \'_factor(x)\' only takes types of \'dict\' and \'int\'.')
	return pf
	
def eulers_totient(x) -> int:
	P = _factor(x)
	res = 1
	for pf, pp in P.items():
		res *= pf**(pp-1) * (pf-1)
	return res
	
def carmichael_num(x) -> int:
	P = _factor(x)
	# Store the carmichael totient values of each prime factor/power of 'x' in a list, denoted as 'tot_vals'.
	tot_vals = [eulers_totient({p: r}) // 2 if p == 2 and r > 2 else eulers_totient({p: r}) for p, r in P.items()]
	if len(tot_vals) == 1: return tot_vals[0]
	return lcm(*tot_vals)

def binomial_coefficient(n: int, m: int) -> int:
	""" Binomial Coefficient
	Returns n!/(m!(n-m)!).  This is used in combinatronics and binomial theorem."""
	from math import factorial as f
	return f(n)/(f(m)*f(n-m))

	

	
