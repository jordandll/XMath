def int_exp_sin(a, w):
	""" Finds the primitive integral of the function, 
		f(t | a, w) := C e^(at + b) sin(wt + p)
	where C, b, and p are arbitrary constants.  Said primitive integral, hencforth 'A_f', has the general form:
	A_f = Ce^(at + b)(A sin(wt + p) + B cos(wt + p)).
	Returns 'A' and 'B' found in the above identity."""
	
	B = - w / (a**2 + w**2)
	A = (1 + w*B) / a
	
	return (A, B)
	
def int_exp_cos (a, w):
	""" Finds the primitive integral of the function, 
		f(t | a, w) := C e^(at + b) cos(wt + p)
	where C, b, and p are arbitrary constants.  Said primitive integral, hencforth 'A_f', has the general form:
	A_f = Ce^(at + b)(A sin(wt + p) + B cos(wt + p)).
	Returns 'A' and 'B' found in the above identity."""
	
	A = w / (a**2 + w**2)
	B = (1 - w*A)/a
	
	return (A, B)
	
