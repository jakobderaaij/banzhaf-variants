# If True, winning coalitions need to have weight strictly greater than the threshold.
# If False, only greater or equal.
STRICT = True 

# If True, uses sp.Rational for exact calculations.
# If False, uses floats (faster but not exact).
EXACT = True 


if STRICT == False: raise NotImplementedError