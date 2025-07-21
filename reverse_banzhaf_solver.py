"""
Reverse Limit Banzhaf Solver

This module provides a simple function to solve the reverse problem:
Given a target Banzhaf power distribution, find the weights that produce it in the limit.

Author: AI Assistant
Date: 2024
"""

import numpy as np
import sympy as sp
from scipy.optimize import root
from helpers import make_distribution

def get_saddle_point(weights, q, guess=1.0):
    """Find the saddle point for the limit Banzhaf computation."""
    weights = make_distribution(weights, False)
    x = sp.symbols('x', positive=True)
    expr = sum(w * x**w / (1 + x**w) for w in weights) - q * sum(weights)
    try:
        return float(sp.nsolve(expr, guess))
    except:
        return guess

def get_limit_banzhaf_powers(weights, q):
    """Compute limit Banzhaf powers for given weights and quota."""
    s = get_saddle_point(weights, q)
    print(s)
    return make_distribution([(1-s**w)/(1+s**w) for w in weights], False)

def solve_reverse_limit_banzhaf(target, quota=3/4, tolerance=1e-8, verbose=False):
    """
    Find weights x such that limit_banzhaf_powers(x, quota) ≈ target.
    
    Parameters:
    - target: Target limit Banzhaf power distribution (will be normalized)
    - quota: Voting quota (default 3/4)
    - tolerance: Convergence tolerance
    - verbose: Whether to print progress information
    
    Returns:
    - weights: The weight vector that produces the target powers in the limit (None if failed)
    - residual: Final residual ||f(weights) - target||
    - success: Whether the solver succeeded
    
    Example:
    >>> weights, residual, success = solve_reverse_limit_banzhaf([1, 2, 3])
    >>> if success:
    >>>     print(f"Weights: {weights}")
    >>>     print(f"Achieved powers: {get_limit_banzhaf_powers(weights, 3/4)}")
    """
    
    # Normalize target to probability distribution
    target_normalized = make_distribution(target, exact=False)
    
    if verbose:
        print(f"Finding weights for target Banzhaf powers: {target_normalized}")
        print(f"Using quota: {quota}")
    
    # Define the function we want to find the preimage of
    f = lambda x: get_limit_banzhaf_powers(x, quota)
    
    # Define residual function for root finding
    def residual_func(x):
        x = np.abs(x)  # Ensure positive domain
        try:
            fx = f(x)
            fx_array = np.array([float(val) for val in fx])
            target_array = np.array([float(val) for val in target_normalized])
            return fx_array - target_array
        except Exception:
            return np.full(len(target_normalized), 1e6)
    
    # Use target as initial guess (with small positive offset)
    initial_guess = np.array(target_normalized) + 1e-6
    
    # Try different numerical methods
    methods = ['lm', 'hybr', 'broyden1']
    
    best_solution = None
    best_residual = float('inf')
    
    for method in methods:
        if verbose:
            print(f"  Trying method: {method}")
        
        try:
            sol = root(residual_func, initial_guess, method=method, tol=tolerance)
            
            if sol.success:
                x_solution = np.abs(sol.x)
                
                # Compute final residual
                final_fx = f(x_solution)
                final_fx_array = np.array([float(val) for val in final_fx])
                target_array = np.array([float(val) for val in target_normalized])
                residual_norm = np.linalg.norm(final_fx_array - target_array)
                
                if residual_norm < best_residual:
                    best_solution = x_solution
                    best_residual = residual_norm
                    
                    if verbose:
                        print(f"    → Success! Residual: {residual_norm:.2e}")
                    
                    # If we found a very good solution, stop early
                    if residual_norm < tolerance * 10:
                        break
                        
            elif verbose:
                print(f"    → Method {method} failed to converge")
                
        except Exception as e:
            if verbose:
                print(f"    → Method {method} failed with error: {str(e)[:50]}...")
    
    success = best_solution is not None
    
    if verbose:
        if success:
            print(f"\n✓ Solution found!")
            print(f"  Weights: {best_solution}")
            print(f"  Final residual: {best_residual:.2e}")
            
            # Verify the solution
            achieved_powers = f(best_solution)
            print(f"  Achieved powers: {[float(val) for val in achieved_powers]}")
            print(f"  Target powers:   {target_normalized}")
        else:
            print(f"\n✗ No solution found. All methods failed.")
    
    return best_solution, best_residual, success

# Additional utility function for batch solving
def solve_multiple_targets(targets_list, quota=3/4, tolerance=1e-8, verbose=True):
    """
    Solve reverse limit Banzhaf problem for multiple targets.
    
    Parameters:
    - targets_list: List of target power distributions
    - quota: Voting quota
    - tolerance: Convergence tolerance
    - verbose: Whether to print progress
    
    Returns:
    - results: Dictionary with results for each target
    """
    if verbose:
        print(f"=== Solving reverse Banzhaf for {len(targets_list)} targets ===\n")
    
    results = {}
    
    for i, target in enumerate(targets_list):
        if verbose:
            print(f"Target {i+1}: {target}")
        
        weights, residual, success = solve_reverse_limit_banzhaf(
            target, quota, tolerance, verbose=False
        )
        
        results[i] = {
            'target': target,
            'weights': weights,
            'residual': residual,
            'success': success
        }
        
        if verbose:
            if success:
                print(f"  ✓ Success! Weights: {weights}")
                print(f"    Residual: {residual:.2e}")
            else:
                print(f"  ✗ Failed")
            print()
    
    return results

if __name__ == "__main__":
    # Quick test
    print("Testing reverse Banzhaf solver...")
    target = [1, 2, 3]
    weights, residual, success = solve_reverse_limit_banzhaf(target, verbose=True)
    
    if success:
        print(f"\nVerification:")
        achieved = get_limit_banzhaf_powers(weights, 3/4)
        print(f"Target:   {make_distribution(target, exact=False)}")
        print(f"Achieved: {[float(x) for x in achieved]}")
        print(f"Error:    {residual:.2e}") 