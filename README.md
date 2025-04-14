# Finding all WVGs
The folder 'storage' contains file all_wvgs.py with all WVGs for n up to (including) 6, and file all_wvgs_at_quota.py with all WVGs for n up to (including) 6 at given quotas. The set of quotas considered is from 0.5 to 0.95 in increments of 0.05.

These files were generated using:
find_all_wvgs.py: Finds all WVGs for n up to (including) 6. It ensures all are found by comparing them to the results from 'Enumeration of Threshold Functions of Eight Variables' by MUROGA, TSUBOI, BAUGH
find_all_wvgs_at_quota.py: Finds all WVGs for n up to (including) 6 and a given quota. It iterates through all WVGs from above and solves an LP to decide whether they can be achieved at the given quota. 

# Biggest Veto Distortion
The file biggest_veto_distortion.py estimates the greatest veto distortion for given n and quota by repeatedly generating random population targets and seeing if the closest banzhaf/shapley vectors induce veto players. The resulting plots get saved in the plots folder.

# Config
EXACT determines whether the code uses exact rational numbers or inexact floats (which make the code faster)
STRICT determines whether the coalition weights needs to be strictly or not strictly greater than the treshold to make the coalition win. I'm not sure if I implemented everything for the non-strict case.
