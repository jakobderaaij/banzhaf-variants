from helpers import *
from config import EXACT, STRICT
from storage.all_wvgs import all_wvgs
import tqdm


CORRECT_NUMBERS = {0:2, 1:3, 2:6, 3:20, 4:150, 5:3_287, 6:244_158, 7:66_291_591} # from Enumeration of Threshold Functions of Eight Variables by MUROGA, TSUBOI, BAUGH


# I found this too late, would have been more useful
# CORRECT_NUMBERS_WITHOUT_PERMUTATIONS = {0:2, 1:3, 2:5, 3:10, 4:27, 5:119, 6:1_113, 7:29_375} # from Enumeration of Threshold Functions of Eight Variables by MUROGA, TSUBOI, BAUGH

def all_ordered_wvgs(n, largest_weight):
    ''' Get all WVGS with n players with weights between 0 and largest weight, with non-decreasing weights.'''
    if n == 1: return [[i] for i in range(largest_weight+1)]
    recursive = all_ordered_wvgs(n-1, largest_weight)
    all_wvgs = []
    for wvg in recursive:
        for weight in range(wvg[-1],largest_weight+1):
            all_wvgs += [wvg + [weight]]
    return all_wvgs




def get_equivalent_players(weights, quota, strict = STRICT):
    n = len(weights)
    p1 = 0
    equivalent_players = {p1: []}
    for i in range(1,n):
        p2 = i
        if weights[p1] == weights[p2]: 
            equivalent_players[p1].append(p2)
            continue
        all_other_indices = powerset([j for j in range(n) if j not in [p1,p2]])
        equal = True

        for indices in all_other_indices:
            if strict:
                if (coalition_weight([p1]+indices,weights)>quota) ^ (coalition_weight([p2]+indices,weights)>quota):
                    equal = False
            else:
                if (coalition_weight([p1]+indices,weights)>=quota) ^ (coalition_weight([p2]+indices,weights)>=quota):
                    equal = False
        
        if equal: 
            equivalent_players[p1].append(p2)
        else:
            p1 = p2
            equivalent_players[p1] = []
    return equivalent_players

if __name__ == '__main__':
    LOWEST_N = 1
    HIGHEST_N = 6
    MAXIMUM_WEIGHTS = {1:1, 2:1, 3:2, 4:3, 5:5, 6:9, 7:20}
    for n in range(LOWEST_N, HIGHEST_N + 1):
        # Only run this if n is not in the list yet
        if n in all_wvgs: continue

        Factorials = {k : math.factorial(k) for k in range(0,n+1)}

        possible_functions_with_wvs = set()
        number_of_function_with_multiplicity = 0
        seen = set()
        for weights in tqdm.tqdm(all_ordered_wvgs(n, MAXIMUM_WEIGHTS[n])): 
            if sum(weights)==0: continue
            for q in range(sum(weights)):
                quota = q + 0.5
                func = to_function(weights, quota, STRICT)
                if func not in seen:
                    seen.add(func)
                    possible_functions_with_wvs.add((func, (tuple(weights), quota)))

                    # Get multiplicacy
                    equivalent_players = get_equivalent_players(weights, q, STRICT)
                    multiplicacy = Factorials[n]
                    for equivalence_set in equivalent_players.values():
                        multiplicacy= multiplicacy/sp.Rational(Factorials[len(equivalence_set) + 1])
                    number_of_function_with_multiplicity += multiplicacy

        print(f'Found {len(possible_functions_with_wvs)} unique ordered WVGs for n={n}. Including permutations, found {number_of_function_with_multiplicity}')
        if CORRECT_NUMBERS[n] != number_of_function_with_multiplicity +2: # +2 because the paper also considers the function where all coalitions (cinluding the empty one) are winning and where no coalition (including the grand coalition) is winning
            print(f"Error: maximum weight {MAXIMUM_WEIGHTS[n]} not big enough for n={n} - not all WVGs found.")
            print(f"Only found {number_of_function_with_multiplicity} of {CORRECT_NUMBERS[n]-2}.")
            continue

        with open(f'storage/all_wvgs.py', "r") as file:
            content = file.readlines()
        last_line = content.pop()
        if last_line != "}": raise ValueError
        with open(f'storage/all_wvgs.py', 'w') as file:
            file.writelines(content)
        with open(f'storage/all_wvgs.py', 'a') as file:
            file.write("\t"+str(n)+": {\n")
            for (func, weights_with_quota) in possible_functions_with_wvs:
                file.write(f'\t\t"{func}": {weights_with_quota} ,\n')
            file.write('\t},\n') 
            file.write('}')
