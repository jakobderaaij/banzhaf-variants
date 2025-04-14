import gurobipy as gp
from helpers import *
from config import EXACT, STRICT
from find_all_wvgs import to_function, get_equivalent_players
from storage.all_wvgs import all_wvgs as ALL_WVGS
ALL_INDICES = dict()

CORRECT_NUMBERS_AT_HALF = {1:1, 2:2, 3:4, 4:12, 5:81, 6:1684, 7:123565} # from Enumeration of Threshold Functions of Eight Variables by MUROGA, TSUBOI, BAUGH

def function_to_winning_coalitions(func, to_loosing = False):
    winning = []
    n = int(math.log2(len(func)))
    if n not in ALL_INDICES: ALL_INDICES[n] = powerset(list(range(n)))
    for i,indices in enumerate(ALL_INDICES[n]):
        if not to_loosing:
            if func[i] == "1":
                winning.append(indices)
        else:
            if func[i] == "0":
                winning.append(indices)
    return winning

def is_subset_of(subset,set_):
    for ele in subset:
        if ele not in set_:
            return False
    return True

def get_minimal_winning_coalitions(winning):
    critical = []
    for coalition in winning:
        is_critical = True
        for other_coalition in winning:
            if other_coalition != coalition and is_subset_of(other_coalition, coalition):
                is_critical = False
        if is_critical: critical.append(coalition)
    return critical

def get_maximal_loosing_coalitions(loosing):
    critical = []
    for coalition in loosing:
        is_critical = True
        for other_coalition in loosing:
            if other_coalition != coalition and is_subset_of(coalition, other_coalition):
                is_critical = False
        if is_critical: critical.append(coalition)
    return critical

if __name__ == "__main__":
    LOWEST_N = 1
    HIGHEST_N = 6
    QUOTA_RANGE = [i/sp.Rational(20) for i in range(10,20)] + [2/sp.Rational(3)]
    N_RANGE = list(range(LOWEST_N, HIGHEST_N + 1))
    all_wvgs_at_quota = dict()
    for n in N_RANGE:
        print(f"n={n}")
        all_wvgs_at_quota[n] = dict()
        for quota in QUOTA_RANGE:
            all_wvgs_at_quota[n][quota] = dict()

            quota_numer = sp.numer(quota)
            quota_denom = sp.denom(quota)

            if n not in ALL_WVGS: raise NotImplementedError
            wvgs = ALL_WVGS[n]
            
            for function, (weights_here, quota_here) in wvgs.items():
                if not EXACT: raise NotImplementedError

                # Check if this function is attainable at this quota:
                minimal_winning_coalitions = get_minimal_winning_coalitions(function_to_winning_coalitions(function))
                maximal_loosing_coalitions = get_maximal_loosing_coalitions(function_to_winning_coalitions(function, to_loosing=True))

                m = gp.Model()
                m.setParam('OutputFlag', 0)  # Suppress output to speed up batch solving
                # print(weights, quota_here/sum(weights))
                # print(banzhaf(weights,quota_here))
                
                slack = m.addVar(name=f"slack", vtype=gp.GRB.CONTINUOUS)
                weight_vars = []
                for player in range(n):
                    weight_vars.append(m.addVar(name=f"weight{player}", vtype=gp.GRB.CONTINUOUS))
                m.update()

                total_weight = gp.LinExpr()
                for weight_var in weight_vars:
                    total_weight += weight_var
                m.addConstr(total_weight == 1, 'total_weight')

                m.addConstr(slack >= 0, 'non-negative_slack')

                for i, winning_coalition in enumerate(minimal_winning_coalitions):
                    winning_set = gp.LinExpr()
                    for player in range(n):
                        if player in winning_coalition:
                            winning_set += weight_vars[player]
                    m.addConstr(quota_denom * winning_set - slack >= quota_numer, "winning_coaliltion_{i}") # TODO: The + 0.000001 is a bit hacky, since guorbi doesn't allow strict inequalities
                
                for i, loosing_coalition in enumerate(maximal_loosing_coalitions):
                    loosing_set = gp.LinExpr()
                    for player in range(n):
                        if player in loosing_coalition:
                            loosing_set += weight_vars[player]
                    m.addConstr(quota_denom * loosing_set + slack <= quota_numer, "loosing_coaliltion_{i}") 
                
                m.setObjective(slack, gp.GRB.MAXIMIZE)

                # Solve the model
                m.optimize()

                # Process results
                found_weights = []
                if m.status == gp.GRB.OPTIMAL:

                    for v in m.getVars():
                        if v.VarName =='slack':
                            slack_value = v.X  
                        else:
                            found_weights.append(v.X)
                    if STRICT and slack_value == 0: 
                        continue

                    if slack_value < 0.001:
                        raise NotImplementedError #this shouldn't happen
                    
                    found_weights.sort()

                    # Sanity Check
                    assert to_function(found_weights, quota) == function
                    
                    all_wvgs_at_quota[n][quota][function] = tuple(found_weights)   
        output_string = "Found: "
        for quota in QUOTA_RANGE:
            output_string+= f'({quota}: {len(all_wvgs_at_quota[n][quota])}), '
        print(output_string)   

    # for sanity check
    Factorials = {k : math.factorial(k) for k in range(0,HIGHEST_N+1)}
    quota = 1/sp.Rational(2)
    for n in N_RANGE:
        if quota in all_wvgs_at_quota[n]:
            number_of_function_with_multiplicity = 0 
            for function, weights in all_wvgs_at_quota[n][quota].items():
                        equivalent_players = get_equivalent_players(weights, quota, STRICT)
                        multiplicacy = Factorials[n]
                        for equivalence_set in equivalent_players.values():
                            multiplicacy= multiplicacy/sp.Rational(Factorials[len(equivalence_set) + 1]) 
                        number_of_function_with_multiplicity += multiplicacy
            assert number_of_function_with_multiplicity == CORRECT_NUMBERS_AT_HALF[n]

    # Write all_wvgs_at_quota to file
    with open(f'storage/all_wvgs_at_quota.py', "w") as file:
        file.write("import sympy as sp\n")
        file.write("all_wvgs_at_quota = {\n")
        for n in N_RANGE:
            file.write(f'\t{n}: '+'{\n')
            for quota in QUOTA_RANGE:
                file.write(f'\t\t{sp.numer(quota)}/sp.Rational({sp.denom(quota)}): '+'{\n')
                for func, weights in all_wvgs_at_quota[n][quota].items():
                    file.write(f'\t\t\t"{func}": {weights},\n')
                file.write('\t\t},\n')
            file.write('\t},\n')
        file.write('}')