import math
import sympy as sp
import numpy as np
from config import EXACT, STRICT

def powerset(list_):
    if len(list_) == 0: return [[]]
    rec = powerset(list_[:-1])
    return rec + [ele+[list_[-1]] for ele in rec] 

def coalition_weight(coalition, weights):
    ''' Coalition is set/list of indices'''
    return sum([weights[i] for i in coalition])

def make_distribution(vector, exact = EXACT):
    tot = sum(vector)
    if exact: return [ele/sp.Rational(tot) for ele in vector]
    return [ele/(tot) for ele in vector]

def round_vector(vector, digits = 3):
    return [round(float(ele),digits) for ele in vector]

def dot_product(vector, vector2):
    return sum([ele1 * ele2 for ele1,ele2 in zip(vector,vector2)])

def distance(vec1, vec2, p = 1, exact = EXACT):
    if exact:
        distances = [abs(sp.Rational(ele1)-sp.Rational(ele2)) for ele1, ele2 in zip(vec1, vec2)]
        p = sp.Rational(p)
    else:
        distances = [abs(ele1-ele2) for ele1, ele2 in zip(vec1, vec2)]

    if p == "infty": return max(distances)
    return sum([ele**p for ele in distances])**(1/p)

def compare_print(benchmark, colored):
    print('[',end='')
    for i, (banzhaf_val, pop_val) in enumerate(zip(round_vector(colored), round_vector(benchmark))):
        if banzhaf_val < pop_val:
            print(f"\033[32m{banzhaf_val}\033[0m", end="")  # Green
        elif banzhaf_val > pop_val:
            print(f"\033[31m{banzhaf_val}\033[0m", end="")  # Red
        else:
            print(f"\033[33m{banzhaf_val}\033[0m", end="")  # Yellow
        if i < len(benchmark)-1:
            print(", ", end="")
    print(']')

def generalized_banzhaf(population, quota, decisivness, normalize = True, exact = EXACT, strict = STRICT):
    n = len(population)
    distribution = [math.comb(n-1, i) * decisivness**i * (1-decisivness)**(n-1-i) for i in range(0,n)]
    return semivalue(population, quota, distribution, normalize, exact, strict)

def semivalue(population, quota, distribution, normalize = True, exact = EXACT, strict = STRICT):
    n = len(population) 
    results = [0]*n
    coalitions = powerset(list(range(n)))
    for coalition in coalitions:
        k = len(coalition)
        weight =  coalition_weight(coalition, population)
        if not strict:
            if weight < quota:
                for person, person_weight in enumerate(population):
                    if person not in coalition:
                        if weight + person_weight >= quota:
                            if exact:
                                results[person] += distribution[k]/sp.Rational(math.comb(n-1, k))
                            else: 
                                results[person] += distribution[k]/(math.comb(n-1, k))
        else:
            if weight <= quota:
                for person, person_weight in enumerate(population):
                    if person not in coalition:
                        if weight + person_weight > quota:
                            if exact:
                                results[person] += distribution[k]/sp.Rational(math.comb(n-1, k))
                            else: 
                                results[person] += distribution[k]/(math.comb(n-1, k))
    index = results

    if normalize: return make_distribution(index, exact)
    return index

    # distribution = make_distribution(distribution)
    # population = [(i,v) for i,v in enumerate(population)]
    # n = len(population) 
    # results = [0]*n
    # coalitions = powerset(population)
    # for coalition in coalitions:
    #     k = len(coalition)
    #     if weight(coalition) < quota:
    #         people_in_coalition = [i for i,_ in coalition]
    #         for person, person_weight in population:
    #             if person not in people_in_coalition:
    #                 if weight(coalition) + person_weight >= quota:
    #                     results[person] += (distribution[k] / (math.comb(n-1, k)))
    # return make_distribution(results)

def banzhaf(population, quota, normalize = True, exact = EXACT, strict = STRICT):
    n = len(population) 
    results = [0]*n
    coalitions = powerset(list(range(n)))
    for coalition in coalitions:
        weight =  coalition_weight(coalition, population)
        if not strict:
            if weight < quota:
                for person, person_weight in enumerate(population):
                    if person not in coalition:
                        if weight + person_weight >= quota:
                            results[person] += 1
        else:
            if weight <= quota:
                for person, person_weight in enumerate(population):
                    if person not in coalition:
                        if weight + person_weight > quota:
                            results[person] += 1
    if exact: index = [ele/sp.Rational(2**(n-1)) for ele in results]
    else: index = [ele/(2**(n-1)) for ele in results]

    if normalize: return make_distribution(index, exact)
    return index

def shapley(population, quota, exact = EXACT, strict = STRICT):
    n = len(population) 
    Factorials = {k : math.factorial(k) for k in range(0,n)}
    results = [0]*n
    coalitions = powerset(list(range(n)))
    for coalition in coalitions:
        k = len(coalition)
        weight =  coalition_weight(coalition, population)
        if not strict:
            if weight < quota:
                for person, person_weight in enumerate(population):
                    if person not in coalition:
                        if weight + person_weight >= quota:
                            results[person] += Factorials[k]*Factorials[n-k-1]
        else: 
            if weight <= quota:
                for person, person_weight in enumerate(population):
                    if person not in coalition:
                        if weight + person_weight > quota:
                            results[person] += Factorials[k]*Factorials[n-k-1]
    if exact: index = [ele/sp.Rational(math.factorial(n)) for ele in results]
    else: index = [ele/(math.factorial(n)) for ele in results]

    return index
   
def print_progress(progress, trials, msg = None):
    progress = progress + 1
    if (progress*1000)%trials == 0:
        output = ""
        if msg: output += f'{msg}: '
        output += str((progress*100)/trials) + "%\r"
        print(output,end='')
    if progress == trials:
        print("Done.    ")
    
def find_closest(target, candidates, p = 1, exact = EXACT):
    min_vec = None
    min_distance = len(target)
    for vec in candidates:
        if distance(vec, target, p, exact) < min_distance:
            min_distance = distance(vec, target, p, exact)
            min_vec = vec
    return min_vec

# def get_pivotal_vectors(weights, quota, exact = EXACT, invert = False, normalize = True):
#     n = len(weights)
#     pivotal_vectors = [[0  for _ in range(n)] for _ in range(n)]
#     coalitions = powerset(list(range(n)))
#     for coalition in coalitions:
#         k = len(coalition)
#         for player in coalition:
#             weight = coalition_weight(coalition, weights)
#             if weight >= quota and weight - weights[player] < quota:
#                 if invert: 
#                     pivotal_vectors[k-1][player] += 1
#                 else:
#                     pivotal_vectors[player][k-1] += 1

#     if not normalize:
#         for player in range(n):
#             for size in range(n):
#                 if invert:
#                     if exact:
#                         pivotal_vectors[size][player] *= 1/sp.Rational(math.comb(n-1, size))
#                     else:
#                         pivotal_vectors[size][player] *= 1/(math.comb(n-1, size))
#                 else: 
#                     if exact:
#                         pivotal_vectors[player][size] *= 1/sp.Rational(math.comb(n-1, size))
#                     else:
#                         pivotal_vectors[player][size] *= 1/(math.comb(n-1, size))
#     else:
#         for size in range(n):
#             if invert: total = sum([pivotal_vectors[size][i] for i in range(n)])
#             else: total = sum([pivotal_vectors[i][size] for i in range(n)])
#             if total == 0: continue
#             for player in range(n):
#                 if invert:
#                     if exact:
#                         pivotal_vectors[size][player] *= 1/sp.Rational(total)
#                     else:
#                         pivotal_vectors[size][player] *= 1/total
#                 else: 
#                     if exact:
#                         pivotal_vectors[player][size] *= 1/sp.Rational(total)
#                     else:
#                         pivotal_vectors[player][size] *= 1/(total)
    
#     return pivotal_vectors

# def get_worst_distance(target, population, quota, p, exact = EXACT):
#     n = len(target)
#     worst_distance = 0
#     worst_distance_weights = None
#     pivotal_vectors = get_pivotal_vectors(population, quota, exact, True)
#     for vector in pivotal_vectors[:(n+1)//2]:
#         if sum(vector) == 0: continue
#         vector = make_distribution(vector, exact)
#         dis = distance(target, vector, p, exact)
#         if dis > worst_distance:
#             worst_distance = dis
#             worst_distance_weights = vector
#     if not worst_distance_weights: print("Something unexpected happened", target, population, vector)
#     return worst_distance, worst_distance_weights


if __name__ == "__main__":
    print("Running tests:")
    population = [3,2,1,1]
    quota = 4
    assert shapley(population, 4, True) == [1/sp.Rational(2), 1/sp.Rational(6), 1/sp.Rational(6), 1/sp.Rational(6)]
    print("Shapley Success")
    population = [4,3,2,1]
    quota = 6
    assert banzhaf(population, quota, False, True) == [5/sp.Rational(8), 3/sp.Rational(8), 3/sp.Rational(8), 1/sp.Rational(8)]
    assert banzhaf(population, quota, True, True) == [5/sp.Rational(12), 3/sp.Rational(12), 3/sp.Rational(12), 1/sp.Rational(12)]
    print("Banzhaf Success")

    # assert get_pivotal_vectors([5,5,3,3,3], 19/sp.Rational(2), True, False, False) == [[0, 1/sp.Rational(4), 1, 1/sp.Rational(4), 0], [0, 1/sp.Rational(4), 1, 1/sp.Rational(4), 0], [0, 0, 2/sp.Rational(3),0, 0], [0, 0, 2/sp.Rational(3),0, 0], [0, 0, 2/sp.Rational(3),0, 0]]
    # assert get_pivotal_vectors([5,3,2,1], 11/sp.Rational(2), True, False, False) == [[0, 1, 1, 0], [0, 1/sp.Rational(3), 1/sp.Rational(3), 0], [0, 1/sp.Rational(3), 1/sp.Rational(3), 0], [0, 1/sp.Rational(3), 1/sp.Rational(3), 0]]
   
    weights = np.random.dirichlet(alpha=[1]*6)
    # assert get_pivotal_vectors(weights, 1/2, True, True, False) == [semivalue(weights, 1/2, [0]*(i) + [1] + [0]*(5-i), False, True) for i in range (6)]
    print("Semivalue success")

    assert banzhaf(weights, 1/sp.Rational(2), True, True) == generalized_banzhaf(weights,1/sp.Rational(2), 1/sp.Rational(2), True, True)
    assert banzhaf(weights, 2/sp.Rational(3), True, True) == generalized_banzhaf(weights,2/sp.Rational(3), 1/sp.Rational(2), True, True)
    assert [0,0,0,0,0] == generalized_banzhaf([1,1,1,1,1],1/sp.Rational(2), 1, False, True)
