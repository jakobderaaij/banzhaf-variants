from storage.all_wvgs_at_quota import all_wvgs_at_quota as ALL_WVGS_AT_QUOTA
from helpers import *
import sympy as sp
import numpy as np
from config import EXACT, STRICT
import matplotlib.pyplot as plt



p = 1
trials = {3:10}#, 4:10000, 5: 10000, 6: 10000}

indices = {'banzhaf': lambda weights, q: banzhaf(weights, q, True), 'shapley': shapley, "no_veto": no_veto_index}#'generalized_banzhaf': lambda weights, q: generalized_banzhaf(weights, q, q, True)}

QUOTA_RANGE = [i/sp.Rational(20) for i in range(10,20)] 


for n in trials.keys():
    print(n)
    distortions = {index_name: dict() for index_name in indices.keys()}
    for quota in QUOTA_RANGE:
        distortion = {index_name: 0 for index_name in indices.keys()}
        all_powers = {index_name: dict() for index_name in indices.keys()}

        for func, weights in ALL_WVGS_AT_QUOTA[n][quota].items():
            for index_name, index in indices.items():
                all_powers[index_name][func] = index(weights, quota)
            
        for progress in range(trials[n]):
            target = np.random.dirichlet(alpha=[1]*n) # [1-quota-0.001] + [(quota+0.001)/(n-1)]*(n-1) 
            target.sort()

            for index_name in indices.keys():
                best_dist = None
                best_func = None
                for func, powers in all_powers[index_name].items():
                    dist = distance(target, powers, p, EXACT)
                    if best_dist == None or best_dist > dist:
                        best_dist = dist
                        best_func = func

                best_weights = ALL_WVGS_AT_QUOTA[n][quota][best_func]
                veto_indices = [i for i in range(len(best_weights)) if best_weights[i] > 1-quota]

                for i in veto_indices:
                    if (1 - quota) - (target[i]) > 0 and index_name =='no_veto': print(i,target, best_weights, quota)
                    distortion[index_name] = max(distortion[index_name], (1 - quota) - (target[i]))
    
            print_progress(progress,trials[n])

        for index_name in indices.keys():
            distortions[index_name][quota] = distortion[index_name]
    for index_name in indices.keys(): 
        print(distortions[index_name])

    for index_name in indices.keys():
        plt.plot(distortions[index_name].keys(), distortions[index_name].values(), 'o-', label=index_name)
    plt.xlabel('Quota')
    plt.ylabel('Distortion')
    plt.title(f'Veto Distortion for n={n}, trials = {trials[n]}')
    plt.grid(True)
    plt.legend()
    # plt.savefig(f'plots/veto_distortion/veto_distortion_n{n}_nv.png')
    plt.close()
    print(f"Plot saved as veto_distortion_n{n}.png")

        
