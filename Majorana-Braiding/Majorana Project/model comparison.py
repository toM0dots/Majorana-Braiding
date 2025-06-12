import numpy as np
import random
import numpy.linalg as la
import matplotlib.pyplot as plt
import math
from collections import Counter
from typing import List



def chain_builder(N, rho):
    """
    having chain be a numpy array that the first column is the ferminonic number of the pair and the second column 
    is the Majorana operator which can be assigned randomly at contact to simulate quantum randomness.
    third column is the pair index
    """
    if N % 2 != 0:
        raise ValueError("N must be even for a Majorana chain.")
    chain = np.empty((N * rho, 3), dtype=object)
    for i in range(N):
        chain[i * rho][0] = 1 # initially all pairs have fermion number 1
        chain[i * rho][1] = True
        chain[i * rho][2] = int(i/2) # index for pair (i,j), where i<j, using i to denote the Majorana j is paired with.
    return chain

def hopping_annihilate(chain):
    N = len(chain)
     # get the index of Majorana operators
   
    # i =random.randint(0, N-1)
    index = np.where(chain[:,1]==True)[0]
    for i in index:
        
        hope_direction = random.choice([-1, 1])
        if  chain[:, 2].tolist().count(None) == N: #if all sites are empty

            return chain, 0
        if chain[i][0] == None: #if the site is empty
            pass

        elif chain[(i+hope_direction) % N][0] == None: #if the site to arrive is empty
            chain[(i+hope_direction) % N][0] = chain[i][0]
            chain[i][0] = None
            chain[(i+hope_direction) % N][1] = chain[i][1]
            chain[i][1] = None
            chain[(i+hope_direction) % N][2] = chain[i][2]
            chain[i][2] = None
            pass

        elif chain[(i+hope_direction) % N][2] == chain[i][2]:  #if paired under periodic boundary condition
            if chain[i][0] == 0 and chain[(i+hope_direction) % N][0] == 0:
                chain[i][0] = None
                chain[(i+hope_direction) % N][0] = None
                chain[i][1] = None
                chain[(i+hope_direction) % N][1] = None
                chain[i][2] = None
                chain[(i+hope_direction) % N][2] = None
                
        elif chain[i][2] !=  chain[(i+hope_direction) % N][2]:
            rows_i = np.where(chain[:, 2] == chain[i][2])[0]
            rows = np.where(chain[:, 2] == chain[(i+hope_direction) % N][2])[0]
            outside_index_i = [x for x in rows_i if x != i][0]
            outside_index = [x for x in rows if x != (i+hope_direction)%N][0]
        
            if chain[i][0] == 0:
                chain[i][0] = None
                chain[i][1] = None
                chain[i][2] = None
                chain[outside_index_i][0] = None
                chain[outside_index_i][1] = None
                chain[outside_index_i][2] = None
                
            if chain[(i+hope_direction) % N][0] == 0:
                chain[(i+hope_direction) % N][0] = None
                chain[(i+hope_direction) % N][1] = None
                chain[(i+hope_direction) % N][2] = None
                chain[outside_index][0] = None
                chain[outside_index][1] = None
                chain[outside_index][2] = None
                
            if chain[i][0] == 1 and chain[(i+hope_direction) % N][0] == 1:
                # get the index of paired Majorana
                

                # calculate the ferminon number of new arcs
                ferminon_number = random.choice([0, 1])
                n_sum = chain[i][0] + chain[(i+hope_direction) % N][0]
                # ferminon_number = (chain[i][1] + chain[(i+hope_direction) % N][1]) % 2
                
                if ferminon_number == 0:
                    chain[i][0] = None
                    chain[(i+hope_direction) % N][0] = None
                    chain[i][1] = None
                    chain[(i+hope_direction) % N][1] = None
                    chain[i][2] = None
                    chain[(i+hope_direction) % N][2] = None
                    chain[outside_index_i][0] = None
                    chain[outside_index][0] = None
                    chain[outside_index_i][1] = None
                    chain[outside_index][1] = None
                    chain[outside_index_i][2] = None
                    chain[outside_index][2] = None

                else:
                    chain[i][0] = ferminon_number
                    chain[(i+hope_direction) % N][0] = ferminon_number
                    chain[outside_index_i][0] = (n_sum-ferminon_number) % 2 # conservation of ferminon number
                    chain[outside_index][0] = (n_sum-ferminon_number) % 2

                    # update the pairing index
                    
                    pair_index = chain[(i+hope_direction) % N][2]
                    pair_index_i = chain[i][2]
                    min_index = min(pair_index, pair_index_i)
                    max_index = max(pair_index, pair_index_i)

                    if  i == min(i, (i+hope_direction) % N, outside_index_i, outside_index) or (i+hope_direction) % N == min(i, (i+hope_direction) % N, outside_index_i, outside_index):
                        chain[i][2] = min_index
                        chain[(i+hope_direction) % N][2] = min_index
                        chain[outside_index_i][2] = max_index
                        chain[outside_index][2] = max_index
                    else:
                        chain[i][2] = max_index
                        chain[(i+hope_direction) % N][2] = max_index
                        chain[outside_index_i][2] = min_index
                        chain[outside_index][2] = min_index
                
    density = N - chain[:,2].tolist().count(None)  # Count the number of empty sites
    return chain, density
    
    
    

def evolution_annihilate(chain, t):
    """
    Simulates the evolution of the Majorana chain over time t.
    The function modifies the chain in place.
    """
    density = []
    times = [i for i in range(t)]
    for _ in range(t):
        
        chain, x = hopping_annihilate(chain)
        density.append(x)
    return times, density

def hopping_classical(chain):
    N = len(chain)
     # get the index of Majorana operators
   
    # i =random.randint(0, N-1)
    index = np.where(chain[:,1]==True)[0]
    for i in index:
        
        hope_direction = random.choice([-1, 1])
        if  chain[:, 2].tolist().count(None) == N: #if all sites are empty

            return chain, 0
        if chain[i][0] == None: #if the site is empty
            pass

        elif chain[(i+hope_direction) % N][0] == None: #if the site to arrive is empty
            chain[(i+hope_direction) % N][0] = chain[i][0]
            chain[i][0] = None
            chain[(i+hope_direction) % N][1] = chain[i][1]
            chain[i][1] = None
            chain[(i+hope_direction) % N][2] = chain[i][2]
            chain[i][2] = None
            pass

        else :  #if paired under periodic boundary condition:
            rows_i = np.where(chain[:, 2] == chain[i][2])[0]
            rows = np.where(chain[:, 2] == chain[(i+hope_direction) % N][2])[0]
            outside_index_i = [x for x in rows_i if x != i][0]
            outside_index = [x for x in rows if x != (i+hope_direction)%N][0]
            
            chain[i][0] = None
            chain[(i+hope_direction) % N][0] = None
            chain[i][1] = None
            chain[(i+hope_direction) % N][1] = None
            chain[i][2] = None
            chain[(i+hope_direction) % N][2] = None
            chain[outside_index_i][0] = None
            chain[outside_index][0] = None
            chain[outside_index_i][1] = None
            chain[outside_index][1] = None
            chain[outside_index_i][2] = None
            chain[outside_index][2] = None

                
    density = N - chain[:,2].tolist().count(None)  # Count the number of empty sites
    return chain, density
    
    

def evolution_classical(chain, t):
    """
    Simulates the evolution of the Majorana chain over time t.
    The function modifies the chain in place.
    """
    density = []
    times = [i for i in range(t)]
    for _ in range(t):
        
        chain, x = hopping_classical(chain)
        density.append(x)
    return times, density

def main():
    chain1 = chain_builder(10000, 4)
    # print(np.where(chain[:, 1] == True)[0])
    # print(500 - chain[:,2].tolist().count(None))
    times1, density1 = evolution_annihilate(chain1, t =40000)
    chain2 = chain_builder(10000, 4)
    times2, density2 = evolution_classical(chain2, t = 40000)
    tim1 = np.array(times1)
    density_ref1 = 4000/(np.sqrt(2 * np.pi* tim1[1:])) 
    tim2 = np.array(times2)
    density_ref2 = 4000/(np.sqrt(8 * np.pi* tim2[1:])) 
    plt.plot(tim1[1:], density_ref1, label='Quantum Evolution ref', color = 'green')
    plt.plot(tim2[1:], density_ref2, label='Classical Evolution ref', color = 'green')
    plt.plot(times1, density1, label='Quantum Evolution')
    plt.plot(times2, density2, label='Classical Evolution')
    plt.xlabel('Time')
    plt.ylabel('Density')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('density_evolution.png')
    plt.title('Density Evolution of Majorana Chain')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
