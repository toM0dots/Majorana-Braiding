import numpy as np
import random
import numpy.linalg as la
import matplotlib.pyplot as plt
import math
from collections import Counter
from typing import List
from collections import defaultdict



def chain_builder(N, rho):
    """
    having chain be a numpy array that the first column is the ferminonic number of the pair and the second column 
    is the Majorana operator which can be assigned randomly at contact to simulate quantum randomness.
    third column is the pair index
    """
    if N % 2 != 0:
        raise ValueError("N must be even for a Majorana chain.")
    chain = np.full((N* rho, 3), -10) 
    for i in range(N):
        chain[i * rho][0] = 1 # initially all pairs have fermion number 1
        chain[i * rho][1] = 47
        chain[i * rho][2] = int(i/2) # index for pair (i,j), where i<j, using i to denote the Majorana j is paired with.
    return chain

def hopping_annihilate(chain):
    N = len(chain)
     # get the index of Majorana operators
   
    # i =random.randint(0, N-1)
    index = np.where(chain[:,1]==47)[0]
    Nt = len(index)
    for k in range(Nt):
        index = np.where(chain[:,1]==47)[0]
        i = random.choice(index)   
        hope_direction = random.choice([-1, 1])
        if  chain[:, 2].tolist().count(-10) == N: #if all sites are empty

            return chain, 0
        if chain[i][0] == -10: #if the site is empty
            pass

        elif chain[(i+hope_direction) % N][0] == -10: #if the site to arrive is empty
            chain[(i+hope_direction) % N][0] = chain[i][0]
            chain[i][0] = -10
            chain[(i+hope_direction) % N][1] = chain[i][1]
            chain[i][1] = -10
            chain[(i+hope_direction) % N][2] = chain[i][2]
            chain[i][2] = -10
            pass

        elif chain[(i+hope_direction) % N][2] == chain[i][2]:  #if paired under periodic boundary condition
            if chain[i][0] == 0 and chain[(i+hope_direction) % N][0] == 0:
                chain[i][0] = -10
                chain[(i+hope_direction) % N][0] = -10
                chain[i][1] = -10
                chain[(i+hope_direction) % N][1] = -10
                chain[i][2] = -10
                chain[(i+hope_direction) % N][2] = -10
                
        elif chain[i][2] !=  chain[(i+hope_direction) % N][2]:
            rows_i = np.where(chain[:, 2] == chain[i][2])[0]
            rows = np.where(chain[:, 2] == chain[(i+hope_direction) % N][2])[0]
            outside_index_i = [x for x in rows_i if x != i][0]
            outside_index = [x for x in rows if x != (i+hope_direction) % N][0]
        

            # calculate the ferminon number of new arcs
            ferminon_number = random.choice([0, 1])
            n_sum = chain[i][0] + chain[(i+hope_direction) % N][0]
            # ferminon_number = (chain[i][1] + chain[(i+hope_direction) % N][1]) % 2

            pair_index = chain[(i+hope_direction) % N][2]
            pair_index_i = chain[i][2]
            min_index = min(pair_index, pair_index_i)
            max_index = max(pair_index, pair_index_i)
            if ferminon_number == 0:
                chain[i][0] = -10
                chain[(i+hope_direction) % N][0] = -10
                chain[i][1] = -10
                chain[(i+hope_direction) % N][1] = -10
                chain[i][2] = -10
                chain[(i+hope_direction) % N][2] = -10

                chain[outside_index_i][0] = (n_sum-ferminon_number) % 2 # conservation of ferminon number
                chain[outside_index][0] = (n_sum-ferminon_number) % 2  
                chain[outside_index_i][2] = min_index
                chain[outside_index][2] = min_index

            else:
                chain[i][0] = ferminon_number
                chain[(i+hope_direction) % N][0] = ferminon_number
                chain[outside_index_i][0] = (n_sum-ferminon_number) % 2 # conservation of ferminon number
                chain[outside_index][0] = (n_sum-ferminon_number) % 2

                # update the pairing index
                
               

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
                
    density = N - chain[:,2].tolist().count(-10)  # Count the number of empty sites
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

def hopping_classical(chain, p):
    N = len(chain)
     # get the index of Majorana operators
   
    # i =random.randint(0, N-1)
    index = np.where(chain[:,1]==47)[0]
    Nt = len(index)
    for _ in range(Nt):
        index = np.where(chain[:,1]==47)[0]
        i = random.choice(index)
        hope_direction = random.choice([-1, 1])
        if  chain[:, 2].tolist().count(-10) == N: #if all sites are empty

            return chain, 0
        if chain[i][0] == -10: #if the site is empty
            pass

        elif chain[(i+hope_direction) % N][0] == -10: #if the site to arrive is empty
            chain[(i+hope_direction) % N][0] = chain[i][0]
            chain[i][0] = -10
            chain[(i+hope_direction) % N][1] = chain[i][1]
            chain[i][1] = -10
            chain[(i+hope_direction) % N][2] = chain[i][2]
            chain[i][2] = -10
            pass

        else:
            r = random.random()
            if r <p :
                chain[i][0] = -10
                chain[(i+hope_direction) % N][0] = -10
                chain[i][1] = -10
                chain[(i+hope_direction) % N][1] = -10
                chain[i][2] = -10
                chain[(i+hope_direction) % N][2] = -10
            else:
                pass
        # elif chain[i][2] !=  chain[(i+hope_direction) % N][2]:
        #   #if paired under periodic boundary condition:
        #     rows_i = np.where(chain[:, 2] == chain[i][2])[0]
        #     rows = np.where(chain[:, 2] == chain[(i+hope_direction) % N][2])[0]
        #     outside_index_i = [x for x in rows_i if x != i][0]
        #     outside_index = [x for x in rows if x != (i+hope_direction)%N][0]
            
        #     chain[i][0] = -10
        #     chain[(i+hope_direction) % N][0] = -10
        #     chain[i][1] = -10
        #     chain[(i+hope_direction) % N][1] = -10
        #     chain[i][2] = -10
        #     chain[(i+hope_direction) % N][2] = -10
        #     chain[outside_index_i][0] = -10
        #     chain[outside_index][0] = -10
        #     chain[outside_index_i][1] = -10
        #     chain[outside_index][1] = -10
        #     chain[outside_index_i][2] = -10
        #     chain[outside_index][2] = -10

                
    density = N - chain[:,2].tolist().count(-10)  # Count the number of empty sites
    return chain, density
    
    

def evolution_classical(chain, t, p =1):
    """
    Simulates the evolution of the Majorana chain over time t.
    The function modifies the chain in place.
    """
    density = []
    times = [i for i in range(t)]
    for _ in range(t):
        
        chain, x = hopping_classical(chain, p)
        density.append(x)
    return times, density

def hopping(chain_):
    chain = chain_.copy()
    N = len(chain)
    index = np.where(chain[:,1] == 47)[0] 
    Nt = len(index)  # Number of Majoranas present
    # i =random.randint(0, N-1)
    for _ in range(Nt):
        index = np.where(chain[:, 1] == 47)[0]  # Update indices of Majoranas that are present
        i = random.choice(index)
        hope_direction = np.random.choice([-1, 1])
        if chain[i][0] == -10: #if the site is empty
            pass

        elif chain[(i+hope_direction) % N][0] == -10: #if the site to arrive is empty
            chain[(i+hope_direction) % N][0] = chain[i][0]
            chain[i][0] = -10
            chain[(i+hope_direction) % N][1] = chain[i][1]
            chain[i][1] = -10
            chain[(i+hope_direction) % N][2] = chain[i][2]
            chain[i][2] = -10
            chain[(i+hope_direction) % N][3] = chain[i][3]
            chain[i][3] = -10
            

        elif chain[(i+hope_direction) % N][2] == chain[i][2]:  #if paired under periodic boundary condition
            pass # do nothing, the Majorana operators are already paired
        
        else:
            #if hope_direction < 0:
            # get the index of paired Majorana
            rows_i = np.where(chain[:, 2] == chain[i][2])[0]
            rows = np.where(chain[:, 2] == chain[(i+hope_direction) % N][2])[0]
            outside_index_i = [x for x in rows_i if x != i][0]
            outside_index = [x for x in rows if x != (i+hope_direction)%N][0]

            # calculate the ferminon number of new arcs
            ferminon_number =  np.random.choice([0, 1])
            n_sum = chain[i][0] + chain[(i+hope_direction) % N][0]

            
            
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
    return chain


def evolution(chain, t, plot=False):
    """
    Simulates the evolution of the Majorana chain over time t.
    The function modifies the chain in place.
    """
    if plot: 
        N = len(np.where(chain[:, 1] == 47)[0])  # Number of Majoranas
        history = np.empty((t, N), dtype=object)  # Initialize history for plotting
        for i in range(t):
            chain = hopping(chain)
            wl = worldline(chain)
            history[i,:] = wl.reshape(-1)
        return chain, history
    else:
        for _ in range(t):
            chain = hopping(chain)
        return chain
    
def worldline(chain):
    """
    Generates a worldline representation of the Majorana chain.
    The function returns a list of pairs (Majorana index, fermion number).
    """
    index = np.where(chain[:, 1] == 47)[0]  # Get indices of Majoranas that are present
    wl = np.zeros((1, len(index)), dtype=object)  # Initialize an empty list for the worldline
    for i in index:
        wl[0][chain[i][3]] = i  # position of ith particle
    return wl


def chain_builder_2L(N, rho):
    """
    having chain be a numpy array that the first column is the ferminonic number of the pair and the second column 
    is the Majorana operator which can be assigned randomly at contact to simulate quantum randomness.
    third column is the pair index
    """
    if N % 2 != 0:
        raise ValueError("N must be even for a Majorana chain.")
    chain = np.full(( N*rho, 3, 2 ), -10) 
    for i in range(N):
        chain[(i * rho) % (N * rho)][0][0] = 0 # initially all pairs have fermion number 1
        chain[(i * rho) % (N * rho)][1][0] = 47
        chain[(i * rho) % (N * rho)][2][0] = int(i/2) # index for pair (i,j), where i<j, using i to denote the Majorana j is paired with.
        chain[(i * rho ) % (N * rho)][0][1] = 1 # initially all pairs have fermion number 1
        chain[(i * rho ) % (N * rho)][1][1] = 47
        chain[(i * rho ) % (N * rho)][2][1] = int(i/2) # index for pair (i,j), where i<j, using i to denote the Majorana j is paired with.
    return chain
    
def hopping_annihilate_2L(chain, p):
    N = chain.shape[0]  # Number of sites in the chain

    index0 = np.where(chain[:,1, 0]==47)[0] # get the index of Majorana operators
    index1 = np.where(chain[:,1, 1]==47)[0]
    Nt0 = len(index0)
    Nt1 = len(index1)
    Nt = Nt0 + Nt1
    for _ in range(Nt):
         if  chain[:, 2, 0].tolist().count(-10) == N or chain[:, 2, 1].tolist().count(-10) == N: #if all sites are empty
    
             return chain, 0, 0, 0 ,0
             
         if random.random() < Nt0/ (Nt0 + Nt1):  
            index0 = np.where(chain[:,0, 0]!=-10)[0]
            i = random.choice(index0)
            hop_direction0 = random.choice([-1, 1])
            
            
            if chain[i][0][0] == -10: #if the site is empty
                pass
    
            elif chain[(i+hop_direction0) % N][0][0] == -10: #if the site to arrive is empty
                
                r = random.random()
                if r < p:
                    if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][1] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[(i+hop_direction0) % N][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                        
                        if partner_j == -10:
                            raise ValueError("Partner Majorana is -10, check the chain structure.")
                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[(i+hop_direction0) % N][0][1] = (chain[(i+hop_direction0) % N][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                    elif hop_direction0 >0 and chain[i][2][1] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[i][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != i][0]
    
                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[i][0][1] = (chain[i][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                chain[(i+hop_direction0) % N][0][0] = chain[i][0][0]
                chain[i][0][0] = -10
                chain[(i+hop_direction0) % N][1][0] = chain[i][1][0]
                chain[i][1][0] = -10
                chain[(i+hop_direction0) % N][2][0] = chain[i][2][0]
                chain[i][2][0] = -10
    
            elif chain[(i+hop_direction0) % N][2][0] == chain[i][2][0]:  #if paired under periodic boundary condition
    
                r = random.random()
                if r < p:
                    if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][1] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[(i+hop_direction0) % N][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
    
                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[(i+hop_direction0) % N][0][1] = (chain[(i+hop_direction0) % N][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                    elif hop_direction0 >0 and chain[i][2][1] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[i][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != i][0]
    
                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[i][0][1] = (chain[i][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                if chain[i][0][0] == 0 and chain[(i+hop_direction0) % N][0][0] == 0:
                    
                    chain[i][0][0] = -10
                    chain[(i+hop_direction0) % N][0][0] = -10
                    chain[i][1][0] = -10
                    chain[(i+hop_direction0) % N][1][0] = -10
                    chain[i][2][0] = -10
                    chain[(i+hop_direction0) % N][2][0] = -10
                else:
                    r = random.random()
                    if r < p:
                        if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][1] != -10: #if the site is empty
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                            pair_j = np.where(chain[:, 2, 1] == chain[(i+hop_direction0) % N][2][1])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
    
                            # flip the fermionic parity
                            chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                            chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                            chain[(i+hop_direction0) % N][0][1] = (chain[(i+hop_direction0) % N][0][1]+ 1) % 2
                            chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                        elif hop_direction0 >0 and chain[i][2][1] != -10:
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                            pair_j = np.where(chain[:, 2, 1] == chain[i][2][1])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != i][0]
    
                            # flip the fermionic parity
                            chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                            chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                            chain[i][0][1] = (chain[i][0][1]+ 1) % 2
                            chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                
            
            elif chain[i][2][0] !=  chain[(i+hop_direction0) % N][2][0] :
    
    
                # get the index of paired Majorana
                rows_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
    
                if len(rows_i) == 0:
                    raise ValueError("No paired Majorana found for site i.")
                
                rows = np.where(chain[:, 2, 0] == chain[(i+hop_direction0) % N][2][0])[0]
                outside_index_i = [x for x in rows_i if x != i][0]
                outside_index = [x for x in rows if x != (i+hop_direction0)%N][0]
                
                # 1st braiding
                r = random.random()
                if r < p:
                    if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][1] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[(i+hop_direction0) % N][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
    
                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i ][0][0] = (chain[partner_i ][0][0] + 1) % 2
                        chain[(i+hop_direction0) % N][0][1] = (chain[(i+hop_direction0) % N][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                    elif hop_direction0 >0 and chain[i][2][1] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[i][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != i][0]
    
                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i ][0][0] = (chain[partner_i ][0][0] + 1) % 2
                        chain[i][0][1] = (chain[i][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                pair_index = chain[(i+hop_direction0) % N][2][0]
                pair_index_i = chain[i][2][0]
                min_index = min(pair_index, pair_index_i)
                max_index = max(pair_index, pair_index_i)

                # repairing at collision after 1st braid and calculate the sum of parity
                n_sum = chain[i][0][0] + chain[(i+hop_direction0) % N][0][0]
                ferminon_number = random.choice([0, 1])
                chain[i][0][0] = ferminon_number
                chain[(i+hop_direction0) % N][0][0] = ferminon_number
    
                if chain[i][0][0] == 0 and chain[(i+hop_direction0) % N][0][0] == 0:
                    # annihilate the collided pair
                    chain[i][0][0] = -10
                    chain[i][1][0] = -10
                    chain[i][2][0] = -10
                    chain[(i+hop_direction0) % N][0][0] = -10
                    chain[(i+hop_direction0) % N][1][0] = -10
                    chain[(i+hop_direction0) % N][2][0] = -10
    
                    # update the partner pair
                    chain[outside_index_i][0][0] = n_sum % 2
                    chain[outside_index_i][1][0] = 47
                    chain[outside_index_i][2][0] = min_index
                    chain[outside_index][0][0] = n_sum % 2
                    chain[outside_index][1][0] = 47
                    chain[outside_index][2][0] = min_index
    
            
                # bound back and 2nd braid
                elif chain[i][0][0] == 1 and chain[(i+hop_direction0) % N][0][0] == 1:
                    r = random.random()
                    if r < p:
                        if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][1] != -10: #if the site is empty
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                            pair_j = np.where(chain[:, 2, 1] == chain[(i+hop_direction0) % N][2][1])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
    
                            # flip the fermionic parity
                            chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                            chain[partner_i ][0][0] = (chain[partner_i ][0][0] + 1) % 2
                            chain[(i+hop_direction0) % N][0][1] = (chain[(i+hop_direction0) % N][0][1]+ 1) % 2
                            chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                        elif hop_direction0 >0 and chain[i][2][1] != -10:
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                            pair_j = np.where(chain[:, 2, 1] == chain[i][2][1])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != i][0]
    
                            # flip the fermionic parity
                            chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                            chain[partner_i ][0][0] = (chain[partner_i ][0][0] + 1) % 2
                            chain[i][0][1] = (chain[i][0][1]+ 1) % 2
                            chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
                    # calculate the ferminon number of new arcs
                    chain[outside_index_i][0][0] = (n_sum - chain[i][0][0]) % 2
                    chain[outside_index][0][0] = (n_sum -  chain[(i+hop_direction0) % N][0][0]) % 2
                    
    
                    if  i == min(i, (i+hop_direction0) % N, outside_index_i, outside_index) or (i+hop_direction0) % N == min(i, (i+hop_direction0) % N, outside_index_i, outside_index):
                        chain[i][2][0] = min_index
                        chain[(i+hop_direction0) % N][2][0] = min_index
                        chain[outside_index_i][2][0] = max_index
                        chain[outside_index][2][0] = max_index
                    else:
                        chain[i][2][0] = max_index
                        chain[(i+hop_direction0) % N][2][0] = max_index
                        chain[outside_index_i][2][0] = min_index
                        chain[outside_index][2][0] = min_index
    
                    
         else:
            index1 = np.where(chain[:,0, 1]!= -10)[0] # get the index of Majorana operators
            j = random.choice(index1)
            hop_direction1 = random.choice([-1, 1])
            if  chain[:, 2, 1].tolist().count(-10) == N: #if all sites are empty
    
                return chain, 0, 0
            
            if chain[j][0][1] == -10: #if the site is empty
                pass
    
            elif chain[(j+hop_direction1) % N][0][1] == -10: #if the site to arrive is empty
    
                q = random.random()
                if q < p:
                    if hop_direction1 >0 and chain[(j+hop_direction1) % N ][2][0] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        pair_i = np.where(chain[:, 2, 0] == chain[(j+hop_direction1) % N][2][0])[0]
                        partner_j = [x for x in pair_j if x != j][0]
                        partner_i = [x for x in pair_i if x != (j+hop_direction1)%N][0]
    
                        # flip the fermionic parity
                        chain[(j+hop_direction1) % N][0][0] = (chain[(j+hop_direction1) % N][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                    elif hop_direction1 <0 and chain[j][2][0] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[j][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        partner_i = [x for x in pair_i if x != j][0]
                        partner_j = [x for x in pair_j if x != j][0]
    
                        # flip the fermionic parity
                        chain[j][0][0] = (chain[j][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
    
                chain[(j+hop_direction1) % N][0][1] = chain[j][0][1]
                chain[j][0][1] = -10
                chain[(j+hop_direction1) % N][1][1] = chain[j][1][1]
                chain[j][1][1] = -10
                chain[(j+hop_direction1) % N][2][1] = chain[j][2][1]
                chain[j][2][1] = -10
                
            elif chain[(j+hop_direction1) % N][2][1] == chain[j][2][1]:  #if paired under periodic boundary condition
    
                q = random.random()
                if q < p:
                    if hop_direction1 >0 and chain[(j+hop_direction1) % N ][2][0] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        pair_i = np.where(chain[:, 2, 0] == chain[(j+hop_direction1) % N][2][0])[0]
                        partner_j = [x for x in pair_j if x != j][0]
                        partner_i = [x for x in pair_i if x != (j+hop_direction1)%N][0]
    
                        # flip the fermionic parity
                        chain[(j+hop_direction1) % N][0][0] = (chain[(j+hop_direction1) % N][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2
    
                    elif hop_direction1 <0 and chain[j][2][0] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[j][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        partner_i = [x for x in pair_i if x != j][0]
                        partner_j = [x for x in pair_j if x != j][0]
    
                        # flip the fermionic parity
                        chain[j][0][0] = (chain[j][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2
    
                if chain[j][0][1] == 0 and chain[(j+hop_direction1) % N][0][1] == 0:       
                    chain[j][0][1] = -10
                    chain[(j+hop_direction1) % N][0][1] = -10
                    chain[j][1][1] = -10
                    chain[(j+hop_direction1) % N][1][1] = -10
                    chain[j][2][1] = -10
                    chain[(j+hop_direction1) % N][2][1] = -10
    
                else:
                    q = random.random()
                    if q < p:
                        if hop_direction1 >0 and chain[(j+hop_direction1) % N ][2][0] != -10: #if the site is empty
                            # get the index of paired Majorana
                            pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                            pair_i = np.where(chain[:, 2, 0] == chain[(j+hop_direction1) % N][2][0])[0]
                            partner_j = [x for x in pair_j if x != j][0]
                            partner_i = [x for x in pair_i if x != (j+hop_direction1)%N][0]
    
                            # flip the fermionic parity
                            chain[(j+hop_direction1) % N][0][0] = (chain[(j+hop_direction1) % N][0][0] +1 ) % 2
                            chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                            chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                            chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2
    
                        elif hop_direction1 <0 and chain[j][2][0] != -10:
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, 0] == chain[j][2][0])[0]
                            pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                            partner_i = [x for x in pair_i if x != j][0]
                            partner_j = [x for x in pair_j if x != j][0]
    
                            # flip the fermionic parity
                            chain[j][0][0] = (chain[j][0][0] +1 ) % 2
                            chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                            chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                            chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2
            
            elif chain[j][2][1] !=  chain[(j+hop_direction1) % N][2][1]:
    
            
                # get the index of paired Majorana
                rows_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                rows = np.where(chain[:, 2, 1] == chain[(j+hop_direction1) % N][2][1])[0]
                outside_index_j = [x for x in rows_j if x != j][0]
                outside_index = [x for x in rows if x != (j+hop_direction1)%N][0]
                
    
                # 1st braiding 
                q = random.random()
                if q < p:
                    if hop_direction1 >0 and chain[(j+hop_direction1) % N ][2][0] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        pair_i = np.where(chain[:, 2, 0] == chain[(j+hop_direction1) % N][2][0])[0]
                        partner_j = [x for x in pair_j if x != j][0]
                        partner_i = [x for x in pair_i if x != (j+hop_direction1)%N][0]
    
                        # flip the fermionic parity
                        chain[(j+hop_direction1) % N][0][0] = (chain[(j+hop_direction1) % N][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2
    
                    elif hop_direction1 <0 and chain[j][2][0] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[j][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        partner_i = [x for x in pair_i if x != j][0]
                        partner_j = [x for x in pair_j if x != j][0]
    
                        # flip the fermionic parity
                        chain[j][0][0] = (chain[j][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2
    
                pair_index1 = chain[(j+hop_direction1) % N][2][1]
                pair_index_j = chain[j][2][1]
                min_index = min(pair_index1, pair_index_j)
                max_index = max(pair_index1, pair_index_j)
                # update the repairing index 

                
                
                # repairing, and calculate the parity sum after braiding 
                n_sum1 = chain[j][0][1] + chain[(j+hop_direction1) % N][0][1]
                ferminon_number1 = random.choice([0, 1])
                chain[j][0][1] = ferminon_number1
                chain[(j+hop_direction1) % N][0][1] = ferminon_number1
    
                if chain[j][0][1] == 0 and chain[(j+hop_direction1) % N][0][1] == 0:
                    chain[j][0][1] = -10
                    chain[j][1][1] = -10
                    chain[j][2][1] = -10
                    chain[(j+hop_direction1) % N][0][1] = -10
                    chain[(j+hop_direction1) % N][1][1] = -10
                    chain[(j+hop_direction1) % N][2][1] = -10
    
    
                    chain[outside_index_j][0][1] = n_sum1 % 2
                    chain[outside_index_j][1][1] = 47
                    chain[outside_index_j][2][1] = min_index
                    chain[outside_index][0][1] = n_sum1 % 2
                    chain[outside_index][1][1] = 47
                    chain[outside_index][2][1] = min_index
    
                # bound back and 2nd braid
                elif chain[j][0][1] == 1 and chain[(j+hop_direction1) % N][0][1] == 1:
                    q = random.random()
                    if q < p:
                        if hop_direction1 >0 and chain[(j+hop_direction1) % N ][2][0] != -10: #if the site is empty
                            # get the index of paired Majorana
                            pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                            pair_i = np.where(chain[:, 2, 0] == chain[(j+hop_direction1) % N][2][0])[0]
                            partner_j = [x for x in pair_j if x != j][0]
                            partner_i = [x for x in pair_i if x != (j+hop_direction1)%N][0]
    
                            # flip the fermionic parity
                            chain[(j+hop_direction1) % N][0][0] = (chain[(j+hop_direction1) % N][0][0] +1 ) % 2
                            chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                            chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                            chain[(j+hop_direction1) % N][0][1] = (chain[(j+hop_direction1) % N][0][1]+ 1) % 2
    
                        elif hop_direction1 <0 and chain[j][2][0] != -10:
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, 0] == chain[j][2][0])[0]
                            pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                            partner_i = [x for x in pair_i if x != j][0]
                            partner_j = [x for x in pair_j if x != j][0]
    
                            # flip the fermionic parity
                            chain[j][0][0] = (chain[j][0][0] +1 ) % 2
                            chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                            chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                            chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2
                    
                    #update the parity
                    chain[outside_index_j][0][1] = (n_sum1-chain[j][0][1]) % 2 # conservation of ferminon number
                    chain[outside_index][0][1] = (n_sum1-chain[(j+hop_direction1) % N][0][1]) % 2
    
                    # update the pairing index
                    
                    if  j == min(j, (j+hop_direction1) % N, outside_index_j, outside_index) or (j+hop_direction1) % N == min(j, (j+hop_direction1) % N, outside_index_j, outside_index):
                        chain[j][2][1] = min_index
                        chain[(j+hop_direction1) % N][2][1] = min_index
                        chain[outside_index_j][2][1] = max_index
                        chain[outside_index][2][1] = max_index
                    else:
                        chain[j][2][1] = max_index
                        chain[(j+hop_direction1) % N][2][1] = max_index
                        chain[outside_index_j][2][1] = min_index
                        chain[outside_index][2][1] = min_index

    
    
    density0 = N - chain[:,2,0].tolist().count(-10) # Count the number of empty sites
    density1 = N - chain[:,2,1].tolist().count(-10) # Count the number of empty sites
    parity0 = chain[:, 0, 0].sum()
    parity1 = chain[:, 0 ,1].sum()
    # density = [density0, density1]
    return chain, density0 , density1, parity0, parity1


def evolution_annihilate_2L(chain, t, p):
    """
    Simulates the evolution of the Majorana chain over time t.
    The function modifies the chain in place.
    """
    density = np.zeros((t, 2), dtype=int)  # Initialize density array for two Majorana operators
    parity = np.zeros((t, 2), dtype =int)
    times = [k for k in range(t)]
    for k in range(t):
        
        chain, density0, density1, parity0, parity1= hopping_annihilate_2L(chain, p)
        density[k][0] = density0
        density[k][1] = density1
        parity[k][0] = parity0
        parity[k][1] = parity1
    return times, density, parity


def hopping_no_annihilate_2L(chain_, p):
    chain = chain_.copy()
    N = chain.shape[0]  # Number of sites in the chain

    index0 = np.where(chain[:,1, 0]==47)[0] # get the index of Majorana operators
    index1 = np.where(chain[:,1, 1]==47)[0]
    Nt0 = len(index0)
    Nt1 = len(index1)
    Nt = Nt0 + Nt1
    for _ in range(Nt):
         if  chain[:, 2, 0].tolist().count(-10) == N or chain[:, 2, 1].tolist().count(-10) == N: #if all sites are empty

             return chain, 0, 0, 0 ,0

         if random.random() < Nt0/ (Nt0 + Nt1):
            index0 = np.where(chain[:,0, 0]!=-10)[0]
            i = random.choice(index0)
            hop_direction0 = random.choice([-1, 1])


            if chain[i][0][0] == -10: #if the site is empty
                pass

            elif chain[(i+hop_direction0) % N][0][0] == -10: #if the site to arrive is empty

                r = random.random()
                if r < p:
                    if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][1] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[(i+hop_direction0) % N][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]

                        if partner_j == -10:
                            raise ValueError("Partner Majorana is -10, check the chain structure.")
                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[(i+hop_direction0) % N][0][1] = (chain[(i+hop_direction0) % N][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2

                    elif hop_direction0 >0 and chain[i][2][1] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[i][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != i][0]

                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[i][0][1] = (chain[i][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2

                chain[(i+hop_direction0) % N][0][0] = chain[i][0][0]
                chain[i][0][0] = -10
                chain[(i+hop_direction0) % N][1][0] = chain[i][1][0]
                chain[i][1][0] = -10
                chain[(i+hop_direction0) % N][2][0] = chain[i][2][0]
                chain[i][2][0] = -10

            elif chain[(i+hop_direction0) % N][2][0] == chain[i][2][0]:  #if paired under periodic boundary condition

                r = random.random()
                if r < p:
                    if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][1] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[(i+hop_direction0) % N][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]

                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[(i+hop_direction0) % N][0][1] = (chain[(i+hop_direction0) % N][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2

                    elif hop_direction0 >0 and chain[i][2][1] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[i][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != i][0]

                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[i][0][1] = (chain[i][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2


                r = random.random()
                if r < p:
                    if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][1] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[(i+hop_direction0) % N][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]

                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[(i+hop_direction0) % N][0][1] = (chain[(i+hop_direction0) % N][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2

                    elif hop_direction0 >0 and chain[i][2][1] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[i][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != i][0]

                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[i][0][1] = (chain[i][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2



            elif chain[i][2][0] !=  chain[(i+hop_direction0) % N][2][0] :


                # get the index of paired Majorana
                rows_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]

                if len(rows_i) == 0:
                    raise ValueError("No paired Majorana found for site i.")

                rows = np.where(chain[:, 2, 0] == chain[(i+hop_direction0) % N][2][0])[0]
                outside_index_i = [x for x in rows_i if x != i][0]
                outside_index = [x for x in rows if x != (i+hop_direction0)%N][0]

                # 1st braiding
                r = random.random()
                if r < p:
                    if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][1] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[(i+hop_direction0) % N][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]

                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i ][0][0] = (chain[partner_i ][0][0] + 1) % 2
                        chain[(i+hop_direction0) % N][0][1] = (chain[(i+hop_direction0) % N][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2

                    elif hop_direction0 >0 and chain[i][2][1] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[i][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != i][0]

                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i ][0][0] = (chain[partner_i ][0][0] + 1) % 2
                        chain[i][0][1] = (chain[i][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2

                pair_index = chain[(i+hop_direction0) % N][2][0]
                pair_index_i = chain[i][2][0]
                min_index = min(pair_index, pair_index_i)
                max_index = max(pair_index, pair_index_i)

                # repairing at collision after 1st braid and calculate the sum of parity
                n_sum = chain[i][0][0] + chain[(i+hop_direction0) % N][0][0]
                ferminon_number = random.choice([0, 1])
                chain[i][0][0] = ferminon_number
                chain[(i+hop_direction0) % N][0][0] = ferminon_number


                r = random.random()
                if r < p:
                    if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][1] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[(i+hop_direction0) % N][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]

                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i ][0][0] = (chain[partner_i ][0][0] + 1) % 2
                        chain[(i+hop_direction0) % N][0][1] = (chain[(i+hop_direction0) % N][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2

                    elif hop_direction0 >0 and chain[i][2][1] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[i][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[i][2][1])[0]
                        partner_i = [x for x in pair_i if x != i][0]
                        partner_j = [x for x in pair_j if x != i][0]

                        # flip the fermionic parity
                        chain[i][0][0] = (chain[i][0][0] +1 ) % 2
                        chain[partner_i ][0][0] = (chain[partner_i ][0][0] + 1) % 2
                        chain[i][0][1] = (chain[i][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2
                # calculate the ferminon number of new arcs
                chain[outside_index_i][0][0] = (n_sum - chain[i][0][0]) % 2
                chain[outside_index][0][0] = (n_sum -  chain[(i+hop_direction0) % N][0][0]) % 2


                if  i == min(i, (i+hop_direction0) % N, outside_index_i, outside_index) or (i+hop_direction0) % N == min(i, (i+hop_direction0) % N, outside_index_i, outside_index):
                    chain[i][2][0] = min_index
                    chain[(i+hop_direction0) % N][2][0] = min_index
                    chain[outside_index_i][2][0] = max_index
                    chain[outside_index][2][0] = max_index
                else:
                    chain[i][2][0] = max_index
                    chain[(i+hop_direction0) % N][2][0] = max_index
                    chain[outside_index_i][2][0] = min_index
                    chain[outside_index][2][0] = min_index


         else:
            index1 = np.where(chain[:,0, 1]!= -10)[0] # get the index of Majorana operators
            j = random.choice(index1)
            hop_direction1 = random.choice([-1, 1])
            if  chain[:, 2, 1].tolist().count(-10) == N: #if all sites are empty

                return chain, 0, 0

            if chain[j][0][1] == -10: #if the site is empty
                pass

            elif chain[(j+hop_direction1) % N][0][1] == -10: #if the site to arrive is empty

                q = random.random()
                if q < p:
                    if hop_direction1 >0 and chain[(j+hop_direction1) % N ][2][0] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        pair_i = np.where(chain[:, 2, 0] == chain[(j+hop_direction1) % N][2][0])[0]
                        partner_j = [x for x in pair_j if x != j][0]
                        partner_i = [x for x in pair_i if x != (j+hop_direction1)%N][0]

                        # flip the fermionic parity
                        chain[(j+hop_direction1) % N][0][0] = (chain[(j+hop_direction1) % N][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2

                    elif hop_direction1 <0 and chain[j][2][0] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[j][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        partner_i = [x for x in pair_i if x != j][0]
                        partner_j = [x for x in pair_j if x != j][0]

                        # flip the fermionic parity
                        chain[j][0][0] = (chain[j][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[partner_j][0][1] = (chain[partner_j][0][1]+ 1) % 2

                chain[(j+hop_direction1) % N][0][1] = chain[j][0][1]
                chain[j][0][1] = -10
                chain[(j+hop_direction1) % N][1][1] = chain[j][1][1]
                chain[j][1][1] = -10
                chain[(j+hop_direction1) % N][2][1] = chain[j][2][1]
                chain[j][2][1] = -10

            elif chain[(j+hop_direction1) % N][2][1] == chain[j][2][1]:  #if paired under periodic boundary condition

                q = random.random()
                if q < p:
                    if hop_direction1 >0 and chain[(j+hop_direction1) % N ][2][0] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        pair_i = np.where(chain[:, 2, 0] == chain[(j+hop_direction1) % N][2][0])[0]
                        partner_j = [x for x in pair_j if x != j][0]
                        partner_i = [x for x in pair_i if x != (j+hop_direction1)%N][0]

                        # flip the fermionic parity
                        chain[(j+hop_direction1) % N][0][0] = (chain[(j+hop_direction1) % N][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2

                    elif hop_direction1 <0 and chain[j][2][0] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[j][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        partner_i = [x for x in pair_i if x != j][0]
                        partner_j = [x for x in pair_j if x != j][0]

                        # flip the fermionic parity
                        chain[j][0][0] = (chain[j][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2


                q = random.random()
                if q < p:
                    if hop_direction1 >0 and chain[(j+hop_direction1) % N ][2][0] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        pair_i = np.where(chain[:, 2, 0] == chain[(j+hop_direction1) % N][2][0])[0]
                        partner_j = [x for x in pair_j if x != j][0]
                        partner_i = [x for x in pair_i if x != (j+hop_direction1)%N][0]

                        # flip the fermionic parity
                        chain[(j+hop_direction1) % N][0][0] = (chain[(j+hop_direction1) % N][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2

                    elif hop_direction1 <0 and chain[j][2][0] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[j][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        partner_i = [x for x in pair_i if x != j][0]
                        partner_j = [x for x in pair_j if x != j][0]

                        # flip the fermionic parity
                        chain[j][0][0] = (chain[j][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2

            elif chain[j][2][1] !=  chain[(j+hop_direction1) % N][2][1]:


                # get the index of paired Majorana
                rows_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                rows = np.where(chain[:, 2, 1] == chain[(j+hop_direction1) % N][2][1])[0]
                outside_index_j = [x for x in rows_j if x != j][0]
                outside_index = [x for x in rows if x != (j+hop_direction1)%N][0]


                # 1st braiding
                q = random.random()
                if q < p:
                    if hop_direction1 >0 and chain[(j+hop_direction1) % N ][2][0] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        pair_i = np.where(chain[:, 2, 0] == chain[(j+hop_direction1) % N][2][0])[0]
                        partner_j = [x for x in pair_j if x != j][0]
                        partner_i = [x for x in pair_i if x != (j+hop_direction1)%N][0]

                        # flip the fermionic parity
                        chain[(j+hop_direction1) % N][0][0] = (chain[(j+hop_direction1) % N][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2

                    elif hop_direction1 <0 and chain[j][2][0] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[j][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        partner_i = [x for x in pair_i if x != j][0]
                        partner_j = [x for x in pair_j if x != j][0]

                        # flip the fermionic parity
                        chain[j][0][0] = (chain[j][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2

                pair_index1 = chain[(j+hop_direction1) % N][2][1]
                pair_index_j = chain[j][2][1]
                min_index = min(pair_index1, pair_index_j)
                max_index = max(pair_index1, pair_index_j)
                # update the repairing index



                # repairing, and calculate the parity sum after braiding
                n_sum1 = chain[j][0][1] + chain[(j+hop_direction1) % N][0][1]
                ferminon_number1 = random.choice([0, 1])
                chain[j][0][1] = ferminon_number1
                chain[(j+hop_direction1) % N][0][1] = ferminon_number1

                q = random.random()
                if q < p:
                    if hop_direction1 >0 and chain[(j+hop_direction1) % N ][2][0] != -10: #if the site is empty
                        # get the index of paired Majorana
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        pair_i = np.where(chain[:, 2, 0] == chain[(j+hop_direction1) % N][2][0])[0]
                        partner_j = [x for x in pair_j if x != j][0]
                        partner_i = [x for x in pair_i if x != (j+hop_direction1)%N][0]

                        # flip the fermionic parity
                        chain[(j+hop_direction1) % N][0][0] = (chain[(j+hop_direction1) % N][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[(j+hop_direction1) % N][0][1] = (chain[(j+hop_direction1) % N][0][1]+ 1) % 2

                    elif hop_direction1 <0 and chain[j][2][0] != -10:
                        # get the index of paired Majorana
                        pair_i = np.where(chain[:, 2, 0] == chain[j][2][0])[0]
                        pair_j = np.where(chain[:, 2, 1] == chain[j][2][1])[0]
                        partner_i = [x for x in pair_i if x != j][0]
                        partner_j = [x for x in pair_j if x != j][0]

                        # flip the fermionic parity
                        chain[j][0][0] = (chain[j][0][0] +1 ) % 2
                        chain[partner_i][0][0] = (chain[partner_i][0][0] + 1) % 2
                        chain[j][0][1] = (chain[j][0][1]+ 1) % 2
                        chain[ partner_j][0][1] = (chain[ partner_j][0][1]+ 1) % 2

                    #update the parity
                    chain[outside_index_j][0][1] = (n_sum1-chain[j][0][1]) % 2 # conservation of ferminon number
                    chain[outside_index][0][1] = (n_sum1-chain[(j+hop_direction1) % N][0][1]) % 2

                    # update the pairing index

                    if  j == min(j, (j+hop_direction1) % N, outside_index_j, outside_index) or (j+hop_direction1) % N == min(j, (j+hop_direction1) % N, outside_index_j, outside_index):
                        chain[j][2][1] = min_index
                        chain[(j+hop_direction1) % N][2][1] = min_index
                        chain[outside_index_j][2][1] = max_index
                        chain[outside_index][2][1] = max_index
                    else:
                        chain[j][2][1] = max_index
                        chain[(j+hop_direction1) % N][2][1] = max_index
                        chain[outside_index_j][2][1] = min_index
                        chain[outside_index][2][1] = min_index



    density0 = N - chain[:,2,0].tolist().count(-10) # Count the number of empty sites
    density1 = N - chain[:,2,1].tolist().count(-10) # Count the number of empty sites
    parity0 = chain[:, 0, 0].sum()
    parity1 = chain[:, 0 ,1].sum()
    # density = [density0, density1]
    return chain, density0 , density1, parity0, parity1

def evolution_no_annihilate_2L(chain, t, p):
    """
    Simulates the evolution of the Majorana chain over time t.
    The function modifies the chain in place.
    """
    density = np.zeros((t, 2), dtype=int)  # Initialize density array for two Majorana operators
    parity = np.zeros((t, 2), dtype =int)
    times = [k for k in range(t)]
    for k in range(t):

        chain, density0, density1, parity0, parity1= hopping_no_annihilate_2L(chain, p)
        density[k][0] = density0
        density[k][1] = density1
        parity[k][0] = parity0
        parity[k][1] = parity1
    return chain, times, density, parity

def calculate_pair_distances_braiding(chain):
    """
    Computes the arc length of Majorana pairs

    Returns:
      pair_data: a list of tuples (pair_index, distance)
      distance_counts: a dictionary mapping distance -> count of pairs with that distance
    """
    H = chain.shape[2]
    L = chain.shape[0]

    pair_positions = [defaultdict(list) for _ in range(H)]
    # distance_counts = {}
    distance_counts = defaultdict(int)
    # Loop over each row in the chain and record its index based on its pairing index (third column)
    for h in range(H):
        for i, row in enumerate(chain[:,:,h]):
            pair_idx = row[2]
            if pair_idx == -10:
                pass
            if pair_idx in pair_positions[h]:
                pair_positions[h][pair_idx].append(i)
            else:
                pair_positions[h][pair_idx] = [i]

        # pair_data = []


        for pair_idx, positions in pair_positions[h].items():
            if len(positions) == 2:
                d = min ( abs(positions[1] - positions[0]) , abs(positions[1]-positions[0]+ L),   abs(positions[1]-positions[0]- L) )
                # pair_data.append((pair_idx, d))
                distance_counts[d] = distance_counts.get(d, 0) +1
            else:
                # If there are not exactly two entries for a pair index, issue a warning.
                print(f"Warning: Pair index {pair_idx} appears {len(positions)} times (expected 2).")

    return distance_counts
    
def chain_builder_NL(N, rho, n):
    """
    having chain be a numpy array that the first column is the ferminonic number of the pair and the second column 
    is the Majorana operator which can be assigned randomly at contact to simulate quantum randomness.
    third column is the pair index
    """
    if N % 2 != 0:
        raise ValueError("N must be even for a Majorana chain.")
    chain = np.full((N * rho, 3, n), -10)
    for j in range (n):
        for i in range(N):
            chain[i * rho][0][j] = 1 # initially all pairs have fermion number 1
            chain[i * rho][1][j] = 47
            chain[i * rho][2][j] = int(i/2) # index for pair (i,j), where i<j, using i to denote the Majorana j is paired with.
    return chain

def hopping_annihilate_NL(chain, p):
    N = chain.shape[0]  # Number of sites in the chain
    H = chain.shape[2]
    Nt = []
    for h in range(H):
        Nt.append( chain[:,1,h].tolist().count(47))  # get the number of Majorana operators
    Nt_total = sum(Nt)
    # print(Nt)
    # prob_bin = []
    # for h in range(H): 
    #     prob_h = p ** (h) * (1- p) ** (H-h-1)

    for _ in range(Nt_total):
        
        particle_prob = random.random()
        cumulative = 0.0
        chosen_layer = -10
        for m in range(H):
            cumulative += Nt[m]/ Nt_total
            if particle_prob <= cumulative :
                chosen_layer= m
                # print(chosen_layer)
                break
        density = []
        for h in range(H):
            density_h = N - chain[:,2,h].tolist().count(-10) # Count the number of empty sites
            density.append(density_h)
            if any(d == 0 for d in density): #if all sites are empty
                return chain, density       
         
        index_m = np.where(chain[:,0, chosen_layer]!=-10)[0]
        i = random.choice(index_m)
        hop_direction0 = random.choice([-1, 1])
        
        
        if chain[i][0][chosen_layer] == -10: #if the site is empty
            pass

        elif chain[(i+hop_direction0) % N][0][chosen_layer] == -10: #if the site to arrive is empty
            
            r = random.random()
            for h in range(H):
                if h < chosen_layer:
                    if r < p:
                        if hop_direction0 >0 and chain[(i+hop_direction0) % N ][2][h] != -10: #if the site is empty
                                        # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[(i+hop_direction0) % N][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                            
                            if partner_j == -10:
                                raise ValueError("Partner Majorana is -10, check the chain structure.")
                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[(i+hop_direction0) % N][0][h] = (chain[(i+hop_direction0) % N][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                        elif hop_direction0 <0 and chain[i][2][h] != -10:
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[i][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != i][0]

                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[i][0][h] = (chain[i][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                if h == chosen_layer:
                    pass
                if h > chosen_layer:
                    if r < p:
                        if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][h] != -10: #if the site is empty
                                        # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[(i+hop_direction0) % N][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                            
                            if partner_j == -10:
                                raise ValueError("Partner Majorana is -10, check the chain structure.")
                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[(i+hop_direction0) % N][0][h] = (chain[(i+hop_direction0) % N][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                        elif hop_direction0 >0 and chain[i][2][h] != -10:
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[i][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != i][0]

                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[i][0][h] = (chain[i][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                    

            chain[(i+hop_direction0) % N][0][chosen_layer] = chain[i][0][chosen_layer]
            chain[i][0][chosen_layer] = -10
            chain[(i+hop_direction0) % N][1][chosen_layer] = chain[i][1][chosen_layer]
            chain[i][1][chosen_layer] = -10
            chain[(i+hop_direction0) % N][2][chosen_layer] = chain[i][2][chosen_layer]
            chain[i][2][chosen_layer] = -10

        elif chain[(i+hop_direction0) % N][2][chosen_layer] == chain[i][2][chosen_layer]:  #if paired under periodic boundary condition

            r = random.random()
            for h in range(H):
                if h < chosen_layer:
                    if r < p:
                        if hop_direction0 >0 and chain[(i+hop_direction0) % N ][2][h] != -10: #if the site is empty
                                        # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[(i+hop_direction0) % N][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                            
                            if partner_j == -10:
                                raise ValueError("Partner Majorana is -10, check the chain structure.")
                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[(i+hop_direction0) % N][0][h] = (chain[(i+hop_direction0) % N][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                        elif hop_direction0 <0 and chain[i][2][h] != -10:
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[i][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != i][0]

                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[i][0][h] = (chain[i][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                if h == chosen_layer:
                    pass
                if h > chosen_layer:
                    if r < p:
                        if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][h] != -10: #if the site is empty
                                        # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[(i+hop_direction0) % N][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                            
                            if partner_j == -10:
                                raise ValueError("Partner Majorana is -10, check the chain structure.")
                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[(i+hop_direction0) % N][0][h] = (chain[(i+hop_direction0) % N][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                        elif hop_direction0 >0 and chain[i][2][h] != -10:
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[i][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != i][0]

                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[i][0][h] = (chain[i][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                    

            if chain[i][0][chosen_layer] == 0 and chain[(i+hop_direction0) % N][0][chosen_layer] == 0:
                
                chain[i][0][chosen_layer] = -10
                chain[(i+hop_direction0) % N][0][chosen_layer] = -10
                chain[i][1][chosen_layer] = -10
                chain[(i+hop_direction0) % N][1][chosen_layer] = -10
                chain[i][2][chosen_layer] = -10
                chain[(i+hop_direction0) % N][2][chosen_layer] = -10
            else:
                r = random.random()
                for h in range(H):
                    if h < chosen_layer:
                        if r < p:
                            if hop_direction0 >0 and chain[(i+hop_direction0) % N ][2][h] != -10: #if the site is empty
                                            # get the index of paired Majorana
                                pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                                pair_j = np.where(chain[:, 2, h] == chain[(i+hop_direction0) % N][2][h])[0]
                                partner_i = [x for x in pair_i if x != i][0]
                                partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                                
                                if partner_j == -10:
                                    raise ValueError("Partner Majorana is -10, check the chain structure.")
                                # flip the fermionic parity
                                chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                                chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                                chain[(i+hop_direction0) % N][0][h] = (chain[(i+hop_direction0) % N][0][h]+ 1) % 2
                                chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                            elif hop_direction0 <0 and chain[i][2][h] != -10:
                                # get the index of paired Majorana
                                pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                                pair_j = np.where(chain[:, 2, h] == chain[i][2][h])[0]
                                partner_i = [x for x in pair_i if x != i][0]
                                partner_j = [x for x in pair_j if x != i][0]

                                # flip the fermionic parity
                                chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                                chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                                chain[i][0][h] = (chain[i][0][h]+ 1) % 2
                                chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                    if h == chosen_layer:
                        pass
                    if h > chosen_layer:
                        if r < p:
                            if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][h] != -10: #if the site is empty
                                            # get the index of paired Majorana
                                pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                                pair_j = np.where(chain[:, 2, h] == chain[(i+hop_direction0) % N][2][h])[0]
                                partner_i = [x for x in pair_i if x != i][0]
                                partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                                
                                if partner_j == -10:
                                    raise ValueError("Partner Majorana is -10, check the chain structure.")
                                # flip the fermionic parity
                                chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                                chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                                chain[(i+hop_direction0) % N][0][h] = (chain[(i+hop_direction0) % N][0][h]+ 1) % 2
                                chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                            elif hop_direction0 >0 and chain[i][2][h] != -10:
                                # get the index of paired Majorana
                                pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                                pair_j = np.where(chain[:, 2, h] == chain[i][2][h])[0]
                                partner_i = [x for x in pair_i if x != i][0]
                                partner_j = [x for x in pair_j if x != i][0]

                                # flip the fermionic parity
                                chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                                chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                                chain[i][0][h] = (chain[i][0][h]+ 1) % 2
                                chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2
                        
                
        
        elif chain[i][2][chosen_layer] !=  chain[(i+hop_direction0) % N][2][chosen_layer] :
    
    
            # get the index of paired Majorana
            rows_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]

            if len(rows_i) == 0:
                raise ValueError("No paired Majorana found for site i.")
            
            rows = np.where(chain[:, 2, chosen_layer] == chain[(i+hop_direction0) % N][2][chosen_layer])[0]
            outside_index_i = [x for x in rows_i if x != i][0]
            outside_index = [x for x in rows if x != (i+hop_direction0)%N][0]

            # braiding
            r = random.random()
            for h in range(H):
                if h < chosen_layer:
                    if r < p:
                        if hop_direction0 >0 and chain[(i+hop_direction0) % N ][2][h] != -10: #if the site is empty
                                        # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[(i+hop_direction0) % N][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                            
                            if partner_j == -10:
                                raise ValueError("Partner Majorana is -10, check the chain structure.")
                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[(i+hop_direction0) % N][0][h] = (chain[(i+hop_direction0) % N][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                        elif hop_direction0 <0 and chain[i][2][h] != -10:
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[i][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != i][0]

                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[i][0][h] = (chain[i][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                if h == chosen_layer:
                    pass
                if h > chosen_layer:
                    if r < p:
                        if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][h] != -10: #if the site is empty
                                        # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[(i+hop_direction0) % N][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                            
                            if partner_j == -10:
                                raise ValueError("Partner Majorana is -10, check the chain structure.")
                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[(i+hop_direction0) % N][0][h] = (chain[(i+hop_direction0) % N][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                        elif hop_direction0 >0 and chain[i][2][h] != -10:
                            # get the index of paired Majorana
                            pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                            pair_j = np.where(chain[:, 2, h] == chain[i][2][h])[0]
                            partner_i = [x for x in pair_i if x != i][0]
                            partner_j = [x for x in pair_j if x != i][0]

                            # flip the fermionic parity
                            chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                            chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                            chain[i][0][h] = (chain[i][0][h]+ 1) % 2
                            chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                    

            pair_index = chain[(i+hop_direction0) % N][2][chosen_layer]
            pair_index_i = chain[i][2][chosen_layer]
            min_index = min(pair_index, pair_index_i)
            max_index = max(pair_index, pair_index_i)

            # repairing at collision after 1st braid and calculate the sum of parity
            n_sum = chain[i][0][chosen_layer] + chain[(i+hop_direction0) % N][0][chosen_layer]
            ferminon_number = random.choice([0, 1])
            chain[i][0][chosen_layer] = ferminon_number
            chain[(i+hop_direction0) % N][0][chosen_layer] = ferminon_number

            if chain[i][0][chosen_layer] == 0 and chain[(i+hop_direction0) % N][0][chosen_layer] == 0:
                # annihilate the collided pair
                chain[i][0][chosen_layer] = -10
                chain[i][1][chosen_layer] = -10
                chain[i][2][chosen_layer] = -10
                chain[(i+hop_direction0) % N][0][chosen_layer] = -10
                chain[(i+hop_direction0) % N][1][chosen_layer] = -10
                chain[(i+hop_direction0) % N][2][chosen_layer] = -10

                # update the partner pair
                chain[outside_index_i][0][chosen_layer] = n_sum % 2
                chain[outside_index_i][1][chosen_layer] = 47
                chain[outside_index_i][2][chosen_layer] = min_index
                chain[outside_index][0][chosen_layer] = n_sum % 2
                chain[outside_index][1][chosen_layer] = 47
                chain[outside_index][2][chosen_layer] = min_index

        

            elif chain[i][0][chosen_layer] == 1 and chain[(i+hop_direction0) % N][0][chosen_layer] == 1:
                r = random.random()
                for h in range(H):
                    if h < chosen_layer:
                        if r < p:
                            if hop_direction0 >0 and chain[(i+hop_direction0) % N ][2][h] != -10: #if the site is empty
                                            # get the index of paired Majorana
                                pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                                pair_j = np.where(chain[:, 2, h] == chain[(i+hop_direction0) % N][2][h])[0]
                                partner_i = [x for x in pair_i if x != i][0]
                                partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                                
                                if partner_j == -10:
                                    raise ValueError("Partner Majorana is -10, check the chain structure.")
                                # flip the fermionic parity
                                chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                                chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                                chain[(i+hop_direction0) % N][0][h] = (chain[(i+hop_direction0) % N][0][h]+ 1) % 2
                                chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                            elif hop_direction0 <0 and chain[i][2][h] != -10:
                                # get the index of paired Majorana
                                pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                                pair_j = np.where(chain[:, 2, h] == chain[i][2][h])[0]
                                partner_i = [x for x in pair_i if x != i][0]
                                partner_j = [x for x in pair_j if x != i][0]

                                # flip the fermionic parity
                                chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                                chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                                chain[i][0][h] = (chain[i][0][h]+ 1) % 2
                                chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                    if h == chosen_layer:
                        pass
                    if h > chosen_layer:
                        if r < p:
                            if hop_direction0 <0 and chain[(i+hop_direction0) % N ][2][h] != -10: #if the site is empty
                                            # get the index of paired Majorana
                                pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                                pair_j = np.where(chain[:, 2, h] == chain[(i+hop_direction0) % N][2][h])[0]
                                partner_i = [x for x in pair_i if x != i][0]
                                partner_j = [x for x in pair_j if x != (i+hop_direction0)%N][0]
                                
                                if partner_j == -10:
                                    raise ValueError("Partner Majorana is -10, check the chain structure.")
                                # flip the fermionic parity
                                chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                                chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                                chain[(i+hop_direction0) % N][0][h] = (chain[(i+hop_direction0) % N][0][h]+ 1) % 2
                                chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2

                            elif hop_direction0 >0 and chain[i][2][h] != -10:
                                # get the index of paired Majorana
                                pair_i = np.where(chain[:, 2, chosen_layer] == chain[i][2][chosen_layer])[0]
                                pair_j = np.where(chain[:, 2, h] == chain[i][2][h])[0]
                                partner_i = [x for x in pair_i if x != i][0]
                                partner_j = [x for x in pair_j if x != i][0]

                                # flip the fermionic parity
                                chain[i][0][chosen_layer] = (chain[i][0][chosen_layer] +1 ) % 2
                                chain[partner_i][0][chosen_layer] = (chain[partner_i][0][chosen_layer] + 1) % 2
                                chain[i][0][h] = (chain[i][0][h]+ 1) % 2
                                chain[partner_j][0][h] = (chain[partner_j][0][h]+ 1) % 2
                            

                # calculate the ferminon number of new arcs
                chain[outside_index_i][0][chosen_layer] = (n_sum - chain[i][0][chosen_layer]) % 2
                chain[outside_index][0][chosen_layer] = (n_sum -  chain[(i+hop_direction0) % N][0][chosen_layer]) % 2
                

                if  i == min(i, (i+hop_direction0) % N, outside_index_i, outside_index) or (i+hop_direction0) % N == min(i, (i+hop_direction0) % N, outside_index_i, outside_index):
                    chain[i][2][chosen_layer] = min_index
                    chain[(i+hop_direction0) % N][2][chosen_layer] = min_index
                    chain[outside_index_i][2][chosen_layer] = max_index
                    chain[outside_index][2][chosen_layer] = max_index
                else:
                    chain[i][2][chosen_layer] = max_index
                    chain[(i+hop_direction0) % N][2][chosen_layer] = max_index
                    chain[outside_index_i][2][chosen_layer] = min_index
                    chain[outside_index][2][chosen_layer] = min_index
    
    
    # density = [density0, density1]
    return chain, density

def evolution_annihilate_NL(chain_, t, p):
    """
    Simulates the evolution of the Majorana chain over time t.
    The function modifies the chain in place.
    """
    chain = chain_.copy()  # Create a copy of the chain to avoid modifying the original
    H = chain.shape[2]  # Number of layers
    densities = np.zeros((t, H), dtype=int)  # Initialize density array for two Majorana operators
    times = []
    for k in range(t):
        times.append(k)
        chain, density= hopping_annihilate_NL(chain, p)
        if any(d ==0 for d in density):
            break
        for h in range(H):
            densities[k][h] = density[h]
        
    return times, densities

def calculate_pair_distances(chain):
    """
    Computes the arc length of Majorana pairs

    Returns:
      pair_data: a list of tuples (pair_index, distance)
      distance_counts: a dictionary mapping distance -> count of pairs with that distance
    """
    pair_positions = {}
    # Loop over each row in the chain and record its index based on its pairing index (third column)
    for i, row in enumerate(chain):
        pair_idx = row[2]
        if pair_idx == -10:
            pass
        if pair_idx in pair_positions:
            pair_positions[pair_idx].append(i)
        else:
            pair_positions[pair_idx] = [i]

    pair_data = []
    distance_counts = defaultdict(int)
    L = len(chain)
    for pair_idx, positions in pair_positions.items():
        if len(positions) == 2:
            d = min ( abs(positions[1] - positions[0]) , abs(positions[1]-positions[0]+ L),   abs(positions[1]-positions[0]- L) ) # Calculate the distance considering periodic boundary conditions
            pair_data.append((pair_idx, d))
            distance_counts[d] = distance_counts.get(d, 0)  +1
        else:
            # If there are not exactly two entries for a pair index, issue a warning.
            print(f"Warning: Pair index {pair_idx} appears {len(positions)} times (expected 2).")

    return pair_data, distance_counts

def linear_fit(x, y):
    """
    Fit a line y = ax + b to the data using least squares.
    
    Parameters:
        x (array-like): Independent variable data.
        y (array-like): Dependent variable data.
        
    Returns:
        a (float): Slope of the line.
        b (float): Intercept of the line.
    """
    x = np.asarray(x)
    y = np.asarray(y)

# Use np.polyfit for degree 1 polynomial
    a, b = np.polyfit(x, y, 1)
    return a, b   

def chain_compact(chain):
    """
    Remove all empty sites 1D.

    Returns:
      compact_chain: numpy array
          A 2D array with the first column as Majorana indices and the second column as fermion numbers.
    """
    index = np.where(chain[:, 1] == 47)[0]  # Get indices of Majoranas that are present
    N = len(index)
    compact_chain = np.zeros((N, 3), dtype=int)
    for i in range(len(index)):
        compact_chain[i][0] = chain[index[i]][0]  # Majorana index
        compact_chain[i][2] = chain[index[i]][2]   # Fermion number
    return compact_chain

def chain_compact_H(chain):
    """
    Remove all empty sites N-layer.

    Returns:
      compact_chain: numpy array
          A 2D array with the first column as Majorana indices and the second column as fermion numbers.
    """
    H = chain.shape[2]
    N = chain.shape[0]
    Nt =[]
    for h in range(H):
      index = np.where(chain[:, 1, h] == 47)[0]  # Get indices of Majoranas that are present
      Nt.append(len(index))
    Nt_short = max(Nt)
    compact_chain = np.zeros((Nt_short, 3, H), dtype=int)
    for h in range(H):
      index = np.where(chain[:, 1, h] == 47)[0]  # Get indices of Majoranas that are present
      Nt = len(index)
      
      for i in range(len(index)):
          compact_chain[i][0][h] = chain[index[i]][0][h]  # Majorana index
          compact_chain[i][2][h] = chain[index[i]][2][h]   # Fermion number
    return compact_chain


def main_compare():
# def main():
    chain1 = chain_builder(1000, 4)
    # print(np.where(chain[:, 1] == 47)[0])
    # print(500 - chain[:,2].tolist().count(-10))
    times1, density1 = evolution_annihilate(chain1, t =100000)
    chain2 = chain_builder(1000, 4)
    times2, density2 = evolution_classical(chain2, t = 100000)
    chain3 = chain_builder_2L(1000, 4)
    times3, density3 = evolution_annihilate_2L(chain3, 100000, 0.5)
    chain4 = chain_builder_2L(1000, 4)
    times4, density4 = evolution_annihilate_2L(chain4, 100000, 1)
    tim1 = np.array(times1)
    density_ref1 = 2/(np.sqrt(4 * np.pi* tim1[1:])) 
    tim2 = np.array(times2)
    density_ref2 = 1/(np.sqrt(4 * np.pi* tim2[1:])) 
    tim3 = np.array(times3)
    density_ref3 =  3/(np.sqrt(4 * np.pi* tim3[1:])) 
    tim4 = np.array(times4)
    density_ref4 =  4/(np.sqrt(4 * np.pi* tim4[1:])) 
    plt.figure(figsize =[10,6], dpi = 100 )
    plt.plot(times1[1:], np.array(density1[1:])/4000, label='Quantum Evolution', color = 'red')
    plt.plot(times2[1:], np.array(density2[1:])/4000, label='Classical Evolution', color = 'blue')
    plt.plot(times3[1:], (density3[1:, 1]+ density3[1:, 0])/4000, label='2 Layer Evolution p=0.5', color = 'orange')
    plt.plot(times4[1:], (density4[1:, 1]+ density4[1:, 0])/4000, label='2 Layer Evolution p=1', color = 'brown')
    plt.plot(tim1[1:], density_ref1, '--', label= r'$\frac{1}{\sqrt{2 \pi D t}}$', color = 'green')
    plt.plot(tim2[1:], density_ref2, '--', label=r'$\frac{1}{\sqrt{8 \pi D t}}$', color = 'green')
    plt.plot(tim3[1:], density_ref3, '--', label=r'$\frac{3}{\sqrt{8 \pi D t}}$', color = 'green')
    plt.plot(tim4[1:], density_ref4, '--', label=r'$\frac{4}{\sqrt{8 \pi D t}}$', color = 'green')
    plt.xlabel('Time')
    plt.ylabel('Density')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Density Decay of Majorana Chain')
    plt.legend()
    plt.savefig('density_evolution.png')
    plt.show()

    a1, b1 = linear_fit(np.log(times1[100:80000]), np.log(density1[100:80000]))
    a2, b2 = linear_fit(np.log(times2[100:80000]), np.log(density2[100:80000]))
    a3, b3 = linear_fit(np.log(times3[100:80000]), np.log(density3[100:80000]))
    
    # slope1 = (np.log(density1[10]) -np.log(density1[10000]) )/ (np.log(times1[10]) - np.log(times1[10000]) )
    # slope2 = (np.log(density2[10]) -np.log(density2[10000]) )/ (np.log(times2[10]) - np.log(times2[10000]) )
    # slope3 = (np.log(density3[10]) -np.log(density3[2000]) )/ (np.log(times3[10]) - np.log(times3[2000]) )
    # print("slope quantum:", slope1)
    # print("slope classical:", slope2)
    # print("slope 2L:", slope3)
    print("slope quantum:", a1, b1)
    print("slope classical:", a2, b2)
    print("slope 2L:", a3, b3)
    

# def main():
def main2L_single():
    chain1 = chain_builder_2L(1000, 4)
    times1, density1 = evolution_annihilate_2L(chain1, 20000, 0)
    chain2 = chain_builder_2L(1000, 4)
    times2, density2 = evolution_annihilate_2L(chain2, 20000, 0.25)
    chain3 = chain_builder_2L(1000, 4)
    times3, density3 = evolution_annihilate_2L(chain3, 20000, 0.5)
    chain4 = chain_builder_2L(1000, 4)
    times4, density4 = evolution_annihilate_2L(chain4, 20000, 0.75)
    chain5 = chain_builder_2L(1000, 4)
    times5, density5 = evolution_annihilate_2L(chain5, 20000, 1)
    
    density_ref1 = 4/(math.sqrt(3 * math.pi) * np.sqrt(times1[1:])) 
    density_ref2 = 3/(math.sqrt(3 * math.pi) * np.sqrt(times1[1:])) 
    # plt.plot(times1, (density1[:, 1]+ density1[:, 0]), label='Total Density', color='blue')
    
    plt.plot(times1, (density1[:, 1]+ density1[:, 0])/4000, label='Total Density p=0', color='blue')
    plt.plot(times2, (density2[:, 1]+ density2[:, 0])/4000, label='Total Density p=0.25', color='purple')
    plt.plot(times3, (density3[:, 1]+ density3[:, 0])/4000, label='Total Density p=0.5', color='orange')
    plt.plot(times4, (density4[:, 1]+ density4[:, 0])/4000, label='Total Density p=0.75', color='green')
    plt.plot(times5, (density5[:, 1]+ density5[:, 0])/4000, label='Total Density p=1', color='red')
    plt.plot(times1[1:], density_ref1, label= r'$\frac{4}{\sqrt{8 \pi D t}}$', color='black', linestyle='--')
    plt.plot(times1[1:], density_ref2, label= r'$\frac{3}{\sqrt{8 \pi D t}}$', color='black', linestyle='--')
    plt.xlabel('Time')
    plt.ylabel('Density')
    plt.xscale('log')
    plt.yscale('log')
    plt.title(f'Evolution of Majorana Chain Density {len(chain1)} ')
    # plt.xlim(0,6)
    plt.legend()
    plt.grid()
    plt.savefig('density decay of multiple exchange rate.png')
    plt.show()
    
def main_test_cl():
    chain1 = chain_builder(4000, 4)
    chain2 = chain_builder(4000, 4)
    chain3 = chain_builder(4000, 4)
    chain4 = chain_builder(4000, 4)
    times1, density1 = evolution_classical(chain1, 40000, 0.4)
    times2, density2 = evolution_classical(chain2, 40000, 0.6)
    times3, density3 = evolution_classical(chain3, 40000, 0.8)
    times4, density4 = evolution_classical(chain4, 40000, 1)
    plt.figure(figsize =[10,6], dpi = 100 )
    plt.plot(times1[1:], np.array(density1[1:]), label='p=0.4', color = 'red')
    plt.plot(times2[1:], np.array(density2[1:]), label='p=0.6', color = 'blue')
    plt.plot(times3[1:], np.array(density3[1:]), label='p=0.8', color = 'orange')
    plt.plot(times4[1:], np.array(density4[1:]), label='p=1', color = 'purple')
    plt.xlabel('Time')
    plt.ylabel('Density')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Classical Density of Majorana Chain')
    plt.legend()
    plt.savefig('classical_density.png')
    plt.show()

def main_17a():
# def main():
    # === Setup ===
    num_seeds = 5
    N = 1000
    rho = 3
    times = [1, 4, 10, 20, 50, 100, 150]
    
    # === Containers for averaging ===
    all_counts = [defaultdict(list) for _ in times]  # one defaultdict per time
    
    
    # === Main loop over random seeds ===
    for seed in range(num_seeds):
        np.random.seed(seed)
        chain = chain_builder_2L(N, rho)
    
        for i, t in enumerate(times):
            evolved_chain = evolution(chain, t)
            compacted = chain_compact(evolved_chain)
            _, distance_counts = calculate_pair_distances(compacted)
    
            for d, count in distance_counts.items():
                all_counts[i][d].append(count)
    
    # === Process average results ===
    avg_counts = []
    
    for i, time_counts in enumerate(all_counts):
        avg = {d: np.sum(c_list)/num_seeds for d, c_list in time_counts.items()}
        avg_counts.append(avg)
    
    # === Plotting ===
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'yellow', 'pink']
    
    for i, (avg_count_dict, color) in enumerate(zip(avg_counts, colors)):
        distances = sorted(avg_count_dict.keys())
        counts = [avg_count_dict[d] / (N * rho) for d in distances]
        plt.plot(distances, counts, label=f't={times[i]}', color=color)

    plt.xlabel('Distance')
    plt.ylabel('Average Count / (N * ρ)')
    plt.title('Averaged Pair Distances over Seeds')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.grid()
    plt.show()
    plt.savefig('17(a).png')


def main_17b():
# def main():
    # Parameters
    num_seeds = 20
    N = 200
    rho = 4
    annihilation_times = [10, 100, 1000]
    colors = ['blue', 'red', 'green', 'orange']
    labels = ['Initial', 't=10', 't=100', 't=10000']
    
    # Accumulators for distance counts at each time point
    count_accumulators = [defaultdict(list) for _ in range(1 + len(annihilation_times))]
    
    # Run simulation across seeds
    for seed in range(num_seeds):
        np.random.seed(seed)
    
        chain = chain_builder(N, rho)
        chain1 = evolution(chain, 2000)
        chain1_compacted = chain_compact(chain1)
        # Initial distances before annihilation
        _, distance_counts = calculate_pair_distances(chain1_compacted)
        for d, c in distance_counts.items():
            count_accumulators[0][d].append(c)
    
        # Distances after annihilation at different times
        for i, t in enumerate(annihilation_times):
            evolved_chain, _, _ = evolution_annihilate(chain1, t)
            compacted = chain_compact(evolved_chain)
            _, distance_counts_ann = calculate_pair_distances(compacted)
            for d, c in distance_counts_ann.items():
                count_accumulators[i + 1][d].append(c)
    
    # Average counts over seeds
    averaged_counts = []
    for dist_dict in count_accumulators:
        avg = {d: np.sum(c_list)/num_seeds for d, c_list in dist_dict.items()}
        averaged_counts.append(avg)
    
    # === Print distances with average count > 1 ===
    print("\nDistances with average count > 1:")
    greater_than_one = []  # for storing
    for i, avg_dict in enumerate(averaged_counts):
        print(f"\n--- {labels[i]} ---")
        gt1 = {d: c for d, c in avg_dict.items() if c > 1.1}
        greater_than_one.append(gt1)
        for d, c in sorted(gt1.items()):
            print(f"Distance: {d}, Avg Count: {c:.2f}")
    
    # === Plotting ===
    for i, avg_dict in enumerate(averaged_counts):
        distances = sorted(avg_dict.keys())
        counts = [avg_dict[d] for d in distances]
        plt.plot(distances, counts, label=labels[i], color=colors[i])
    
    
    plt.xlabel('Distance')
    plt.ylabel('Average Count')
    plt.title(f'Averaged P(s) vs s over {num_seeds} seeds')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.show()
    plt.savefig('17(b).png')

    
# def main_test2L():
def main():

    num_trials = 40  # number of seeds/runs to average
    N = 400
    rho = 4
    layers = 4
    Tmax = 10000
    # p_values = [0, 0.25, 0.5, 0.75, 1]  # different probabilities for your parameter
    p_values = [0, 0.5, 1]

    # For storing average results
    all_results = {}

    for p in p_values:
        densities = []
        parities = []
    
        for _ in range(num_trials):
            chain = chain_builder_NL(N, rho, layers)
            times, density = evolution_annihilate_NL(chain, Tmax, p)
            total_density = density.sum(axis =1)
            # total_parity = parity[:, 0] + parity[:, 1]
            densities.append(total_density)
            # parities.append(total_parity)
    
        # Convert to numpy array and compute average and std
        densities = np.array(densities)/( N * rho) # shape: (num_trials, time_steps)
        mean_density = np.mean(densities, axis=0)
        std_density = np.std(densities, axis=0)
        mean_parity = np.mean(parities, axis=0)
        # Save for plotting
        all_results[p] = {
            'times': times,
            'mean': mean_density,
            'std': std_density,
            'parity' : mean_parity
        }

    # Reference density
    density_ref_4 =  (2* layers) / (math.sqrt(4 * math.pi) * np.sqrt(times[1:]))
    density_ref_3 =  (layers +1) / (math.sqrt(4 * math.pi) * np.sqrt(times[1:]))
    density_ref_2 =  (layers +2) / (math.sqrt(4 * math.pi) * np.sqrt(times[1:]))
    
    
    
    # Plotting
    # for p, color in zip(p_values, ['purple', 'blue', 'orange', 'red', 'green' ]):
    for p, color in zip(p_values, ['purple', 'blue', 'orange' ]):
        result = all_results[p]
        plt.plot(result['times'], result['mean'], label=f'Total Density p={p}', color=color)
        # plt.plot(result['times'], result['parity'], label=f'Total Density p={p}', color=color)
        # Optional: error band
        # plt.fill_between(result['times'], result['mean'] - result['std'], result['mean'] + result['std'],
        #                  color=color, alpha=0.2)
    
    plt.plot(times[1:], density_ref_4, label=r'Reference density $\frac{2*N}{\sqrt{4 \pi t}}$', color='red', linestyle='--')
    plt.plot(times[1:], density_ref_3, label=r'Reference density $\frac{N+1}{\sqrt{4 \pi t}}$', color='red', linestyle='--')
    plt.plot(times[1:], density_ref_2, label=r'Reference density $\frac{N+2}{\sqrt{4 \pi t}}$', color='red', linestyle='--')
    
    plt.xlabel('Time')
    # plt.ylabel('Density')
    plt.ylabel('parity')
    plt.xscale('log')
    plt.yscale('log')
    plt.title(f'Averaged Majorana Chain Density (L = {N} * {rho}, Trials = {num_trials}), with separated initial spacing of two layers')
    plt.legend()
    plt.grid()
    plt.savefig(f'{layers}-Layer averaged (1)', dpi = 100)
    plt.show()

if __name__ == "__main__":
    main()
