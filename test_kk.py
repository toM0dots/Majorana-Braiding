import numpy as np
import random
import numpy.linalg as la
import matplotlib.pyplot as plt
import math
from collections import Counter
from typing import List

def chain_builder_classical(N, rho, r1):
    """
    having chain be a numpy array that the first column is the ferminonic number of the pair and the second column 
    is the Majorana operator which can be assigned randomly at contact to simulate quantum randomness.
    third column is the pair index
    """
    if N % 2 != 0:
        raise ValueError("N must be even for a Majorana chain.")
    
    chain = np.full((N* rho, 3), -10) 
       
    for i in range(N):
            
        rho_middle = int(rho * r1)
        # print(rho_middle)
            
        chain[(i * rho + rho_middle) % (N * rho)][0] = 1 # initially all pairs have fermion number 1
        chain[(i * rho + rho_middle) % (N * rho)][1] = 47
        chain[(i * rho + rho_middle) % (N * rho)][2] = int(i/2)
        chain[((i+1) * rho) % (N * rho)][0] = 1 # initially all pairs have fermion number 1
        chain[((i+1) * rho) % (N * rho)][1] = 47
        chain[((i+1) * rho) % (N * rho)][2] = int(i/2)
            
        rho_middle = int(rho * r1)
        # print(rho_middle)
            
        chain[(i * rho + rho_middle) % (N * rho)][0] = 1 # initially all pairs have fermion number 1
        chain[(i * rho + rho_middle) % (N * rho)][1] = 47
        chain[(i * rho + rho_middle) % (N * rho)][2] = int(i/2)
        chain[((i+1) * rho) % (N * rho)][0] = 1 # initially all pairs have fermion number 1
        chain[((i+1) * rho) % (N * rho)][1] = 47
        chain[((i+1) * rho) % (N * rho)][2] = int(i/2)


    return chain

def hopping_classical(chain):
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
            chain[i][0] = -10
            chain[(i+hope_direction) % N][0] = -10
            chain[i][1] = -10
            chain[(i+hope_direction) % N][1] = -10
            chain[i][2] = -10
            chain[(i+hope_direction) % N][2] = -10

                
    density = N - chain[:,2].tolist().count(-10)  # Count the number of empty sites
    return chain, density
    
    

def evolution_classical(chain, t, fraction = False):
    """
    Simulates the evolution of the Majorana chain over time t.
    The function modifies the chain in place.
    """
    density = []
    times = [i for i in range(t)]
    N = chain.shape[0] 
    fractions = np.zeros((t, 2))  # To store even and odd fractions if needed
    for _ in range(t):
        
        chain, x = hopping_classical(chain)
        density.append(x)
        if fraction:
            positions = np.where(chain[:,1]==47)[0]
            distances = [j - i for i, j in zip(positions, positions[1:])]
    
            wrap_dist = (positions[0] + N) - positions[-1]
            distances.append(wrap_dist)
            total_gap = sum(distances)
            even_gap_sum = sum(distances[0::2])  # gaps at indices 0,2,4…
            odd_gap_sum  = sum(distances[1::2])  # gaps at indices 1,3,5…
            even_fraction = even_gap_sum / total_gap if total_gap != 0 else 0
            odd_fraction = odd_gap_sum / total_gap if total_gap != 0 else 0
            fractions[_][0] = even_fraction
            fractions[_][1] = odd_fraction

    return chain, times, density, fractions if fraction else None





def chain_builder(N, rho, parity_ratio =1):
    """
    having chain be a numpy array that the first column is the ferminonic number of the pair and the second column 
    is the Majorana operator which can be assigned randomly at contact to simulate quantum randomness.
    third column is the pair index
    """
    if N % 2 != 0:
        raise ValueError("N must be even for a Majorana chain.")
    chain = np.full((N* rho, 3), -10)
    for i in range(N):
        if i % 2 == 0: # if i is even, then the pair is a fermion pair
            if random.random() < parity_ratio:
                parity = 0
            else:
                parity = 1
        chain[i * rho][0] = parity # initially all pairs have fermion number 1
        chain[i * rho][1] = 47
        chain[i * rho][2] = int(i/2) # index for pair (i,j), where i<j, using i to denote the Majorana j is paired with.
    return chain

def chain_builder_0(N, rho, parity):
    """
    having chain be a numpy array that the first column is the ferminonic number of the pair and the second column 
    is the Majorana operator which can be assigned randomly at contact to simulate quantum randomness.
    third column is the pair index
    """
    if N % 2 != 0:
        raise ValueError("N must be even for a Majorana chain.")
    chain = np.full((N* rho, 3), -10)
    for i in range(N):
        chain[i * rho][0] = parity # initially all pairs have fermion number 1
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
        if  chain[:, 2].tolist().count(-10) == N: #if all sites are empty

            return chain, 0
        
        index = np.where(chain[:,1]==47)[0]
        i = random.choice(index)

        hope_direction = random.choice([-1, 1])
        
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
    return chain, times, density

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
    distance_counts = {} 
    
    for pair_idx, positions in pair_positions.items():
        if len(positions) == 2:
            d = abs(positions[1] - positions[0])
            pair_data.append((pair_idx, d))
            distance_counts[d] = distance_counts.get(d, 0) + 1
        else:
            # If there are not exactly two entries for a pair index, issue a warning.
            print(f"Warning: Pair index {pair_idx} appears {len(positions)} times (expected 2).")
    
    return pair_data, distance_counts
     

# def main_quantum():
def main():
    N      = 1200
    rho    = 4
    t_max  = 16000
    seeds  = range(15)
    params = [0, 0.25, 0.5, 0.75, 0.85, 0.95, 1.0]
    # params = [1.0]

    # Prepare a dict to collect density arrays for each param
    densities = {p: [] for p in params}

    for seed in seeds:
        # make everything reproducible
        random.seed(seed)
        np.random.seed(seed)

        # build chains for this seed
        # chains = {
        #     0.0 : chain_builder(N, rho, 0.0),
        #     0.25: chain_builder(N, rho, 0.25),
        #     0.5 : chain_builder(N, rho, 0.5),
        #     0.75: chain_builder(N, rho, 0.75),
        #     0.85: chain_builder(N, rho, 0.85),
        #     0.95 : chain_builder(N, rho, 0.95),
        #     1.0 : chain_builder(N, rho, 1.0),
        # }

        chains = {
            p : chain_builder(N, rho, p) for p in params
            
        }
        for chain in chains.values():
            pass
            # chain[0*rho][0] = 1
            # chain[1*rho][0] = 1
            # chain[2*rho][0] = 1
            # chain[3*rho][0] = 1
            # chain[4*rho][0] = 1
            # chain[5*rho][0] = 1
            # chain[6*rho][0] = 1
            # chain[7*rho][0] = 1
            # chain[8*rho][0] = 1
        # run evolution and collect densities
        for p, chain in chains.items():
            _, times, density = evolution_annihilate(chain, t_max)
            densities[p].append(density)

    # Now compute the average density vs time for each parameter
    avg_densities = {}
    for p in params:
        arr = np.stack(densities[p], axis=0)   # shape = (100, len(times))
        avg = arr.mean(axis=0)                 # shape = (len(times),)
        avg_densities[p] = avg

    # Example: plot the 5 averaged curves
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12,9))
    for p in params:
        plt.plot(times, avg_densities[p]/(N*rho), label=f"parity-0 ratio ={p}")
    plt.plot(times[1:], 2/(np.sqrt(4 * np.pi * np.array(times[1:]))), linestyle='--', label = r'$\frac{2}{\sqrt{4 \pi t}}$')
    plt.plot(times[1:], 1.5/(np.sqrt(4 * np.pi * np.array(times[1:]))), linestyle='--', label = r'$\frac{1.5}{\sqrt{4 \pi t}}$')
    plt.xlabel("Time")
    plt.ylabel("⟨density⟩")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()

    plt.figure(figsize=(12,9))
    for p in params:
        plt.plot(times, avg_densities[p]*(np.sqrt(4 * np.pi * np.array(times)))/(N*rho), label=f"parity-0 ratio={p}")
    plt.plot(times, np.full(t_max, 2 ), linestyle='--', label = '2')
    plt.plot(times, np.full(t_max, 1.5),linestyle='--', label = '1.5')
    # plt.axvline(x = 0.125*(N*rho)**2, color='r', linestyle='--', label = 'D L^2')
    plt.xlabel("time")
    plt.ylabel("⟨density⟩")
    plt.xscale('log')
    # plt.yscale('log')
    plt.legend()
    plt.show()

# def main():
def main_classical():
    N      = 300
    rho    = 4
    t_max  = 3000
    trails = 30
    seeds  = range(trails)
    # params = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0]
    params = [0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5]
    # Prepare a dict to collect density arrays for each param
    densities = {p: [] for p in params}
    densities = []
    fractions = []

    for seed in seeds:
        # make everything reproducible
        random.seed(seed)
        np.random.seed(seed)

        # build chains for this seed
        chain = np.vstack((chain_builder_classical(N, rho, 0.25), chain_builder_classical(N, rho, 0.5)))
        _, times, density, frac = evolution_classical(chain, t_max, fraction=True)
        densities.append(density)
        fractions.append(frac)
        # chains = {
        #     0.0 : chain_builder_classical(N, rho, 0.0),
        #     0.1 : chain_builder_classical(N, rho, 0.1),
        #     0.2: chain_builder_classical(N, rho, 0.2),
        #     0.25: chain_builder_classical(N, rho, 0.25),
        #     0.3 : chain_builder_classical(N, rho, 0.3),
        #     0.4 : chain_builder_classical(N, rho, 0.4),
        #     0.5 : chain_builder_classical(N, rho, 0.5),
        #     # 0.6: chain_builder_classical(N, rho, 0.6),
        #     # 0.8: chain_builder_classical(N, rho, 0.8),
        #     # 0.9 : chain_builder_classical(N, rho, 0.9),
        #     # 1.0 : chain_builder_classical(N, rho, 1.0),
        #     # 0.8: chain_builder_classical(N, 2*rho, 0.2),
        #     # 0.9 : chain_builder_classical(N, 2*rho, 0.4),
        #     # 1.0: chain_builder_classical(N, 2*rho, 0.6),
        # }

        # run evolution and collect densities
        # for p, chain in chains.items():
        #     _, times, density = evolution_classical(chain, t_max)
        #     densities[p].append(density)

    # Now compute the average density vs time for each parameter
    avg_densities = {}
    for p in params:
        arr = np.stack(densities[p], axis=0)   # shape = (100, len(times))
        avg = arr.mean(axis=0)                 # shape = (len(times),)
        avg_densities[p] = avg

    # Example: plot the 5 averaged curves
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12,9))
    for p in params:
        plt.plot(times, avg_densities[p]/(N*rho), label=f"Dual magnetization {p}s, {1-p}s")
        plt.plot(times[1:], 4*(1-p)*p/(np.sqrt(4 * np.pi * np.array(times[1:]))), linestyle='--', label = fr'\frac{{4 {p: .2f} (1-{p: .2f})}}{{\sqrt{{4 \pi t}} }}')
    # plt.plot(times[1:], 0.75/(np.sqrt(4 * np.pi * np.array(times[1:]))), linestyle='--', label = '0.75/sqrt(4πt)')
    plt.title("Classical density vs time", fontsize=15, fontfamily='Times New Roman')
    plt.xlabel("Time",fontsize=16, fontfamily='Times New Roman')
    plt.ylabel(r"$\rho(t)$",fontsize=16, fontfamily='Times New Roman')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig("density_vs_time_classical.png")
    plt.show()

    plt.figure(figsize=(12,9))
    for p in params:
        plt.plot(times, avg_densities[p]*(np.sqrt(4 * np.pi * np.array(times)))/(N*rho), label=f"Spin up/down ratio={p: .2f}")
        plt.plot(times, np.full(t_max, 4*(1-p)*p) , linestyle='--', label = fr'$4 {p: .2f} (1-{p : .2f}) $')
    # plt.plot(times, np.full(t_max, 1 ), linestyle='--', label = '1')
    # plt.plot(times, np.full(t_max, 0.75),linestyle='--', label = '0.75')
    plt.title("Classical density Conserved Quantity", fontsize=20, fontfamily='Times New Roman')
    plt.xlabel("Time", fontsize=16, fontfamily='Times New Roman')
    plt.ylabel(r"\rho(t) \sqrt{8 \pi D t}$", fontsize=16, fontfamily='Times New Roman')
    plt.xscale('log')
    # plt.yscale('log')
    plt.legend()
    plt.savefig(f"density_vs_time_classical_conserved_quantity_{N}*{rho}_trails={trails}.png")
    plt.show()

def main():
# def main_ising():
    N      = 400
    rho    = 4
    t_max  = 10000
    seeds  = range(100)
    # params = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0]
    params = [0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5]
    # Prepare a dict to collect density arrays for each param
    densities = {p: [] for p in params}
    densities = []
    fractions = []

    for seed in seeds:
        # make everything reproducible
        random.seed(seed)
        np.random.seed(seed)

        # build chains for this seed
        chain = np.vstack((chain_builder_classical(N, rho, 0.25), chain_builder_classical(N, rho, 0.75)))
        chain = chain[:-2, :]
        _, times, density, frac = evolution_classical(chain, t_max, fraction=True)
        densities.append(density)
        fractions.append(frac)
        # chains = {
        #     0.0 : chain_builder_classical(N, rho, 0.0),
        #     0.1 : chain_builder_classical(N, rho, 0.1),
        #     0.2: chain_builder_classical(N, rho, 0.2),
        #     0.25: chain_builder_classical(N, rho, 0.25),
        #     0.3 : chain_builder_classical(N, rho, 0.3),
        #     0.4 : chain_builder_classical(N, rho, 0.4),
        #     0.5 : chain_builder_classical(N, rho, 0.5),
        # }

        # run evolution and collect densities
        # for p, chain in chains.items():
        #     _, times, density = evolution_classical(chain, t_max)
        #     densities[p].append(density)

    # Now compute the average density vs time for each parameter
    avg_densities = np.mean(densities, axis=0)  # shape = (100, len(times))
    avg_frac = np.mean(fractions, axis=0)  # shape = (100, len(times), 2)

    # Example: plot the 5 averaged curves
    import matplotlib.pyplot as plt

    plt.figure(figsize=(12,9))
    plt.plot(times, avg_densities/(2*N*rho), label=f"alternating initial conditions")
    plt.plot(times[1:], 1/(np.sqrt(4 * np.pi * np.array(times[1:]))), linestyle='--', label = f'1/sqrt(4πt)')
    plt.plot(times[1:], 15/(16*np.sqrt(4 * np.pi * np.array(times[1:]))), linestyle='--', label = f'15/16sqrt(4πt)')
    # plt.plot(times[1:], 0.75/(np.sqrt(4 * np.pi * np.array(times[1:]))), linestyle='--', label = '0.75/sqrt(4πt)')
    plt.xlabel("time")
    plt.ylabel("⟨density⟩")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    # plt.savefig("density_vs_time_classical.png")
    plt.show()

    avg_frac = np.array(avg_frac)
    plt.figure(figsize=(12,9))
    plt.plot(times, 4*avg_frac[:,0]*avg_frac[:,1], label=f"fraction 4x(1-x)")
    plt.axhline(y=15/16, color='r', linestyle='--', label='1')
    # plt.plot(times, np.full(t_max, 1 ), linestyle='--', label = '1')
    # plt.plot(times, np.full(t_max, 0.75),linestyle='--', label = '0.75')
    plt.xlabel("time")
    plt.ylabel("constant")
    plt.xscale('log')
    # plt.yscale('log')
    plt.legend()
    # plt.savefig("density_vs_time_classical(1).png")
    plt.show()

if __name__ == "__main__":
    main()