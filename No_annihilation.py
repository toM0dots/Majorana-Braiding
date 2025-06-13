import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import random
import math

# Define the Pauli matrices
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
I = np.eye(2, dtype=complex)



def chain_builder(N, rho):
    """
    having chain be a numpy array that the first column is the ferminonic number of the pair and the second column 
    is the Majorana operator which can be assigned randomly at contact to simulate quantum randomness.
    third column is the pair index
    """
    if N % 2 != 0:
        raise ValueError("N must be even for a Majorana chain.")
    chain = np.empty((N * rho, 4), dtype=object)
    for i in range(N):
        chain[i * rho][0] = 1 # initially all pairs have fermion number 1
        chain[i * rho][1] = True
        chain[i * rho][2] = int(i/2) # index for pair (i,j), where i<j, using i to denote the Majorana j is paired with.
        chain[i * rho][3] = i
    return chain

def hopping(chain_):
    chain = chain_.copy()
    N = len(chain)
    index = np.where(chain[:,1] == True)[0] 
    # i =random.randint(0, N-1)
    for k in range(len(index)):
        i = random.choice(index)
        hope_direction = np.random.choice([-1, 1])
        if chain[i][0] == None: #if the site is empty
            pass

        elif chain[(i+hope_direction) % N][0] == None: #if the site to arrive is empty
            chain[(i+hope_direction) % N][0] = chain[i][0]
            chain[i][0] = None
            chain[(i+hope_direction) % N][1] = chain[i][1]
            chain[i][1] = None
            chain[(i+hope_direction) % N][2] = chain[i][2]
            chain[i][2] = None
            chain[(i+hope_direction) % N][3] = chain[i][3]
            chain[i][3] = None
            

        elif chain[(i+hope_direction) % N][2] == chain[i][2]:  #if paired under periodic boundary condition
            # if chain[i][0] == 0 and chain[(i+hope_direction) % N][0] == 0:
            #     chain[i][2] = None
            #     chain[(i+hope_direction) % N][2] = None
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
            """       
            else:
                
                # get the index of pairing Majorana
                rows, cols = np.where(chain == chain[i][2])
                rows_right, cols_left = np.where(chain == chain[(i+hope_direction) % N][2])
                outside_index = [x for x in rows if x != i][0]
                right_outside_index = [x for x in rows_right if x != (i+hope_direction)%N][0]

                # calculate the ferminon number of new arcs
                chain[i][1] = np.random.choice([0, 1])
                chain[(i+hope_direction) % N][1] = np.random.choice([0, 1])
                ferminon_number = (chain[i][1] + chain[(i+hope_direction) % N][1]) %2
                
                chain[i][0] = ferminon_number
                chain[(i+hope_direction) % N][0] = ferminon_number
                chain[outside_index][0] = (ferminon_number+1) % 2 # conservation of ferminon number
                chain[right_outside_index][0] = (ferminon_number+1) % 2

                # update the pairing index
                pair_index_right = chain[(i+hope_direction) % N][2]
                pair_index_left = chain[i][2]
                if  i == min(i, (i+hope_direction) % N, outside_index, right_outside_index) or (i+hope_direction) % N == min(i, (i+hope_direction) % N, outside_index, right_outside_index):
                    chain[i][2] = pair_index_left
                    chain[(i+hope_direction) % N][2] = pair_index_left
                    chain[outside_index][2] = pair_index_right
                    chain[right_outside_index][2] = pair_index_right
                else:
                    chain[i][2] = pair_index_right
                    chain[(i+hope_direction) % N][2] = pair_index_right
                    chain[outside_index][2] = pair_index_left
                    chain[right_outside_index][2] = pair_index_left

                """
    return chain
    

    
def worldline(chain):
    """
    Generates a worldline representation of the Majorana chain.
    The function returns a list of pairs (Majorana index, fermion number).
    """
    index = np.where(chain[:, 1] == True)[0]  # Get indices of Majoranas that are present
    wl = np.zeros((1, len(index)), dtype=object)  # Initialize an empty list for the worldline
    for i in index:
        wl[0][chain[i][3]] = i  # position of ith particle
    return wl
    
def evolution(chain, t, plot=False):
    """
    Simulates the evolution of the Majorana chain over time t.
    The function modifies the chain in place.
    """
    if plot: 
        N = len(np.where(chain[:, 1] == True)[0])  # Number of Majoranas
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
        if pair_idx is None:
            pass
        if pair_idx in pair_positions:
            pair_positions[pair_idx].append(i)
        else:
            pair_positions[pair_idx] = [i]
    
    pair_data = []    
    distance_counts = {} 
    L = len(chain)
    for pair_idx, positions in pair_positions.items():
        if len(positions) == 2:
            d = min ( abs(positions[1] - positions[0]) , abs(positions[1]-positions[0]+ L),   abs(positions[1]-positions[0]- L) ) # Calculate the distance considering periodic boundary conditions
            pair_data.append((pair_idx, d))
            distance_counts[d] = distance_counts.get(d, 0) + 1
        else:
            # If there are not exactly two entries for a pair index, issue a warning.
            print(f"Warning: Pair index {pair_idx} appears {len(positions)} times (expected 2).")
    
    return pair_data, distance_counts
     

def S_A(chain, R):
    """
    Calculates entanglement entropy S_A based on the distances of Majorana pairs.
    
    Returns:
      S_A: float
          The computed value of S_A.
    """
    S_A = 0.0
    x0 = int((len(chain)-R)/2)
    for i in range(x0, x0+R):
        ids = np.where(chain[:, 2] == chain[i][2])[0]
        other_majorana = [x for x in ids if x != i][0]
        if not(x0 < other_majorana <= x0 + R):  # Check if the paired Majorana is outside the region A
          # Check if the pair exists
            S_A += 1
        # else:
        #     S_A += 1/2
    return S_A

        
def main():
    chain = chain_builder(2000, 4)
    chain1 = evolution(chain, t =5000)
    chain2 = evolution(chain1, t = 8000)
    chain3 = evolution(chain2, t = 8000)
    chain4 = evolution(chain3, t = 30000)

    pair_data1, distance_counts1 = calculate_pair_distances(chain)
    pair_data2, distance_counts2 = calculate_pair_distances(chain1)
    pair_data3, distance_counts3 = calculate_pair_distances(chain2)
    pair_data4, distance_counts4 = calculate_pair_distances(chain3)
    pair_data5, distance_counts5 = calculate_pair_distances(chain4)

    distances1 = list(distance_counts1.keys())
    counts1 = [distance_counts1[d] for d in distances1]
    distances2 = list(distance_counts2.keys())
    counts2 = [distance_counts2[d] for d in distances2]
    distances3 = list(distance_counts3.keys())
    counts3 = [distance_counts3[d] for d in distances3]
    distances4 = list(distance_counts4.keys())
    counts4 = [distance_counts4[d] for d in distances4]

    pairs1 = sorted(zip(distances1, counts1), key=lambda pair: pair[0])
    distances1_sorted, counts1_sorted = map(list, zip(*pairs1))
    pairs2 = sorted(zip(distances2, counts2), key=lambda pair: pair[0])
    distances2_sorted, counts2_sorted = map(list, zip(*pairs2))
    pairs3 = sorted(zip(distances3, counts3), key=lambda pair: pair[0])
    distances3_sorted, counts3_sorted = map(list, zip(*pairs3))
    pairs4 = sorted(zip(distances4, counts4), key=lambda pair: pair[0])
    distances4_sorted, counts4_sorted = map(list, zip(*pairs4))
    plt.scatter(distances1_sorted, counts1_sorted, color= 'red', s=10)  # adjust size if needed
    plt.scatter(distances2_sorted, counts2_sorted, color= 'blue', s=10)  # adjust size if needed
    plt.scatter(distances3_sorted, counts3_sorted, color= 'green', s=10)  # adjust size if needed
    plt.scatter(distances4_sorted, counts4_sorted, color= 'purple', s=10)  # adjust size if needed
    # plt.scatter(distances5, counts5, color= 'orange', s=10)  # adjust size if needed
    plt.xlabel(r'$\rho \ell$')
    plt.ylabel(r'P($\ell$)')
    # plt.xlim(0, 50)
    # plt.ylim(0, 75)
    plt.yscale('log')
    plt.xscale('log')
    plt.title("Correlation: Distance vs Number of Majorana Pairs")
    plt.grid(True)
    plt.savefig("No_annihilation.png", dpi=300)
    plt.show()
