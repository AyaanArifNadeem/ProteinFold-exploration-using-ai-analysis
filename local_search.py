import random
import math
import time
import numpy as np
from environment import get_energy, plot_protein, PROTEIN_SEQ



# basically flipping about corners to priduce a new structure with lower energy
def get_v_flip_neighbor(coords):
    """Flips existing corners."""
    if len(coords) < 3: return None
    indices = list(range(1, len(coords) - 1))
    random.shuffle(indices)
    
    for idx in indices:
        p_prev = np.array(coords[idx-1])
        p_curr = np.array(coords[idx])
        p_next = np.array(coords[idx+1])
        
        dist_sq = np.sum((p_next - p_prev)**2)
        if dist_sq == 2: 
            new_pos = tuple((p_prev + p_next - p_curr).tolist())
            if new_pos not in coords:
                new_coords = list(coords)
                new_coords[idx] = new_pos
                return new_coords 
    return None 


# shifting position of tail to produce new combinations
def get_end_move_neighbor(coords):
    new_coords = list(coords)
    end_idx = random.choice([0, -1])
    anchor_idx = 1 if end_idx == 0 else -2
    anchor_pos = coords[anchor_idx]
    
    possible_moves = [
        (anchor_pos[0]+1, anchor_pos[1]), (anchor_pos[0]-1, anchor_pos[1]),
        (anchor_pos[0], anchor_pos[1]+1), (anchor_pos[0], anchor_pos[1]-1)
    ]
    
    valid_moves = [m for m in possible_moves if m not in coords]
    
    if valid_moves:
        new_coords[end_idx] = random.choice(valid_moves)
        return new_coords
    return None



def simulated_annealing(sequence, initial_coords, iterations=20000):
    start_time = time.time()
    current_coords = initial_coords
    current_energy = get_energy(current_coords, sequence)
    best_coords = current_coords
    best_energy = current_energy
    
    temp = 100.0
    cooling_rate = (0.01 / 100.0) ** (1.0 / iterations)
    
    print(f"Starting Local Search... Initial Energy: {current_energy}")
    
    for i in range(iterations):
        if random.random() < 0.5:
            neighbor = get_v_flip_neighbor(current_coords)
        else:
            neighbor = get_end_move_neighbor(current_coords)
            
        if neighbor:
            neighbor_energy = get_energy(neighbor, sequence)
            delta_e = neighbor_energy - current_energy
            
            if delta_e <= 0 or random.random() < math.exp(-delta_e / temp):
                current_coords = neighbor
                current_energy = neighbor_energy
                
                if current_energy < best_energy:
                    best_energy = current_energy
                    best_coords = current_coords
                    
        temp *= cooling_rate

    end_time = time.time()
    print(f"\nLocal Search took {end_time - start_time:.4f} seconds")
    print(f"Initial Energy: {get_energy(initial_coords, sequence)} | Best Energy Found: {best_energy}")
    
    return best_coords, best_energy





#block to run the code for LOCAL SEARCH
if __name__ == "__main__":
    initial_fold = [(i, 0) for i in range(len(PROTEIN_SEQ))] 
    
    final_coords, final_energy = simulated_annealing(PROTEIN_SEQ, initial_fold)
    plot_protein(final_coords, PROTEIN_SEQ)