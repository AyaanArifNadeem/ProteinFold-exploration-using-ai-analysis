import numpy as np
import matplotlib.pyplot as plt
import random

# 'H' = Hydrophobic (Blue), 'P' = Polar (Red)
PROTEIN_SEQ = "HHPPHHPHPHPHHPHPHPHHPHHPP" 
ACTIONS = ['U', 'D', 'L', 'R']
MOVE_MAP = {'U': (0, 1), 'D': (0, -1), 'L': (-1, 0), 'R': (1, 0)}


# Returns legal moves
def get_valid_actions(current_coords):
    last_x, last_y = current_coords[-1]
    valid_moves = []
    
    for move in ACTIONS:
        dx, dy = MOVE_MAP[move]
        next_step = (last_x + dx, last_y + dy)
        
        if next_step not in current_coords:
            valid_moves.append(move)
            
    return valid_moves


# Returns energy value of protein struc (the lower the better)
def get_energy(coords, sequence):
    energy = 0
    
    for i in range(len(coords)):
        for j in range(i + 2, len(coords)):
            if sequence[i] == 'H' and sequence[j] == 'H':
                dist = abs(coords[i][0] - coords[j][0]) + abs(coords[i][1] - coords[j][1])
                if dist == 1: 
                    energy -= 1
    return energy


#
class ProteinState:
    def __init__(self, sequence, coords=None):
        self.sequence = sequence
        # Start at origin if no coordinates are provided
        self.coords = coords if coords else [(0, 0)]

    def is_goal(self):
        return len(self.coords) == len(self.sequence)

    def generate_successors(self):
        successors = []
        valid_moves = get_valid_actions(self.coords)
        
        for move in valid_moves:
            dx, dy = MOVE_MAP[move]
            new_x, new_y = self.coords[-1][0] + dx, self.coords[-1][1] + dy
            new_coords = self.coords + [(new_x, new_y)]
            successors.append(ProteinState(self.sequence, new_coords))
            
        return successors


# Visualizer
def plot_protein(coords, sequence):
    coords = np.array(coords)
    plt.figure(figsize=(6,6))
    
    plt.plot(coords[:,0], coords[:,1], color='black', zorder=1)
    
    for i, (x, y) in enumerate(coords):
        color = 'blue' if sequence[i] == 'H' else 'red'
        plt.scatter(x, y, color=color, s=200, edgecolors='black', zorder=2)
        plt.text(x, y, f"{i}", fontsize=12, ha='center', va='center', color='white')

    plt.grid(True)
    plt.title(f"Protein Fold Visualization (Sequence: {sequence})")
    plt.show()


# # Uncomment this to test enviornment using random fold genereation
# # Warning recomment this to allow other files to work properly
# current_state = ProteinState(PROTEIN_SEQ)

# print("Starting Random Fold Simulation...")

# while not current_state.is_goal():
#     successors = current_state.generate_successors()

#     if not successors:
#         print("Agent got stuck! Dead end reached.")
#         break
        
#     current_state = random.choice(successors)

# if current_state.is_goal():
#     print(f"Goal Reached! Final Energy: {get_energy(current_state.coords, PROTEIN_SEQ)}")
#     plot_protein(current_state.coords, PROTEIN_SEQ)