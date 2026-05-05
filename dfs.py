from environment import ProteinState, get_energy, PROTEIN_SEQ, plot_protein



# Main DFS ALGO
def depth_first_search(initial_sequence):
    start_state = ProteinState(initial_sequence)
    frontier = [start_state]
    best_state = None
    lowest_energy = float('inf')
    states_explored = 0
    
    print(f"Starting DFS for sequence: {initial_sequence}...")
    
    while frontier:
        current_state = frontier.pop()
        states_explored += 1
        
        if current_state.is_goal():
            current_energy = get_energy(current_state.coords, initial_sequence)
            
            if current_energy < lowest_energy:
                lowest_energy = current_energy
                best_state = current_state
                print(f"New best fold found! Energy: {lowest_energy}")
                
            continue
            
        successors = current_state.generate_successors()
        
        for successor in successors:
            frontier.append(successor)
            
    print(f"\nSearch Complete! Total states explored: {states_explored}")
    return best_state, lowest_energy



# Code block to test DFS
# TO change protein sequence go back to environment file and change sequence

best_fold, best_energy = depth_first_search(PROTEIN_SEQ)

if best_fold:
    print(f"Absolute Minimum Energy: {best_energy}")
    plot_protein(best_fold.coords, PROTEIN_SEQ)
else:
    print("No valid folds possible (the sequence might be trapped).")