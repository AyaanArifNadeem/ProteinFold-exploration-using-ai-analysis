import heapq
import time
from environment import ProteinState, get_energy, PROTEIN_SEQ, plot_protein

#Change string in environment.py to change the string the algo runs on


def calculate_heuristic(state):
    current_length = len(state.coords)
    remaining_sequence = state.sequence[current_length:]
    
    h_count = remaining_sequence.count('H')
    return -1 * h_count



def a_star_search(initial_sequence):
    frontier = []
    start_state = ProteinState(initial_sequence)
    
    # g(n) is current energy, h(n) is heuristic, f(n) = g(n) + h(n)
    start_g = 0 
    start_h = calculate_heuristic(start_state)
    start_f = start_g + start_h
    
    tie_breaker = 0
    heapq.heappush(frontier, (start_f, tie_breaker, start_state))
    explored = set()
    states_explored = 0
    print(f"Starting A* Search for sequence: {initial_sequence}...")
    start_time = time.time()
    
    while frontier:
        current_f, _, current_state = heapq.heappop(frontier)
        states_explored += 1
        
        if current_state.is_goal():
            end_time = time.time()
            print(f"\nSearch took {end_time - start_time:.4f} seconds")
            print(f"Total states explored: {states_explored}")
            return current_state, get_energy(current_state.coords, initial_sequence)
            
        state_shape = tuple(current_state.coords)
        if state_shape in explored:
            continue 
        explored.add(state_shape)
        
        successors = current_state.generate_successors()
        for successor in successors:
            g_n = get_energy(successor.coords, initial_sequence)
            h_n = calculate_heuristic(successor)
            f_n = g_n + h_n
            tie_breaker += 1
            heapq.heappush(frontier, (f_n, tie_breaker, successor))
            
    return None, 0




# Code block to run A*
if __name__ == '__main__':
    best_fold, best_energy = a_star_search(PROTEIN_SEQ)

    if best_fold:
        print(f"Absolute Minimum Energy Found: {best_energy}")
        plot_protein(best_fold.coords, PROTEIN_SEQ)
    else:
        print("No valid folds possible.")