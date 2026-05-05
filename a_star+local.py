import time
from a_star import a_star_search
from local_search import simulated_annealing
from environment import plot_protein, get_energy

PROTEIN_SEQ = 'HHPHPHPHPHPHHPHPHPHPHPHPHHPHP'

print("Running A* Search")
a_star_start_time = time.time()
a_star_result, a_star_energy = a_star_search(PROTEIN_SEQ)
a_star_end_time = time.time()
a_star_time_taken = a_star_end_time - a_star_start_time

if a_star_result:
    print(f"A* found a fold with energy: {a_star_energy}")
    print(f"A* Time Taken: {a_star_time_taken:.4f} seconds")
    
    
    # Create an initial fold (mostly straight)
    SA_initial_coords = [(i, 0) for i in range(len(PROTEIN_SEQ))]
    
    # Bend the very last bead up by 90 degrees to create a corner
    # this allows SA to run slightly better
    last_x, _ = SA_initial_coords[-2]
    SA_initial_coords[-1] = (last_x, 1) 
    SA_energy = get_energy(SA_initial_coords, PROTEIN_SEQ)
        
    print("\nRunning simulated annealing on straight line")
    sa_start_time = time.time()
    refined_coords, refined_energy = simulated_annealing(
        PROTEIN_SEQ, 
        SA_initial_coords, 
        iterations=50000 
    )
    sa_end_time = time.time()
    sa_time_taken = sa_end_time - sa_start_time
    
    print(f"\nRefinement Complete!")
    print(f"Original Straight Line Energy: {SA_energy} | Refined SA Energy: {refined_energy}")
    print(f"Simulated Annealing Time Taken: {sa_time_taken:.4f} seconds")
    