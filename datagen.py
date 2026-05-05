import csv
import random
from a_star import a_star_search
from environment import get_energy


#Ensures at least 40% of the beads are 'H'
def generate_solvable_sequence(length=15, h_ratio=0.4):
    while True:
        seq = "".join(random.choice("HP") for _ in range(length))
        if seq.count('H') / length >= h_ratio:
            return seq



def create_dataset(num_samples=200):
    with open('protein_data.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["sequence", "energy", "coords"])
        
        samples_done = 0
        while samples_done < num_samples:
            length = random.randint(12, 18)
            seq = generate_solvable_sequence(length)
            
            print(f"Sample {samples_done + 1}: Attempting {seq}...")
            best_state, energy = a_star_search(seq)
            
            if best_state:
                writer.writerow([seq, energy, best_state.coords])
                samples_done += 1
                print(f"Success! Energy: {energy}")
            else:
                print("Sequence too complex, skipping...")




if __name__ == "__main__":
    create_dataset(200)