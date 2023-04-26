#First extract the "Bonds" section from data file and remove header (save as csv file).
#Only keep two columns for the two atom ids for each bond. Then run atom_parser_v2.py
#Next, delete all columns above and below the segment you need in the output.csv.
#Next run this script. Make sure to have the lammpstrj file in this directory.

import csv
import math

# Set up the parameters
target_atom_type = 1 # Atom type of target diffusing atom(s)
num_atoms = # Number of atoms in your system.
trj_file = # Name of your LAMMPS trajectory file.

# Function to read the .lammpstrj file
def read_lammpstrj(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    rows = []
    for line in lines[9:]:
        row = line.split()
        if row and len(row) >= 6 and row[0].isdigit():  # Check if the row has at least 6 elements and the first item is a number (assumes this is an atom entry)
            atom_id = int(row[0])
            atom_type = int(row[1])
            coords = tuple(float(val) for val in row[3:6])
            rows.append((atom_id, atom_type) + coords)
    print(len(rows))
    return rows

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

# Read in the .lammpstrj file
rows = read_lammpstrj(trj_file)

# Read in the csv file
chains = {}
with open('output.csv', 'r') as f:
    reader = csv.reader(f)
    for i, row in enumerate(reader):
        chains[f'chain_{i+1}'] = [int(x) for x in row]

# Initialize the dictionary to hold the distances of each target species
target_species_distances = {}

# Initialize the dictionary to hold the previous positions of each target species
previous_positions = {}

row_num = 0
x = 1 # Analyze every xth timestep in your simulation.
y = 201 # Max number of timesteps to analyze
atom_num = int(num_atoms + 9)

for loop_num in range(0, y, x):
    print(f"Processing loop {loop_num}")  # Progress update
    print(f"Processing from row {row_num}")  # Progress update

    for row in rows[row_num:row_num+atom_num]:
        if len(row) < 3:  # Skip the row if it has less than 3 elements
            continue
        atom_id, atom_type, x, y, z = row
        if atom_type == target_atom_type:
            species_id = atom_id
            coords = (x, y, z)

            # If species_id is not in target_species_distances, initialize it with distance 0
            if species_id not in target_species_distances:
                target_species_distances[species_id] = 0

            # If species_id is not in previous_positions, store its current position
            if species_id not in previous_positions:
                previous_positions[species_id] = coords
            else:
                # Calculate distance moved since last position
                distance = math.sqrt(sum([(a - b)**2 for a, b in zip(coords, previous_positions[species_id])]))
                if species_id == 113:
                    print(f"Species {species_id} - Old position: {previous_positions[species_id]} - New position: {coords} - Distance: {distance}")  # Print intermediate results
                target_species_distances[species_id] += distance
                previous_positions[species_id] = coords

    # Increment the row number to skip to the next loop
    row_num += atom_num

# Sort the distances in descending order and create a list of tuples
sorted_distances = sorted(target_species_distances.items(), key=lambda x: x[1], reverse=True)

# Write the sorted distances to a csv file
with open('sorted_distances.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Species ID', 'Total Distance'])
    for species_id, distance in sorted_distances:
        writer.writerow([species_id, distance])

print("Done!")  # Final progress update
