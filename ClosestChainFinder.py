import csv
import math

# Set up the parameters
target_species = # Atom ID of target (diffusing) atom for analysis
type_ids = # Vector containing atom ID types corresponding to forcefields of the atoms on the polymer chains
           # that typically coordinate with the diffusing atom.
num_neighbors = # Specify how many of the closest atoms (of type ID specified above) to the target atom will be considered
num_atoms = # Number of atoms in your system.
trj_file = # Name of your LAMMPS trajectory file.

# Function to read the .lammpstrj file
def read_lammpstrj(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    rows = []
    for line in lines[9:]:
        row = line.split()
        if row and row[0].isdigit():  # Check if the first item is a number (assumes this is an atom entry)
            numeric_row = [float(val) for val in row if is_float(val)]
            rows.append(numeric_row)
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

# Initialize the dictionary to hold the closest chains for each loop
closest_chains = {}
row_num = 0
x = 1 # Analyze every xth timestep in your simulation.
y = 201 # Max number of timesteps to analyze
atom_num = num_atoms+9

for loop_num in range(0, y, x):
    print(f"Processing loop {loop_num}")  # Progress update
    print(f"Processing from row {row_num}")  # Progress update

    # Find the row for the target species
    target_row = None
    search_start = max(0, row_num - atom_num)
    search_end = min(row_num + 2 * atom_num, len(rows))
    for row in rows[search_start:search_end]:
        if row[0] == target_species:
            target_row = row
            break

    if target_row is None:
        print(f"Target species {target_species} not found in loop {loop_num} (row {row_num}-{row_num+atom_num-1})")
        continue


    # Find the coordinates of the target species
    target_coords = tuple(target_row[3:6])

    # Find the distances to all other species with the specified type ids
    distances = []
    for row in rows[row_num:row_num+atom_num]:
        if len(row) < 2:  # Skip the row if it has less than 2 elements
            continue
        if row[1] in type_ids and row[0] != target_species:
            coords = tuple(row[3:6])
            distance = math.sqrt(sum([(a - b)**2 for a, b in zip(target_coords, coords)]))
            distances.append((row[0], distance))

    # Sort the distances and keep only the closest num_neighbors
    distances = sorted(distances, key=lambda x: x[1])
    closest = distances[:num_neighbors]

    # Find the chains that contain the closest species
    closest_chains_this_loop = []
    for species_id, distance in closest:
        for chain, species_list in chains.items():
            if species_id in species_list:
                closest_chains_this_loop.append(chain)
                break

    # Add this loop's closest chains to the dictionary
    closest_chains[f'loop_{loop_num}'] = closest_chains_this_loop

    # Increment the row number to skip to the next loop
    row_num = row_num+atom_num*x
    

# Write the closest chains to a csv file
with open('closest_chains.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Loop'] + [f'Closest {i+1}' for i in range(num_neighbors)])
    for loop, chains in closest_chains.items():
        row = [loop] + chains + ['molecule_' + str(i+1) for i in range(num_neighbors - len(chains))]
        writer.writerow(row)
        
print("Done!")  # Final progress update

##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################
##################################################################################################################################################


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_color(loop_label):
    chain = closest_chains[loop_label][0]
    if chain in color_mapping:
        return np.array([plt.cm.get_cmap("tab10")(color_mapping[chain])])
    else:
        return np.array([0, 0, 0])

# Create a color mapping for the chains
color_mapping = {chain: i for i, chain in enumerate(sorted(list(chains)))}
print(chains)

# Extract the target species positions
# Extract the target species positions
target_positions = []
for i in range(len(rows)):
    if i % (atom_num+9 * x) == 0 and i // (atom_num+9 * x) < y:
        target_positions.append(rows[i][3:6])

# Create a list of loop numbers for which we have data in closest_chains
loop_numbers = [int(loop_name.split("_")[1]) for loop_name in closest_chains.keys()]

# Plot the target species positions
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")


# Plot the target species positions with colors based on the closest chain
for i, (loop_num, pos) in enumerate(zip(loop_numbers, target_positions)):
    label = f"loop_{loop_num}"
    ax.scatter(*pos, c=get_color(label), marker="o", s=50, edgecolors="k", linewidths=0.5, depthshade=True)
    
    # Connect the points with a black line
    if i > 0:
        ax.plot([prev_pos[0], pos[0]], [prev_pos[1], pos[1]], [prev_pos[2], pos[2]], c='k')
        
    prev_pos = pos

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.grid(False)  # removes the grid lines
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
plt.show()










