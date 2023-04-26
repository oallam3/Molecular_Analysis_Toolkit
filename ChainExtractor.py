import networkx as nx
import csv

# Load data from CSV file
with open('bonds.csv', 'r') as f:
    reader = csv.reader(f)
    data = [row for row in reader]

# Create graph
G = nx.Graph()
for row in data:
    G.add_edge(str(row[0]), str(row[1]))

# Get connected components
chains = list(nx.connected_components(G))

# Create dictionary
chain_dict = {}
chain_count = 1
mol_count = 1
for chain in chains:
    if len(chain) < 10:
        chain_name = 'molecule_' + str(mol_count)
        mol_count += 1
    else:
        chain_name = 'chain_' + str(chain_count)
        chain_count += 1
    chain_dict[chain_name] = [int(node) for node in chain]

# Write to new CSV file
with open('output.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Chain'] + ['Numbers'])
    for key, value in chain_dict.items():
        writer.writerow([key] + value)

