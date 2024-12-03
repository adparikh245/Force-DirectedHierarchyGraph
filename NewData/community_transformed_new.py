import json
import random

# Load edges to gather unique node IDs
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_edges_cleaned.txt', 'r') as f:
    edges = [line.strip().split() for line in f.readlines()]

# Extract unique node IDs
node_ids = set()
for edge in edges:
    for node in edge:
        clean_node = ''.join(filter(str.isdigit, node))  # Extract digits from node
        if clean_node:
            node_ids.add(int(clean_node))  # Convert to integer

# Define the number of communities
num_clusters = 15  # Match with numClusters in the config file

# Randomly assign users to communities
community_data = {f"Com.{i}": [] for i in range(num_clusters)}
for node_id in sorted(node_ids):
    # Assign each user to a random community
    random_community = f"Com.{random.randint(0, num_clusters - 1)}"
    community_data[random_community].append(node_id)

# Save the reduced community data to a JSON file
output_path = '/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/c_lastfm_asia_community_reduced.json'
with open(output_path, 'w', encoding='utf-8') as f:
    json.dump(community_data, f, indent=4, ensure_ascii=False)
