import json

# Load edges to gather unique node IDs
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_edges.txt', 'r') as f:
    edges = [line.strip().split() for line in f.readlines()]

# Extract unique node IDs
node_ids = set()
for edge in edges:
    for node in edge:
        clean_node = ''.join(filter(str.isdigit, node))  # Extract digits from node
        if clean_node:
            node_ids.add(int(clean_node))  # Convert to integer

# Create dummy names for nodes
node_data = {node_id: f"User {node_id}" for node_id in sorted(node_ids)}

# Save the Python dictionary to a new Python file
output_file = '/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/node_data.py'
with open(output_file, 'w') as f:
    f.write("# This file contains the node_data dictionary\n")
    f.write("node_data = {\n")
    for node_id, user_name in node_data.items():
        f.write(f"    {node_id}: '{user_name}',\n")
    f.write("}\n")
