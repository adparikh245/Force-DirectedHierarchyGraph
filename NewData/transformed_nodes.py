import json

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

# Create dummy names for nodes
node_data = {node_id: f"User {node_id}" for node_id in sorted(node_ids)}

# Save to file with integer keys
output_file = '/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/node_user_mapping.txt'
with open(output_file, 'w') as f:
    f.write("{\n")
    for i, (key, value) in enumerate(node_data.items()):
        f.write(f"  {key}: \"{value}\"")  # Write key as an integer
        if i < len(node_data) - 1:
            f.write(",\n")
    f.write("\n}")

