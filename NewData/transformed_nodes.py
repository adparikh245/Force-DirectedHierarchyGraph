import json

# Load edges to gather unique node IDs
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_edges.txt', 'r') as f:
    edges = [line.strip().split() for line in f.readlines()]

# Extract unique node IDs
node_ids = set()
for edge in edges:
    for node in edge:
        clean_node = ''.join(filter(str.isdigit, node))
        if clean_node:
            node_ids.add(clean_node)

# Create dummy names for nodes
node_data = {str(node_id): f"User {node_id}" for node_id in sorted(node_ids, key=int)}

# Save to JSON file
with open('c_lastfm_asia_user.json', 'w') as f:
    json.dump(node_data, f, indent=4)
print("Node file created as 'c_lastfm_asia_user.json'")
