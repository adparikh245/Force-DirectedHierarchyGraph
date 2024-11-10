import json

# Load edges to gather unique node IDs
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_edges.txt', 'r') as f:
    edges = [line.strip().split() for line in f.readlines()]

# Extract unique node IDs from edges and clean them
node_ids = set()
for edge in edges:
    for node in edge:
        # Remove any non-numeric characters from node ID and ensure it's valid
        clean_node = ''.join(filter(str.isdigit, node))
        if clean_node:  # Only add if it has digits
            node_ids.add(clean_node)

# Create node mapping with dummy names
node_data = {str(node_id): f"User {node_id}" for node_id in sorted(node_ids, key=int)}

# Save as JSON in the correct format
with open('c_lastfm_asia_user.json', 'w') as f:
    json.dump(node_data, f, indent=4)

# Verify edge formatting and clean edges
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_edges.txt', 'r') as f:
    edges = [line.strip() for line in f if line.strip()]

# Save cleaned edges back to file
with open('lastfm_asia_edges_cleaned.txt', 'w') as f:
    for edge in edges:
        f.write(f"{edge}\n")

# Load the lastfm_asia_target.json file
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_target.json', 'r') as f:
    target_data = json.load(f)

# Convert targets to community format
community_data = {}
for node_id, community_id in target_data.items():
    community_key = f"Com.{community_id}"
    if community_key not in community_data:
        community_data[community_key] = []
    community_data[community_key].append(node_id)

# Save the new community structure
with open('c_lastfm_asia_community.json', 'w') as f:
    json.dump(community_data, f, indent=4)
