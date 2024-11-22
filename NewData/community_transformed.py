import json
from node_data import node_data  # Import node_data directly as a Python dictionary

# Load target data to form communities (assumes edges file is plain text, not JSON)
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_edges_cleaned.txt', 'r', encoding='utf-8') as f:
    target_data = [line.strip().split() for line in f.readlines()]  # Process as a list of edges

# Convert edges to community format with user names
community_data = {}

for edge in target_data:
    if len(edge) == 2:  # Ensure the edge contains two nodes
        node_id_1, node_id_2 = map(int, edge)  # Convert both IDs to integers

        # Example logic: Assign both nodes to the same community
        community_key = f"Com.{min(node_id_1, node_id_2)}"  # Use the smaller ID as the community key

        if community_key not in community_data:
            community_data[community_key] = []

        # Map node IDs to user names and add them to the community
        for node_id in (node_id_1, node_id_2):
            user_name = node_data.get(node_id, f"User {node_id}")
            if user_name not in community_data[community_key]:  # Avoid duplicates
                community_data[community_key].append(user_name)

# Save community data to JSON in a compact format
output_path = '/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/c_lastfm_asia_community.json'

with open(output_path, 'w', encoding='utf-8') as f:
    # Write the JSON with proper formatting
    json.dump(community_data, f, indent=4, ensure_ascii=False)

print(f"Community file created successfully at: {output_path}")

