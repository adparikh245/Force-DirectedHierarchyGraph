import json

# Load the id2name mapping from the JSON file
id2name_path = '/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/id2name.json'
with open(id2name_path, 'r', encoding='utf-8') as f:
    id2name_mapping = json.load(f)  # This loads names as keys and IDs as values

# Reverse the id2name mapping to get a name-to-id dictionary
name_to_id_mapping = {value: key for key, value in id2name_mapping.items()}

# Load target data to form communities
edges_path = '/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_edges_cleaned.txt'
with open(edges_path, 'r', encoding='utf-8') as f:
    edges = [line.strip().split() for line in f.readlines()]  # Process as a list of edges

# Initialize a dictionary to store communities
community_data = {}

# Process edges and group nodes into communities
for edge in edges:
    if len(edge) == 2:  # Ensure the edge contains exactly two nodes
        # Convert both node IDs to integers
        node_id_1, node_id_2 = map(int, edge)

        # Define a community key
        community_key = f"Com.{min(node_id_1, node_id_2)}"

        # Initialize the community key if it doesn't exist
        if community_key not in community_data:
            community_data[community_key] = []

        # Map node IDs to user names and add them to the community
        for node_id in (node_id_1, node_id_2):
            user_name = name_to_id_mapping.get(node_id, f"User {node_id}")  # Retrieve user name by ID
            if user_name not in community_data[community_key]:  # Avoid duplicates
                community_data[community_key].append(user_name)

# Save the resulting community data to a JSON file
output_path = '/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/c_lastfm_asia_community.json'

with open(output_path, 'w', encoding='utf-8') as f:
    # Write the JSON file with indentation for readability
    json.dump(community_data, f, indent=4, ensure_ascii=False)
