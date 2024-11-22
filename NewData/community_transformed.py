import json

# Load target data to form communities
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_target.json', 'r', encoding='utf-8') as f:
    target_data = json.load(f)

# Load the user mapping (node_id to user name)
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/c_lastfm_asia_user.json', 'r', encoding='utf-8') as f:
    user_mapping = json.load(f)

# Convert targets to community format with user names
community_data = {}

for entry in target_data:
    node_id = entry.get("id")
    community_id = entry.get("target")

    if node_id is not None and community_id is not None:
        community_key = f"Com.{community_id}"
        if community_key not in community_data:
            community_data[community_key] = []
        # Map the node_id to the user name
        user_name = user_mapping.get(str(node_id), f"User {node_id}")
        community_data[community_key].append(user_name)

# Save community data to JSON in a compact format
output_path = '/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/c_lastfm_asia_community.json'

with open(output_path, 'w', encoding='utf-8') as f:
    # Write the JSON with proper formatting
    json.dump(community_data, f, indent=4, ensure_ascii=False)

print(f"Community file created successfully at: {output_path}")
