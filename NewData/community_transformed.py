import json

# Load target data to form communities
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_target.json', 'r', encoding='utf-8') as f:
    target_data = json.load(f)

# Convert targets to community format
community_data = {}

# Assume target_data is a list of dictionaries like [{ "node_id": "community_id" }, ...]
for entry in target_data:
    # Assuming each entry has keys "node_id" and "community_id"
    node_id = entry.get("node_id")
    community_id = entry.get("community_id")
    
    # Ensure we have both values
    if node_id is not None and community_id is not None:
        community_key = f"Com.{community_id}"
        if community_key not in community_data:
            community_data[community_key] = []
        community_data[community_key].append(node_id)

# Save community data
with open('c_lastfm_asia_community.json', 'w') as f:
    json.dump(community_data, f, indent=4)
print("Community file created as 'c_lastfm_asia_community.json'")
