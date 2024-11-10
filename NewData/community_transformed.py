import json

# Load target data to form communities
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_target.json', 'r', encoding='utf-8') as f:
    target_data = json.load(f)

# Convert targets to community format
community_data = {}

for entry in target_data:
    node_id = entry.get("id")
    community_id = entry.get("target")
    
    if node_id is not None and community_id is not None:
        community_key = f"Com.{community_id}"
        if community_key not in community_data:
            community_data[community_key] = []
        community_data[community_key].append(node_id)

# Write to file in chunks
with open('c_lastfm_asia_community.json', 'w') as f:
    f.write("{\n")
    for i, (community_key, nodes) in enumerate(community_data.items()):
        json.dump({community_key: nodes}, f)
        if i < len(community_data) - 1:
            f.write(",\n")
    f.write("\n}")

print("Community file created as 'c_lastfm_asia_community.json'")
