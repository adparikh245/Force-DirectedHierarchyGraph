import json

# Load the original file
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_features.json', 'r') as f:
    original_data = json.load(f)

# Transform the data to the desired format
new_data = {}
for node_id, features in original_data.items():
    # Create a new key with the "Com." prefix and add the list of features
    community_key = f"Com.{node_id}"
    new_data[community_key] = features

# Save the new data format to a new file
with open('transformed_lastfm_asia_features.json', 'w') as f:
    json.dump(new_data, f, indent=4)

print("Transformation complete. Check 'transformed_lastfm_asia_features.json' for the new format.")
