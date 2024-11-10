import json

# Load the original file
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_features.json', 'r') as f:
    original_data = json.load(f)

# Transform the data to an array of node objects
new_data = []
for node_id, features in original_data.items():
    # Create a new node object for each entry
    node_object = {
        "id": int(node_id),
        "features": features  # Add the list of features as a property
    }
    new_data.append(node_object)

# Save the new data format to a new file
with open('transformed_lastfm_asia_nodes.json', 'w') as f:
    json.dump(new_data, f, indent=4)

print("Transformation complete. Check 'transformed_lastfm_asia_nodes.json' for the new format.")
