import json
import random

# Load the existing community file
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/c_lastfm_asia_community.json', 'r', encoding='utf-8') as f:
    original_community_data = json.load(f)

# Combine all users into a single list
all_users = []
for users in original_community_data.values():
    all_users.extend(users)

# Shuffle the list of users for redistribution
random.shuffle(all_users)

# Set the number of clusters and maximum nodes per community
num_clusters = 14
max_nodes_per_community = 150

# Initialize the reduced community structure
reduced_community_data = {f"Com.{i}": [] for i in range(num_clusters)}

# Redistribute users into clusters with a maximum node limit
for user in all_users:
    # Find the first cluster with available space
    for i in range(num_clusters):
        if len(reduced_community_data[f"Com.{i}"]) < max_nodes_per_community:
            reduced_community_data[f"Com.{i}"].append(user)
            break

# Save the reduced community file
output_path = '/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/c_lastfm_asia_community_reduced.json'
with open(output_path, 'w', encoding='utf-8') as f:
    json.dump(reduced_community_data, f, indent=4, ensure_ascii=False)