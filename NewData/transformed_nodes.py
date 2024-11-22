import json

# Load edges to gather unique node IDs
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_edges_cleaned.txt', 'r') as f:
    edges = [line.strip().split() for line in f.readlines()]

# Extract unique node IDs
node_ids = set()
for edge in edges:
    for node in edge:
        clean_node = ''.join(filter(str.isdigit, node))  # Extract digits from node
        if clean_node:
            node_ids.add(int(clean_node))  # Convert to integer

# Create the id-to-name mapping with names as keys and IDs as values
id2name_mapping = {f"User {node_id}": node_id for node_id in sorted(node_ids)}

# Save the mapping to a JSON file
output_path = '/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/id2name.json'
with open(output_path, 'w', encoding='utf-8') as f:
    json.dump(id2name_mapping, f, indent=4, ensure_ascii=False)

print(f"ID-to-name mapping saved successfully to: {output_path}")
