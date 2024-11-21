# Clean and save the edges
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_edges.txt', 'r') as f:
    edges = [line.strip() for i, line in enumerate(f) if i > 0 and line.strip()]  # Skip the first line

# Save cleaned edges back to a new file
with open('lastfm_asia_edges_cleaned.txt', 'w') as f:
    for edge in edges:
        f.write(f"{edge}\n")
print("Edge file created as 'lastfm_asia_edges_cleaned.txt'")
