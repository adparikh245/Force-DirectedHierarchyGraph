# Clean and save the edges
with open('/Users/ananyaparikh/Documents/Coding/Force-DirectedHierarchyGraph/NewData/lastfm_asia_edges.txt', 'r') as f:
    edges = [line.strip() for line in f if line.strip()]

# Save cleaned edges back to a new file
with open('lastfm_asia_edges_cleaned.txt', 'w') as f:
    for edge in edges:
        f.write(f"{edge}\n")
print("Edge file created as 'c_lastfm_asia_edges.txt'")
