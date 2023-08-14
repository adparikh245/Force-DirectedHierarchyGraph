// import { readFileSync, writeFile, readFile } from 'fs'

export function readTxt(url){
  fetch(url)
    .then((response) => response.text())
    .then((fileContent) => {
      return fileContent
    })
    .catch((error) => {
      console.log(
        'Error occurred while reading the file:',
        error
      );
    });
};

export function parseGraphDataset(dataset) {
  // Define an empty graph object
  const graph = {
    nodes: [],
    links: [],
  };
  

  // Split the dataset into individual lines
  const lines = dataset.split('\n');

  // Iterate over each line and parse the node and weight information
  lines.forEach((line) => {
    const [node1, node2, weight] = line
      .split(' ')
      .map(Number);
    if (
      typeof node1 !== 'undefined' &&
      typeof node2 !== 'undefined'
    ) {
      const node1Index = graph.nodes.findIndex(
        (node) => node.id === node1
      );
      if (node1Index === -1) {
        graph.nodes.push({ id: node1, group: 1 });
      }
      const node2Index = graph.nodes.findIndex(
        (node) => node.id === node2
      );
      if (node2Index === -1) {
        graph.nodes.push({ id: node2, group: 2 });
      }

      graph.links.push({
        source: node1,
        target: node2,
        value: weight,
      });
    }
  });

  return graph;
}

export function parseGraphDatasetLabel(dataset) {
  // Define an empty graph object
  const mapIdToLabel = {};

  // Split the dataset into individual lines
  const lines = dataset.split('\n');

  // Iterate over each line and parse the node and weight information
  lines.forEach((line) => {
    const [id, label] = line
      .split(' ')
      .map(Number);
    if (
      typeof id !== 'undefined' &&
      typeof label !== 'undefined'
    ) {
      if (!(id in mapIdToLabel)) {
        mapIdToLabel[id] = label;
      }
    }
  });

  return mapIdToLabel;
}

export function parseGraphDatasetSNAP(
  dataset,
  id2Name=null,
  mapIdToLabel=null,
) {
  const graph = {
    nodes: [],
    links: [],
  };
  const processedLinks = new Set(); // Create a Set to track processed link IDs

  const lines = dataset.split('\n');

  for (const line of lines) {
    const [node1, node2] = line.split(/[ \t]+/);
    if (/^\d+\s+\d+$/.test(line) || /^\d+\t+\d+$/.test(line)) {

      if (node1 !== '' && node2 !== '') {
        const node1ID = id2Name ? id2Name[node1.toString()] : node1.toString();
        const node2ID = id2Name ? id2Name[node2.toString()] : node2.toString();

        let node1Index = graph.nodes.findIndex((node) => node.id === node1ID);
        if (node1Index === -1) {
          const node1Data = { id: node1ID };
          if (mapIdToLabel && mapIdToLabel[node1ID]) {
            node1Data.group = mapIdToLabel[node1ID];
          }
          graph.nodes.push(node1Data);
        }

        let node2Index = graph.nodes.findIndex((node) => node.id === node2ID);
        if (node2Index === -1) { 
          const node2Data = { id: node2ID };
          if (mapIdToLabel && mapIdToLabel[node2ID]) {
            node2Data.group = mapIdToLabel[node2ID];
          }
          graph.nodes.push(node2Data);
        }

        const linkId1 = node1ID + '-' + node2ID;
        const linkId2 = node2ID + '-' + node1ID;

        if (!processedLinks.has(linkId1) && !processedLinks.has(linkId2)) {
          graph.links.push({ source: node1ID, target: node2ID });
          processedLinks.add(linkId1);
          processedLinks.add(linkId2);
        }
      }
    }
    }

  return graph;
}

export function parseGraphFromJSON(jsonText) {
  const jsonData = JSON.parse(jsonText);
  
  const graph = {
    nodes: [],
    links: [],
  };

  for (const node of jsonData.nodes) {
    graph.nodes.push({id: node.id, group: node.group});
  }

  for (const link of jsonData.links) {
    const source = link.source;
    const target = link.target;
    graph.links.push({ source, target });
  }

  return graph;
}

export function selectParseFormat(text, fileContent, id2Name=null){
  console.log(text);
  console.log(text.endsWith('.json'));
  let graph;
  if (text.endsWith('.json')) {
    // console.log('here');
    graph = parseGraphFromJSON(fileContent);
  } 
  else if (text.endsWith('.txt')) {
    	if (id2Name) {
      	graph = parseGraphDatasetSNAP(fileContent, id2Name);
      }
      else {
        // console.log('here');
        graph = parseGraphDatasetSNAP(fileContent);
      }
  } else {
    console.error("Error: Unsupported format");
  }
  console.log(graph);
  return graph;
}

export function swapKeysAndValues(jsonData) {
  try {
    const swappedData = {};

    for (const key in jsonData) {
      if (jsonData.hasOwnProperty(key)) {
        const value = jsonData[key];
        swappedData[value.toString()] = key;
      }
    }

    return swappedData;
  } catch (error) {
    console.error('Invalid JSON data:', error);
    return null;
  }
}

// Simple Random Walk-Based Node Embedding
export function nodeEmbedding(
  graph,
  numWalks,
  walkLength,
  dimensions,
  seed
) {
  const rng = seededRandom(seed); // Initialize the RNG with the seed.

  // Initialize embeddings for each node
  const embeddings = {};

  // Perform random walks
  for (let walk = 0; walk < numWalks; walk++) {
    graph.forEachNode((node) => {
      const walkPath = [node];

      // Generate random walk path using seeded RNG
      for (let step = 1; step < walkLength; step++) {
        const neighbors = graph.outNeighbors(node);
        if (neighbors.length === 0) break; // If no neighbors, end the walk.
        const nextNode = neighbors[Math.floor(rng() * neighbors.length)]; // Use the seeded RNG.
        walkPath.push(nextNode);
        node = nextNode;
      }

      // Update node embeddings
      for (const n of walkPath) {
        if (!embeddings[n]) {
          embeddings[n] = new Array(dimensions).fill(0);
        }
        embeddings[n][walk % dimensions] += 1;
      }
    });
  }

  // Normalize embeddings
  for (const node in embeddings) {
    const embedding = embeddings[node];
    const sum = embedding.reduce((acc, val) => acc + val, 0);
    embeddings[node] = embedding.map((val) => val / sum);
  }

  return embeddings;
}

export function kMeansClustering(
  embeddings,
  numClusters,
  maxIterations,
  seed
) {
  const nodes = Object.keys(embeddings);
  const embeddingsMatrix = nodes.map((node) => embeddings[node]);

  // Initialize centroids using seeded random number generator
  const rng = seededRandom(seed);
  let centroids = embeddingsMatrix.slice(0, numClusters).map(() => {
    const randomIndex = Math.floor(rng() * embeddingsMatrix.length);
    return embeddingsMatrix[randomIndex];
  });

  let clusters = [];

  // Iterate until convergence or maximum iterations reached
  for (let iteration = 0; iteration < maxIterations; iteration++) {
    // Assign nodes to clusters
    clusters = new Array(numClusters).fill().map(() => []);

    for (const [nodeIndex, embedding] of embeddingsMatrix.entries()) {
      let minDistance = Infinity;
      let assignedCluster = -1;

      for (const [centroidIndex, centroid] of centroids.entries()) {
        const distance = euclideanDistance(embedding, centroid);

        if (distance < minDistance) {
          minDistance = distance;
          assignedCluster = centroidIndex;
        }
      }

      clusters[assignedCluster].push(nodes[nodeIndex]);
    }

    // Update centroids
    const newCentroids = [];

    for (const cluster of clusters) {
      if (cluster.length === 0) {
        // Use seeded random number generator to pick a random centroid
        const randomIndex = Math.floor(rng() * embeddingsMatrix.length);
        newCentroids.push(embeddingsMatrix[randomIndex]);
      } else {
        const clusterEmbeddings = cluster.map((node) => embeddings[node]);
        const clusterCentroid = calculateCentroid(clusterEmbeddings);
        newCentroids.push(clusterCentroid);
      }
    }

    // Check convergence
    if (areCentroidsEqual(centroids, newCentroids)) {
      break;
    }

    centroids = newCentroids;
  }

  // Create a dictionary with nodes as keys and cluster counts as values
  const communityDict = {};
  const nodeClusterDict = {};
  clusters.forEach((cluster, index) => {
    let parentId = 'Com.'+index;
    communityDict[parentId] = cluster;
    cluster.forEach((node) => {
      nodeClusterDict[node] = index;
    });
  });
  
  return [clusters, communityDict, nodeClusterDict];
}

function findFarthestNode(embeddings, centroids) {
  let farthestNode = null;
  let maxDistance = -1;

  for (const node in embeddings) {
    const nodeEmbedding = embeddings[node];
    let minDistance = Infinity;

    for (const centroid of centroids) {
      const distance = euclideanDistance(nodeEmbedding, centroid);

      if (distance < minDistance) {
        minDistance = distance;
      }
    }

    if (minDistance > maxDistance) {
      maxDistance = minDistance;
      farthestNode = node;
    }
  }

  return farthestNode;
}

// Helper functions
function euclideanDistance(vector1, vector2) {
  let sum = 0;
  for (let i = 0; i < vector1.length; i++) {
    sum += Math.pow(vector1[i] - vector2[i], 2);
  }
  return Math.sqrt(sum);
}

function calculateCentroid(embeddings) {
  const numDimensions = embeddings[0].length;
  const centroid = new Array(numDimensions).fill(
    0
  );

  for (const embedding of embeddings) {
    for (let i = 0; i < numDimensions; i++) {
      centroid[i] += embedding[i];
    }
  }

  for (let i = 0; i < numDimensions; i++) {
    centroid[i] /= embeddings.length;
  }

  return centroid;
}

function areCentroidsEqual(
  centroids1,
  centroids2
) {
  for (let i = 0; i < centroids1.length; i++) {
    const centroid1 = centroids1[i];
    const centroid2 = centroids2[i];

    for (let j = 0; j < centroid1.length; j++) {
      if (centroid1[j] !== centroid2[j]) {
        return false;
      }
    }
  }

  return true;
}


// helper function
function cntBTWClust(clusters, edgeList) {
  const clusterEdgesCount = {};

  for (let i = 0; i < clusters.length; i++) {
    for (let j = i + 1; j < clusters.length; j++) {
      const clusterPair = `${i}-${j}`;
      clusterEdgesCount[clusterPair] = 0;
    }
  }

  for (const edge of edgeList) {
    const sourceClusterIndex = clusters.findIndex(cluster => cluster.includes(edge.source.toString()));
    const targetClusterIndex = clusters.findIndex(cluster => cluster.includes(edge.target.toString()));

    if (sourceClusterIndex !== targetClusterIndex) {
      const clusterPair = `${Math.min(sourceClusterIndex, targetClusterIndex)}-${Math.max(sourceClusterIndex, targetClusterIndex)}`;
      clusterEdgesCount[clusterPair]++;
    }
  }
  
  return clusterEdgesCount
};

class Graph {
  constructor() {
    this.nodes = new Map(); // Map to store nodes and their edges
  }

  addNode(nodeId) {
    if (!this.nodes.has(nodeId)) {
      this.nodes.set(nodeId, new Set());
    }
  }

  addEdge(sourceId, targetId) {
    // Add source node if not already present
    if (!this.nodes.has(sourceId)) {
      this.addNode(sourceId);
    }

    // Add target node if not already present
    if (!this.nodes.has(targetId)) {
      this.addNode(targetId);
    }

    // Add edge from source to target
    this.nodes.get(sourceId).add(targetId);
  }

  outNeighbors(nodeId) {
    if (this.nodes.has(nodeId)) {
      return [...this.nodes.get(nodeId)];
    } else {
      return [];
    }
  }

  forEachNode(callback) {
    for (let nodeId of this.nodes.keys()) {
      callback(nodeId);
    }
  }
}

export function constructGraph(nodeEdges) {
  // // Example usage
  const graph = new Graph();
  // ... Build and populate the graph with nodes and edges
  // Add nodes to the graph
  nodeEdges.nodes.forEach((node) => {
    graph.addNode(node.id);
  });

  // Add edges to the graph
  nodeEdges.links.forEach((link) => {
    graph.addEdge(link.source, link.target);
  });
  return graph;
}

export function constGraphFromClust(clusters, edgeList, colorPalette){
  const groupColors = {};
  let clusterEdgesCount = cntBTWClust(clusters, edgeList);
  let commNodeEdges = {nodes: [], links: []};
	for (let i = 0; i < clusters.length; i++) {
    const colorIndex = i % colorPalette.length;
    const color = colorPalette[colorIndex];
    groupColors[i] = color;
    
    // Create nodes of the first community layer
  	let currentId = 'Com.' + i;
    commNodeEdges.nodes.push({id: currentId, group: i, groupColor: groupColors[i], numEle:clusters[i].length})};
  	
  
  // Create edges of the first community layer
  Object.entries(clusterEdgesCount).forEach(([key, value]) => {
    if(value > 0){
      let commPair = key.split('-');
      let currentSrce = 'Com.' + commPair[0];
      let currentTGT  = 'Com.' + commPair[1];
      commNodeEdges.links.push({source: currentSrce, target: currentTGT, value});
  	}
  });
  return [commNodeEdges, groupColors];
};


export function constSubGFromClust(clust, edgeList, numCluster, groupColors, nodeClusterDict, config){
  const maxNumNodePerCluster = config[0];
  const maxNumEdgePerCluster = config[1];
  const maxNumAdjNodePerCluster = config[2];
  const maxNumAdjEdgePerCluster = config[3];
  
  // Generate edges
  const clusterSet = new Set(clust);
  let data = {};
  
  // Generate nodes and links
  data.nodes = Array.from(clusterSet, (clus) => ({ id: clus, group: numCluster, groupColor: groupColors[numCluster], numEle:0}));
  data.links = edgeList.filter((edge) => clusterSet.has(edge.source) && clusterSet.has(edge.target));
  if (data.nodes.length > maxNumNodePerCluster) {
  	data = filterDataByTopNodes(data, maxNumNodePerCluster);
  }
  
  
  data.links = data.links.slice(0, maxNumEdgePerCluster);

  let numAdjEdges = 0;
  const existingNodes = new Set(data.nodes.map((node) => node.id)); // Create a Set of existing node IDs

  edgeList.forEach((edge) => {
    // Edges of the current cluster
    const isEdgeInDataLinks = data.links.some((existingEdge) => (
      (existingEdge.source === edge.source && existingEdge.target === edge.target) ||
      (existingEdge.source === edge.target && existingEdge.target === edge.source)
    ));

    const isSourceInDataNodes = existingNodes.has(edge.source); // Check if source node already exists
    const isTargetInDataNodes = existingNodes.has(edge.target); // Check if target node already exists

    // Edges of adjacent clusters
    if (!isEdgeInDataLinks && numAdjEdges <= maxNumAdjEdgePerCluster && (clusterSet.has(edge.source) || clusterSet.has(edge.target))) {
      data.links.push(edge);
      numAdjEdges++;

      if (!isSourceInDataNodes) {
        let numClusterS = nodeClusterDict[edge.source];
        data.nodes.push({
          id: edge.source,
          group: numClusterS,
          groupColor: groupColors[numClusterS],
          numEle:0
        });
        existingNodes.add(edge.source); // Add source node ID to the existingNodes set
      }
      if (!isTargetInDataNodes) {
        let numClusterT = nodeClusterDict[edge.target];
        data.nodes.push({
          id: edge.target,
          group: numClusterT,
          groupColor: groupColors[numClusterT],
          numEle:0
        });
        existingNodes.add(edge.target); // Add target node ID to the existingNodes set
      }
    }
  });
  return data;
};



export function filterDataByTopNodes(data, numTopNodes) {
  const edgeCount = {};
  data.links.forEach(link => {
    const sourceNodeID = link.source;
    const targetNodeID = link.target;

    edgeCount[sourceNodeID] = (edgeCount[sourceNodeID] || 0) + 1;
    edgeCount[targetNodeID] = (edgeCount[targetNodeID] || 0) + 1;
  });

  const topNodes = Object.entries(edgeCount)
    .sort((a, b) => b[1] - a[1])
    .slice(0, numTopNodes)
    .map(([nodeID]) => nodeID);

  const filteredNodes = data.nodes.filter(node => topNodes.includes(node.id));
  const filteredEdges = data.links.filter(link =>
    topNodes.includes(link.source) && topNodes.includes(link.target)
  );
  
  return  {nodes: filteredNodes, links: filteredEdges,}
}

// -------------  Utility functions

let levelsData = {};
// Update data for a specific level
export function updateLevelData(level, data) {
  levelsData[level] = data;
  console.log(levelsData);
}

// Access data for a specific level
export function getLevelData(level) {
  return levelsData[level];
}

// Custom deterministic random number generator
function seededRandom(seed) {
  let m_w = seed;
  let m_z = 987654321;

  // Returns a random number between 0 and 1
  return function () {
    m_z = (36969 * (m_z & 65535) + (m_z >> 16)) & 0xffffffff;
    m_w = (18000 * (m_w & 65535) + (m_w >> 16)) & 0xffffffff;
    let result = ((m_z << 16) + m_w) & 0xffffffff;
    result /= 4294967296;
    return result + 0.5;
  };
}
