(function (d3) {
  'use strict';

  // import { readFileSync, writeFile, readFile } from 'fs'

  function parseGraphDatasetSNAP(
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

  function parseGraphFromJSON(jsonText) {
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

  function selectParseFormat(text, fileContent, id2Name=null){
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

  function swapKeysAndValues(jsonData) {
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
  function nodeEmbedding(
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

  function kMeansClustering(
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
  }
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

  function constructGraph(nodeEdges) {
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

  function constGraphFromClust(clusters, edgeList, colorPalette){
    const groupColors = {};
    let clusterEdgesCount = cntBTWClust(clusters, edgeList);
    let commNodeEdges = {nodes: [], links: []};
  	for (let i = 0; i < clusters.length; i++) {
      const colorIndex = i % colorPalette.length;
      const color = colorPalette[colorIndex];
      groupColors[i] = color;
      
      // Create nodes of the first community layer
    	let currentId = 'Com.' + i;
      commNodeEdges.nodes.push({id: currentId, group: i, groupColor: groupColors[i], numEle:clusters[i].length});}  	
    
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
  }

  function constSubGFromClust(clust, edgeList, numCluster, groupColors, nodeClusterDict, config){
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
  }


  function filterDataByTopNodes(data, numTopNodes) {
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
  function updateLevelData(level, data) {
    levelsData[level] = data;
    console.log(levelsData);
  }

  // Access data for a specific level
  function getLevelData(level) {
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

  // import { readFileSync, writeFile, readFile } from 'fs'


  const width = window.innerWidth;
  const height = window.innerHeight;

  // Create an SVG container for the graph
  const svg = d3.select('body')
    .append('svg')
    .attr('width', width)
    .attr('height', height);

  // Create a <g> element as the container for the graph elements
  var container = svg.append('g');

  function zoomed(event) {
    container.attr('transform', event.transform);
  }

  // Create a zoom behavior
  const zoom = d3.zoom()
    .scaleExtent([0.1, 4]) // Set the zoom scale range
    .on('zoom', zoomed); // Specify the zoom event handler

  // Apply the zoom behavior to the SVG container
  svg.call(zoom);

  ////////// Step 1: Read data from a URL 

  async function fetchData() {
    try {
      
      // Fetch config.json
      const configT = await fetch('config.json');
      const configText = await configT.text();
      const config = JSON.parse(configText);
      
      const id2NameUrl = config.id2nameURL;
      const nodesEdgesUrl = config.edgesURL;
      const commUrl = config.communityFile;
      // id2NameUrl = 'https://raw.githubusercontent.com/Khaninsi/D3-visualization/main/data/user.json';
      // let nodesEdgesUrl = 'https://raw.githubusercontent.com/Khaninsi/D3-visualization/main/data/email-Eu-core.txt';
      // let nodesEdgesUrl = 'https://raw.githubusercontent.com/Khaninsi/D3-visualization/main/data/enron_edges.txt';

      let nodeLabelObj;
      let nodeEdges;
      let complexGraph;
      let id2NameDict;
      let communityDict;;
      
      // Fetch nodesEdgesUrl
      const response1 = await fetch(nodesEdgesUrl);
      const fileContent1 = await response1.text();
      
      // Fetch id2NameUrl if given
      if (id2NameUrl && id2NameUrl.length>0) {
        const response2 = await fetch(id2NameUrl);
        const response2Text = await response2.text();
        const response2JS = JSON.parse(response2Text);
        id2NameDict =  swapKeysAndValues(response2JS);
        nodeEdges = selectParseFormat(nodesEdgesUrl, fileContent1, id2NameDict);
      }
      else {
  			nodeEdges = selectParseFormat(nodesEdgesUrl, fileContent1);
      }

      if (commUrl.length > 1){
        const response3 = await fetch(commUrl);
        const response3Text = await response3.text();
        communityDict = JSON.parse(response3Text);
        // console.log('communityDict');
        // console.log(communityDict);
      }
      
      
      
      // Continue with the next section or perform additional operations
      console.log('Completed fetching and processing data.');
      return {
        nodeEdges,
        config,
        communityDict,
      };
    } catch (error) {
      console.log('Error occurred while reading the file:', error);
    }
  }

  // Call the fetchData function to initiate the execution
  fetchData()
  	.then(({nodeEdges, config, communityDict}) => {
      // Example: Perform operations using nodeEdges
      console.log('First node edge:', nodeEdges.nodes[0]);
    
    	// Limit number of the total nodes by the configuration if exceeds
    	if (nodeEdges.nodes.length > config.maxNumNodes){
      	nodeEdges = filterDataByTopNodes(nodeEdges, config.maxNumNodes);
        console.log('Select topn nodes');
      }
    	const graphObj = constructGraph(nodeEdges);

      // TODO: determine these 3 parameters dynamically
      const numWalks = config.numWalks; 
      const walkLength = config.walkLength;
      const dimensions = config.dimensions;
    	const numClusters = config.numClusters;
    	const maxIterations = config.maxIterations;
      const colorPalette = ['#f94144', '#f3722c', '#f9844a', 
      '#f9c74f', '#90be6d', '#43aa8b', '#4d908e', '#577590', '#277da1', 
      '#0f8b8d', '#a8201a', '#06d6a0', '#585123', '#eec170', 
      '#f2a65a', '#f58549', '#772f1a', '#4c5760', '#93a8ac', '#a59e8c',
      '#66635b', '#002500', '#0fa3b1', '#eddea4', '#f7a072', '#ff9b42', '#fe5f55', '#f0b67f'];
    
  		// TODO: Improve this performance
      const embeddings = nodeEmbedding(graphObj, numWalks, walkLength, dimensions, 0);
    
    	let clusters;
    	// let communityDict;
    	let nodeClusterDict;
    
    	if (config.communityFile.length > 1){
        	nodeClusterDict = {};
      		// Read clusters and communityDict from the file instead
        	clusters = Object.values(communityDict);
          nodeClusterDict = {};
          clusters.forEach((cluster, index) => {
            cluster.forEach((node) => {
              nodeClusterDict[node] = index;
            });
    			});
      }
    	else {
        const clustersOutput = kMeansClustering(embeddings, numClusters, maxIterations, 0);
        clusters = clustersOutput[0];
        communityDict = clustersOutput[1];
        nodeClusterDict = clustersOutput[2];
      }    	
    	// Create nodes and edges of community from the clusters
    	const output = constGraphFromClust(clusters, nodeEdges.links, colorPalette);
    	const commNodeEdges = output[0];
    	const groupColors   = output[1];
    
    	const graphRelation = [];
    
  		const configNumEdgesConst = [config.numNodePerCluster, config.numEdgePerCluster, config.numAdjNodePerCluster, config.numAdjEdgePerCluster];
      
      let i = 0;
      for (const [key, value] of Object.entries(communityDict)) {
          graphRelation.push({parentId: key, children:constSubGFromClust(value, nodeEdges.links, 
                                                                        i, groupColors, nodeClusterDict, 
                                                                        configNumEdgesConst)});
        i += 1;
  		}    
      
    	// Create nodes and edges of second layers from the clusters
    	// for(let i = 0; i < clusters.length; i++) {
    	// let clus = clusters[i];
    	// let parentId = 'Com.'+i;
    	// graphRelation.push({parentId: parentId, children:constSubGFromClust(clus, nodeEdges.links, 
    	// i, groupColors, nodeClusterDict, 
    	// configNumEdgesConst)})
    	// };
    	console.log(graphRelation);

      const graph = commNodeEdges;

      // Implement drillUp function
      let currentLevel = 0;
      let currentNode;

      function drillUp() {
        container.selectAll('*').remove();
        if (currentLevel > 0) {
          // Get the previous graph data based on the current level
          currentLevel--;
          currentNode = getLevelData(currentLevel);
          renderGraph(currentNode);
        }
      }

      function reset() {
        container.selectAll('*').remove();
        currentLevel = 0;
        currentNode = getLevelData(currentLevel);
        renderGraph(currentNode);
      }

      // Add drill-up button
      const button = svg
        .append('g')
        .attr('class', 'button')
        .attr(
          'transform',
          'translate(' + (width - 100) + ', 20)'
        ) // Adjust the position as needed
        .on('click', drillUp);

      button
        .append('rect')
        .attr('width', 80)
        .attr('height', 30)
        .attr('rx', 5)
        .attr('ry', 5)
        .style('fill', 'lightblue');

      button
        .append('text')
        .attr('x', 40)
        .attr('y', 20)
        .attr('text-anchor', 'middle')
        .attr('alignment-baseline', 'middle')
        .text('<- Back');

      // Add another button below the first button
      const button2 = svg
        .append('g')
        .attr('class', 'button')
        .attr(
          'transform',
          'translate(' + (width - 100) + ', 60)'
        )
        .on('click', reset);

      button2
        .append('rect')
        .attr('width', 80)
        .attr('height', 30)
        .attr('rx', 5)
        .attr('ry', 5)
        .style('fill', 'lightgreen');

      button2
        .append('text')
        .attr('x', 40)
        .attr('y', 20)
        .attr('text-anchor', 'middle')
        .attr('alignment-baseline', 'middle')
        .text('Reset');

      updateLevelData(currentLevel, graph);
      renderGraph(graph);

      function renderGraph(graph) {
        container.selectAll('.link').remove();
        container.selectAll('.node').remove();
        container.selectAll('*').remove();

        const nodes = graph.nodes;
        const links = graph.links;

        const x = d3.scaleOrdinal().range([
          20,
          width - 20,
        ]);

        // Create the force simulation
        const simulation = d3.forceSimulation(nodes)
          .force(
            'link',
            d3.forceLink(links).id((d) => d.id)
          )
          .force(
            'charge',
            d3.forceManyBody().strength(-200)
          )
          .force(
            'center',
            d3.forceCenter(width / 2, height / 2)
          );

        // Initialize the positions of the nodes
        nodes.forEach((node) => {
          node.x =
            Math.random() * (width - 20) + 10;
          node.y =
            Math.random() * (height - 20) + 10;
        });

        const strokeWidthScale = d3.scaleLinear()
          .domain([0, d3.max(links, (d) => d.value)]) // Define the input domain
          .range([0.1, 3]); // Define the output range for stroke-width

        const strokeOpacity = d3.scaleLinear()
          .domain([0, d3.max(links, (d) => d.value)]) // Define the input domain
          .range([0.1, 0.4]); // Define the output range for stroke-width

        // Create the links
        const link = container
          .selectAll('.link')
          .data(links)
          .join('line')
          .attr('class', 'link')
          .style('opacity', 0.4)
          .style('stroke', '#aaa')
          .style('stroke-width', (d) =>
            strokeWidthScale(d.value || 0.1)
          );

        const groupColors = {};
        // Make a hash table for each group
        const colorPalette = ['#f94144', '#f3722c', '#f8961e', '#f9844a', 
        '#f9c74f', '#90be6d', '#43aa8b', '#4d908e', '#577590', '#0f8b8d', '#a8201a', 
        '#ef476f', '#ffd166', '#06d6a0', '#073b4c', '#585123', '#eec170', '#f2a65a', 
        '#f58549', '#772f1a', '#4c5760', '#93a8ac', '#a59e8c',
        '#66635b', '#002500', '#0fa3b1', '#eddea4', '#f7a072', '#ff9b42', '#fe5f55', '#f0b67f'];
        nodes.forEach(function (node, index) {
          if (!(node.group in groupColors)) {
            const colorIndex = index % colorPalette.length;
            const color = colorPalette[colorIndex];
            groupColors[
              node.group
            ] = color;
          }
        });
        
        const radiusScale = d3.scaleLinear()
        .domain([0, d3.max(nodes, (d) => d.numEle)]) // Define the input domain
        .range([7, 25]); // Define the output range for stroke-width

        const node = container
          .selectAll('.nodeGroup')
          .data(nodes)
          .join('g')
          .attr('class', 'node')
          .on('click', handleClick);

        node
          .append('circle')
          .attr("r", d => radiusScale(d.numEle))
          .attr('class', 'nodeCircle')
          .attr(
            'fill',
            (d) => d.groupColor
          );

        node
          .append('text')
          .attr('class', 'nodeName')
          .attr('dx', 10)
          .attr('dy', 5)
          .attr('font-size', 10)
          .attr('font-family', 'Share Tech')
          .text((d) => d.id)
          .style('visibility', (d) => {
            if (nodes.length < 20) {
              return 'visible';
            } else {
              return 'hidden';
            }
          });

        let selectedNode = null;

        function handleClick(event, d) {
          if (selectedNode === d) {
            // Second click on the same node, reset selection
            selectedNode = null;
            d3.selectAll('.nodeName').style(
              'visibility',
              (d) => {
                if (nodes.length < 10) {
                  return 'visible';
                } else {
                  return 'hidden';
                }
              }
            );
            resetData();
          } else {
            // First click or click on a different node
            d3.selectAll('.nodeName').style(
              'visibility',
              'hidden'
            );
            d3.select(this)
              .select('.nodeName')
              .style('visibility', 'visible');
            selectedNode = d;

            // Render children nodes
            graphRelation.forEach((rel) => {
              let currentChr = rel.children;
              if (rel.parentId === d.id) {
                currentLevel++;
                updateLevelData(
                  currentLevel,
                  currentChr
                );
                let currentData = getLevelData(
                  currentLevel
                );

                renderGraph(currentData);
              }

              // If there's no children nodes, click to highlight the node and its adjacent
              updateData(d);
              updateGraph();
            });
          }
        }

        function resetData() {
          // Reset selected node and update data
          // currentLevel = 0;
          selectedNode = null;
          updateData();
          updateGraph();
        }

        function updateData(selectedNode) {
          // Reset opacity for all nodes and links
          node.style('opacity', 1);
          link.style('opacity', 0.4);

          if (selectedNode) {
            const clickedNodeId = selectedNode.id;
            const neighborNodes = new Set();

            // Find neighbor nodes
            links.forEach((link) => {
              if (
                link.source.id === clickedNodeId
              ) {
                neighborNodes.add(link.target.id);
              } else if (
                link.target.id === clickedNodeId
              ) {
                neighborNodes.add(link.source.id);
              }
            });

            // Filter nodes and links
            node.style('opacity', (node) =>
              node.id === clickedNodeId ||
              neighborNodes.has(node.id)
                ? 1
                : 0.1
            );
            link.style('opacity', (link) =>
              link.source.id === clickedNodeId ||
              link.target.id === clickedNodeId
                ? 1
                : 0.1
            );
            // show the labels of neighbors
            node
              .select('.nodeName')
              .style('visibility', (node) =>
                node.id === clickedNodeId ||
                neighborNodes.has(node.id)
                  ? 'visible'
                  : 'hidden'
              );
          }
        }

        function updateGraph() {
          // Update the visual elements of the graph
          link
            .attr('x1', (d) => d.source.x)
            .attr('y1', (d) => d.source.y)
            .attr('x2', (d) => d.target.x)
            .attr('y2', (d) => d.target.y)
            .attr('stroke', 'black');

          node.attr(
            'transform',
            (d) => `translate(${d.x}, ${d.y})`
          );
        }

        // Add drag behavior to the nodes
        node.call(
          d3.drag()
            .on('start', dragStart)
            .on('drag', dragging)
            .on('end', dragEnd)
        );

        // Define the drag functions
        function dragStart(event, d) {
          if (!event.active)
            simulation.alphaTarget(0.3).restart();
          d.fx = d.x;
          d.fy = d.y;
        }

        function dragging(event, d) {
          d.fx = event.x;
          d.fy = event.y;
        }

        function dragEnd(event, d) {
          if (!event.active)
            simulation.alphaTarget(0);
          d.fx = null;
          d.fy = null;
        }

        // Update the node and link positions on each tick of the simulation
        simulation.on('tick', updateGraph);
      }
    })
    .catch((error) => {
      console.log('Error:', error);
    });

}(d3));

//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiaW5kZXguanMiLCJzb3VyY2VzIjpbImdyYXBoVXRpbHMuanMiLCJpbmRleC5qcyJdLCJzb3VyY2VzQ29udGVudCI6WyIvLyBpbXBvcnQgeyByZWFkRmlsZVN5bmMsIHdyaXRlRmlsZSwgcmVhZEZpbGUgfSBmcm9tICdmcydcblxuZXhwb3J0IGZ1bmN0aW9uIHJlYWRUeHQodXJsKXtcbiAgZmV0Y2godXJsKVxuICAgIC50aGVuKChyZXNwb25zZSkgPT4gcmVzcG9uc2UudGV4dCgpKVxuICAgIC50aGVuKChmaWxlQ29udGVudCkgPT4ge1xuICAgICAgcmV0dXJuIGZpbGVDb250ZW50XG4gICAgfSlcbiAgICAuY2F0Y2goKGVycm9yKSA9PiB7XG4gICAgICBjb25zb2xlLmxvZyhcbiAgICAgICAgJ0Vycm9yIG9jY3VycmVkIHdoaWxlIHJlYWRpbmcgdGhlIGZpbGU6JyxcbiAgICAgICAgZXJyb3JcbiAgICAgICk7XG4gICAgfSk7XG59O1xuXG5leHBvcnQgZnVuY3Rpb24gcGFyc2VHcmFwaERhdGFzZXQoZGF0YXNldCkge1xuICAvLyBEZWZpbmUgYW4gZW1wdHkgZ3JhcGggb2JqZWN0XG4gIGNvbnN0IGdyYXBoID0ge1xuICAgIG5vZGVzOiBbXSxcbiAgICBsaW5rczogW10sXG4gIH07XG4gIFxuXG4gIC8vIFNwbGl0IHRoZSBkYXRhc2V0IGludG8gaW5kaXZpZHVhbCBsaW5lc1xuICBjb25zdCBsaW5lcyA9IGRhdGFzZXQuc3BsaXQoJ1xcbicpO1xuXG4gIC8vIEl0ZXJhdGUgb3ZlciBlYWNoIGxpbmUgYW5kIHBhcnNlIHRoZSBub2RlIGFuZCB3ZWlnaHQgaW5mb3JtYXRpb25cbiAgbGluZXMuZm9yRWFjaCgobGluZSkgPT4ge1xuICAgIGNvbnN0IFtub2RlMSwgbm9kZTIsIHdlaWdodF0gPSBsaW5lXG4gICAgICAuc3BsaXQoJyAnKVxuICAgICAgLm1hcChOdW1iZXIpO1xuICAgIGlmIChcbiAgICAgIHR5cGVvZiBub2RlMSAhPT0gJ3VuZGVmaW5lZCcgJiZcbiAgICAgIHR5cGVvZiBub2RlMiAhPT0gJ3VuZGVmaW5lZCdcbiAgICApIHtcbiAgICAgIGNvbnN0IG5vZGUxSW5kZXggPSBncmFwaC5ub2Rlcy5maW5kSW5kZXgoXG4gICAgICAgIChub2RlKSA9PiBub2RlLmlkID09PSBub2RlMVxuICAgICAgKTtcbiAgICAgIGlmIChub2RlMUluZGV4ID09PSAtMSkge1xuICAgICAgICBncmFwaC5ub2Rlcy5wdXNoKHsgaWQ6IG5vZGUxLCBncm91cDogMSB9KTtcbiAgICAgIH1cbiAgICAgIGNvbnN0IG5vZGUySW5kZXggPSBncmFwaC5ub2Rlcy5maW5kSW5kZXgoXG4gICAgICAgIChub2RlKSA9PiBub2RlLmlkID09PSBub2RlMlxuICAgICAgKTtcbiAgICAgIGlmIChub2RlMkluZGV4ID09PSAtMSkge1xuICAgICAgICBncmFwaC5ub2Rlcy5wdXNoKHsgaWQ6IG5vZGUyLCBncm91cDogMiB9KTtcbiAgICAgIH1cblxuICAgICAgZ3JhcGgubGlua3MucHVzaCh7XG4gICAgICAgIHNvdXJjZTogbm9kZTEsXG4gICAgICAgIHRhcmdldDogbm9kZTIsXG4gICAgICAgIHZhbHVlOiB3ZWlnaHQsXG4gICAgICB9KTtcbiAgICB9XG4gIH0pO1xuXG4gIHJldHVybiBncmFwaDtcbn1cblxuZXhwb3J0IGZ1bmN0aW9uIHBhcnNlR3JhcGhEYXRhc2V0TGFiZWwoZGF0YXNldCkge1xuICAvLyBEZWZpbmUgYW4gZW1wdHkgZ3JhcGggb2JqZWN0XG4gIGNvbnN0IG1hcElkVG9MYWJlbCA9IHt9O1xuXG4gIC8vIFNwbGl0IHRoZSBkYXRhc2V0IGludG8gaW5kaXZpZHVhbCBsaW5lc1xuICBjb25zdCBsaW5lcyA9IGRhdGFzZXQuc3BsaXQoJ1xcbicpO1xuXG4gIC8vIEl0ZXJhdGUgb3ZlciBlYWNoIGxpbmUgYW5kIHBhcnNlIHRoZSBub2RlIGFuZCB3ZWlnaHQgaW5mb3JtYXRpb25cbiAgbGluZXMuZm9yRWFjaCgobGluZSkgPT4ge1xuICAgIGNvbnN0IFtpZCwgbGFiZWxdID0gbGluZVxuICAgICAgLnNwbGl0KCcgJylcbiAgICAgIC5tYXAoTnVtYmVyKTtcbiAgICBpZiAoXG4gICAgICB0eXBlb2YgaWQgIT09ICd1bmRlZmluZWQnICYmXG4gICAgICB0eXBlb2YgbGFiZWwgIT09ICd1bmRlZmluZWQnXG4gICAgKSB7XG4gICAgICBpZiAoIShpZCBpbiBtYXBJZFRvTGFiZWwpKSB7XG4gICAgICAgIG1hcElkVG9MYWJlbFtpZF0gPSBsYWJlbDtcbiAgICAgIH1cbiAgICB9XG4gIH0pO1xuXG4gIHJldHVybiBtYXBJZFRvTGFiZWw7XG59XG5cbmV4cG9ydCBmdW5jdGlvbiBwYXJzZUdyYXBoRGF0YXNldFNOQVAoXG4gIGRhdGFzZXQsXG4gIGlkMk5hbWU9bnVsbCxcbiAgbWFwSWRUb0xhYmVsPW51bGwsXG4pIHtcbiAgY29uc3QgZ3JhcGggPSB7XG4gICAgbm9kZXM6IFtdLFxuICAgIGxpbmtzOiBbXSxcbiAgfTtcbiAgY29uc3QgcHJvY2Vzc2VkTGlua3MgPSBuZXcgU2V0KCk7IC8vIENyZWF0ZSBhIFNldCB0byB0cmFjayBwcm9jZXNzZWQgbGluayBJRHNcblxuICBjb25zdCBsaW5lcyA9IGRhdGFzZXQuc3BsaXQoJ1xcbicpO1xuXG4gIGZvciAoY29uc3QgbGluZSBvZiBsaW5lcykge1xuICAgIGNvbnN0IFtub2RlMSwgbm9kZTJdID0gbGluZS5zcGxpdCgvWyBcXHRdKy8pO1xuICAgIGlmICgvXlxcZCtcXHMrXFxkKyQvLnRlc3QobGluZSkgfHwgL15cXGQrXFx0K1xcZCskLy50ZXN0KGxpbmUpKSB7XG5cbiAgICAgIGlmIChub2RlMSAhPT0gJycgJiYgbm9kZTIgIT09ICcnKSB7XG4gICAgICAgIGNvbnN0IG5vZGUxSUQgPSBpZDJOYW1lID8gaWQyTmFtZVtub2RlMS50b1N0cmluZygpXSA6IG5vZGUxLnRvU3RyaW5nKCk7XG4gICAgICAgIGNvbnN0IG5vZGUySUQgPSBpZDJOYW1lID8gaWQyTmFtZVtub2RlMi50b1N0cmluZygpXSA6IG5vZGUyLnRvU3RyaW5nKCk7XG5cbiAgICAgICAgbGV0IG5vZGUxSW5kZXggPSBncmFwaC5ub2Rlcy5maW5kSW5kZXgoKG5vZGUpID0+IG5vZGUuaWQgPT09IG5vZGUxSUQpO1xuICAgICAgICBpZiAobm9kZTFJbmRleCA9PT0gLTEpIHtcbiAgICAgICAgICBjb25zdCBub2RlMURhdGEgPSB7IGlkOiBub2RlMUlEIH07XG4gICAgICAgICAgaWYgKG1hcElkVG9MYWJlbCAmJiBtYXBJZFRvTGFiZWxbbm9kZTFJRF0pIHtcbiAgICAgICAgICAgIG5vZGUxRGF0YS5ncm91cCA9IG1hcElkVG9MYWJlbFtub2RlMUlEXTtcbiAgICAgICAgICB9XG4gICAgICAgICAgZ3JhcGgubm9kZXMucHVzaChub2RlMURhdGEpO1xuICAgICAgICB9XG5cbiAgICAgICAgbGV0IG5vZGUySW5kZXggPSBncmFwaC5ub2Rlcy5maW5kSW5kZXgoKG5vZGUpID0+IG5vZGUuaWQgPT09IG5vZGUySUQpO1xuICAgICAgICBpZiAobm9kZTJJbmRleCA9PT0gLTEpIHsgXG4gICAgICAgICAgY29uc3Qgbm9kZTJEYXRhID0geyBpZDogbm9kZTJJRCB9O1xuICAgICAgICAgIGlmIChtYXBJZFRvTGFiZWwgJiYgbWFwSWRUb0xhYmVsW25vZGUySURdKSB7XG4gICAgICAgICAgICBub2RlMkRhdGEuZ3JvdXAgPSBtYXBJZFRvTGFiZWxbbm9kZTJJRF07XG4gICAgICAgICAgfVxuICAgICAgICAgIGdyYXBoLm5vZGVzLnB1c2gobm9kZTJEYXRhKTtcbiAgICAgICAgfVxuXG4gICAgICAgIGNvbnN0IGxpbmtJZDEgPSBub2RlMUlEICsgJy0nICsgbm9kZTJJRDtcbiAgICAgICAgY29uc3QgbGlua0lkMiA9IG5vZGUySUQgKyAnLScgKyBub2RlMUlEO1xuXG4gICAgICAgIGlmICghcHJvY2Vzc2VkTGlua3MuaGFzKGxpbmtJZDEpICYmICFwcm9jZXNzZWRMaW5rcy5oYXMobGlua0lkMikpIHtcbiAgICAgICAgICBncmFwaC5saW5rcy5wdXNoKHsgc291cmNlOiBub2RlMUlELCB0YXJnZXQ6IG5vZGUySUQgfSk7XG4gICAgICAgICAgcHJvY2Vzc2VkTGlua3MuYWRkKGxpbmtJZDEpO1xuICAgICAgICAgIHByb2Nlc3NlZExpbmtzLmFkZChsaW5rSWQyKTtcbiAgICAgICAgfVxuICAgICAgfVxuICAgIH1cbiAgICB9XG5cbiAgcmV0dXJuIGdyYXBoO1xufVxuXG5leHBvcnQgZnVuY3Rpb24gcGFyc2VHcmFwaEZyb21KU09OKGpzb25UZXh0KSB7XG4gIGNvbnN0IGpzb25EYXRhID0gSlNPTi5wYXJzZShqc29uVGV4dCk7XG4gIFxuICBjb25zdCBncmFwaCA9IHtcbiAgICBub2RlczogW10sXG4gICAgbGlua3M6IFtdLFxuICB9O1xuXG4gIGZvciAoY29uc3Qgbm9kZSBvZiBqc29uRGF0YS5ub2Rlcykge1xuICAgIGdyYXBoLm5vZGVzLnB1c2goe2lkOiBub2RlLmlkLCBncm91cDogbm9kZS5ncm91cH0pO1xuICB9XG5cbiAgZm9yIChjb25zdCBsaW5rIG9mIGpzb25EYXRhLmxpbmtzKSB7XG4gICAgY29uc3Qgc291cmNlID0gbGluay5zb3VyY2U7XG4gICAgY29uc3QgdGFyZ2V0ID0gbGluay50YXJnZXQ7XG4gICAgZ3JhcGgubGlua3MucHVzaCh7IHNvdXJjZSwgdGFyZ2V0IH0pO1xuICB9XG5cbiAgcmV0dXJuIGdyYXBoO1xufVxuXG5leHBvcnQgZnVuY3Rpb24gc2VsZWN0UGFyc2VGb3JtYXQodGV4dCwgZmlsZUNvbnRlbnQsIGlkMk5hbWU9bnVsbCl7XG4gIGNvbnNvbGUubG9nKHRleHQpO1xuICBjb25zb2xlLmxvZyh0ZXh0LmVuZHNXaXRoKCcuanNvbicpKTtcbiAgbGV0IGdyYXBoO1xuICBpZiAodGV4dC5lbmRzV2l0aCgnLmpzb24nKSkge1xuICAgIC8vIGNvbnNvbGUubG9nKCdoZXJlJyk7XG4gICAgZ3JhcGggPSBwYXJzZUdyYXBoRnJvbUpTT04oZmlsZUNvbnRlbnQpO1xuICB9IFxuICBlbHNlIGlmICh0ZXh0LmVuZHNXaXRoKCcudHh0JykpIHtcbiAgICBcdGlmIChpZDJOYW1lKSB7XG4gICAgICBcdGdyYXBoID0gcGFyc2VHcmFwaERhdGFzZXRTTkFQKGZpbGVDb250ZW50LCBpZDJOYW1lKTtcbiAgICAgIH1cbiAgICAgIGVsc2Uge1xuICAgICAgICAvLyBjb25zb2xlLmxvZygnaGVyZScpO1xuICAgICAgICBncmFwaCA9IHBhcnNlR3JhcGhEYXRhc2V0U05BUChmaWxlQ29udGVudCk7XG4gICAgICB9XG4gIH0gZWxzZSB7XG4gICAgY29uc29sZS5lcnJvcihcIkVycm9yOiBVbnN1cHBvcnRlZCBmb3JtYXRcIik7XG4gIH1cbiAgY29uc29sZS5sb2coZ3JhcGgpO1xuICByZXR1cm4gZ3JhcGg7XG59XG5cbmV4cG9ydCBmdW5jdGlvbiBzd2FwS2V5c0FuZFZhbHVlcyhqc29uRGF0YSkge1xuICB0cnkge1xuICAgIGNvbnN0IHN3YXBwZWREYXRhID0ge307XG5cbiAgICBmb3IgKGNvbnN0IGtleSBpbiBqc29uRGF0YSkge1xuICAgICAgaWYgKGpzb25EYXRhLmhhc093blByb3BlcnR5KGtleSkpIHtcbiAgICAgICAgY29uc3QgdmFsdWUgPSBqc29uRGF0YVtrZXldO1xuICAgICAgICBzd2FwcGVkRGF0YVt2YWx1ZS50b1N0cmluZygpXSA9IGtleTtcbiAgICAgIH1cbiAgICB9XG5cbiAgICByZXR1cm4gc3dhcHBlZERhdGE7XG4gIH0gY2F0Y2ggKGVycm9yKSB7XG4gICAgY29uc29sZS5lcnJvcignSW52YWxpZCBKU09OIGRhdGE6JywgZXJyb3IpO1xuICAgIHJldHVybiBudWxsO1xuICB9XG59XG5cbi8vIFNpbXBsZSBSYW5kb20gV2Fsay1CYXNlZCBOb2RlIEVtYmVkZGluZ1xuZXhwb3J0IGZ1bmN0aW9uIG5vZGVFbWJlZGRpbmcoXG4gIGdyYXBoLFxuICBudW1XYWxrcyxcbiAgd2Fsa0xlbmd0aCxcbiAgZGltZW5zaW9ucyxcbiAgc2VlZFxuKSB7XG4gIGNvbnN0IHJuZyA9IHNlZWRlZFJhbmRvbShzZWVkKTsgLy8gSW5pdGlhbGl6ZSB0aGUgUk5HIHdpdGggdGhlIHNlZWQuXG5cbiAgLy8gSW5pdGlhbGl6ZSBlbWJlZGRpbmdzIGZvciBlYWNoIG5vZGVcbiAgY29uc3QgZW1iZWRkaW5ncyA9IHt9O1xuXG4gIC8vIFBlcmZvcm0gcmFuZG9tIHdhbGtzXG4gIGZvciAobGV0IHdhbGsgPSAwOyB3YWxrIDwgbnVtV2Fsa3M7IHdhbGsrKykge1xuICAgIGdyYXBoLmZvckVhY2hOb2RlKChub2RlKSA9PiB7XG4gICAgICBjb25zdCB3YWxrUGF0aCA9IFtub2RlXTtcblxuICAgICAgLy8gR2VuZXJhdGUgcmFuZG9tIHdhbGsgcGF0aCB1c2luZyBzZWVkZWQgUk5HXG4gICAgICBmb3IgKGxldCBzdGVwID0gMTsgc3RlcCA8IHdhbGtMZW5ndGg7IHN0ZXArKykge1xuICAgICAgICBjb25zdCBuZWlnaGJvcnMgPSBncmFwaC5vdXROZWlnaGJvcnMobm9kZSk7XG4gICAgICAgIGlmIChuZWlnaGJvcnMubGVuZ3RoID09PSAwKSBicmVhazsgLy8gSWYgbm8gbmVpZ2hib3JzLCBlbmQgdGhlIHdhbGsuXG4gICAgICAgIGNvbnN0IG5leHROb2RlID0gbmVpZ2hib3JzW01hdGguZmxvb3Iocm5nKCkgKiBuZWlnaGJvcnMubGVuZ3RoKV07IC8vIFVzZSB0aGUgc2VlZGVkIFJORy5cbiAgICAgICAgd2Fsa1BhdGgucHVzaChuZXh0Tm9kZSk7XG4gICAgICAgIG5vZGUgPSBuZXh0Tm9kZTtcbiAgICAgIH1cblxuICAgICAgLy8gVXBkYXRlIG5vZGUgZW1iZWRkaW5nc1xuICAgICAgZm9yIChjb25zdCBuIG9mIHdhbGtQYXRoKSB7XG4gICAgICAgIGlmICghZW1iZWRkaW5nc1tuXSkge1xuICAgICAgICAgIGVtYmVkZGluZ3Nbbl0gPSBuZXcgQXJyYXkoZGltZW5zaW9ucykuZmlsbCgwKTtcbiAgICAgICAgfVxuICAgICAgICBlbWJlZGRpbmdzW25dW3dhbGsgJSBkaW1lbnNpb25zXSArPSAxO1xuICAgICAgfVxuICAgIH0pO1xuICB9XG5cbiAgLy8gTm9ybWFsaXplIGVtYmVkZGluZ3NcbiAgZm9yIChjb25zdCBub2RlIGluIGVtYmVkZGluZ3MpIHtcbiAgICBjb25zdCBlbWJlZGRpbmcgPSBlbWJlZGRpbmdzW25vZGVdO1xuICAgIGNvbnN0IHN1bSA9IGVtYmVkZGluZy5yZWR1Y2UoKGFjYywgdmFsKSA9PiBhY2MgKyB2YWwsIDApO1xuICAgIGVtYmVkZGluZ3Nbbm9kZV0gPSBlbWJlZGRpbmcubWFwKCh2YWwpID0+IHZhbCAvIHN1bSk7XG4gIH1cblxuICByZXR1cm4gZW1iZWRkaW5ncztcbn1cblxuZXhwb3J0IGZ1bmN0aW9uIGtNZWFuc0NsdXN0ZXJpbmcoXG4gIGVtYmVkZGluZ3MsXG4gIG51bUNsdXN0ZXJzLFxuICBtYXhJdGVyYXRpb25zLFxuICBzZWVkXG4pIHtcbiAgY29uc3Qgbm9kZXMgPSBPYmplY3Qua2V5cyhlbWJlZGRpbmdzKTtcbiAgY29uc3QgZW1iZWRkaW5nc01hdHJpeCA9IG5vZGVzLm1hcCgobm9kZSkgPT4gZW1iZWRkaW5nc1tub2RlXSk7XG5cbiAgLy8gSW5pdGlhbGl6ZSBjZW50cm9pZHMgdXNpbmcgc2VlZGVkIHJhbmRvbSBudW1iZXIgZ2VuZXJhdG9yXG4gIGNvbnN0IHJuZyA9IHNlZWRlZFJhbmRvbShzZWVkKTtcbiAgbGV0IGNlbnRyb2lkcyA9IGVtYmVkZGluZ3NNYXRyaXguc2xpY2UoMCwgbnVtQ2x1c3RlcnMpLm1hcCgoKSA9PiB7XG4gICAgY29uc3QgcmFuZG9tSW5kZXggPSBNYXRoLmZsb29yKHJuZygpICogZW1iZWRkaW5nc01hdHJpeC5sZW5ndGgpO1xuICAgIHJldHVybiBlbWJlZGRpbmdzTWF0cml4W3JhbmRvbUluZGV4XTtcbiAgfSk7XG5cbiAgbGV0IGNsdXN0ZXJzID0gW107XG5cbiAgLy8gSXRlcmF0ZSB1bnRpbCBjb252ZXJnZW5jZSBvciBtYXhpbXVtIGl0ZXJhdGlvbnMgcmVhY2hlZFxuICBmb3IgKGxldCBpdGVyYXRpb24gPSAwOyBpdGVyYXRpb24gPCBtYXhJdGVyYXRpb25zOyBpdGVyYXRpb24rKykge1xuICAgIC8vIEFzc2lnbiBub2RlcyB0byBjbHVzdGVyc1xuICAgIGNsdXN0ZXJzID0gbmV3IEFycmF5KG51bUNsdXN0ZXJzKS5maWxsKCkubWFwKCgpID0+IFtdKTtcblxuICAgIGZvciAoY29uc3QgW25vZGVJbmRleCwgZW1iZWRkaW5nXSBvZiBlbWJlZGRpbmdzTWF0cml4LmVudHJpZXMoKSkge1xuICAgICAgbGV0IG1pbkRpc3RhbmNlID0gSW5maW5pdHk7XG4gICAgICBsZXQgYXNzaWduZWRDbHVzdGVyID0gLTE7XG5cbiAgICAgIGZvciAoY29uc3QgW2NlbnRyb2lkSW5kZXgsIGNlbnRyb2lkXSBvZiBjZW50cm9pZHMuZW50cmllcygpKSB7XG4gICAgICAgIGNvbnN0IGRpc3RhbmNlID0gZXVjbGlkZWFuRGlzdGFuY2UoZW1iZWRkaW5nLCBjZW50cm9pZCk7XG5cbiAgICAgICAgaWYgKGRpc3RhbmNlIDwgbWluRGlzdGFuY2UpIHtcbiAgICAgICAgICBtaW5EaXN0YW5jZSA9IGRpc3RhbmNlO1xuICAgICAgICAgIGFzc2lnbmVkQ2x1c3RlciA9IGNlbnRyb2lkSW5kZXg7XG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgY2x1c3RlcnNbYXNzaWduZWRDbHVzdGVyXS5wdXNoKG5vZGVzW25vZGVJbmRleF0pO1xuICAgIH1cblxuICAgIC8vIFVwZGF0ZSBjZW50cm9pZHNcbiAgICBjb25zdCBuZXdDZW50cm9pZHMgPSBbXTtcblxuICAgIGZvciAoY29uc3QgY2x1c3RlciBvZiBjbHVzdGVycykge1xuICAgICAgaWYgKGNsdXN0ZXIubGVuZ3RoID09PSAwKSB7XG4gICAgICAgIC8vIFVzZSBzZWVkZWQgcmFuZG9tIG51bWJlciBnZW5lcmF0b3IgdG8gcGljayBhIHJhbmRvbSBjZW50cm9pZFxuICAgICAgICBjb25zdCByYW5kb21JbmRleCA9IE1hdGguZmxvb3Iocm5nKCkgKiBlbWJlZGRpbmdzTWF0cml4Lmxlbmd0aCk7XG4gICAgICAgIG5ld0NlbnRyb2lkcy5wdXNoKGVtYmVkZGluZ3NNYXRyaXhbcmFuZG9tSW5kZXhdKTtcbiAgICAgIH0gZWxzZSB7XG4gICAgICAgIGNvbnN0IGNsdXN0ZXJFbWJlZGRpbmdzID0gY2x1c3Rlci5tYXAoKG5vZGUpID0+IGVtYmVkZGluZ3Nbbm9kZV0pO1xuICAgICAgICBjb25zdCBjbHVzdGVyQ2VudHJvaWQgPSBjYWxjdWxhdGVDZW50cm9pZChjbHVzdGVyRW1iZWRkaW5ncyk7XG4gICAgICAgIG5ld0NlbnRyb2lkcy5wdXNoKGNsdXN0ZXJDZW50cm9pZCk7XG4gICAgICB9XG4gICAgfVxuXG4gICAgLy8gQ2hlY2sgY29udmVyZ2VuY2VcbiAgICBpZiAoYXJlQ2VudHJvaWRzRXF1YWwoY2VudHJvaWRzLCBuZXdDZW50cm9pZHMpKSB7XG4gICAgICBicmVhaztcbiAgICB9XG5cbiAgICBjZW50cm9pZHMgPSBuZXdDZW50cm9pZHM7XG4gIH1cblxuICAvLyBDcmVhdGUgYSBkaWN0aW9uYXJ5IHdpdGggbm9kZXMgYXMga2V5cyBhbmQgY2x1c3RlciBjb3VudHMgYXMgdmFsdWVzXG4gIGNvbnN0IGNvbW11bml0eURpY3QgPSB7fTtcbiAgY29uc3Qgbm9kZUNsdXN0ZXJEaWN0ID0ge307XG4gIGNsdXN0ZXJzLmZvckVhY2goKGNsdXN0ZXIsIGluZGV4KSA9PiB7XG4gICAgbGV0IHBhcmVudElkID0gJ0NvbS4nK2luZGV4O1xuICAgIGNvbW11bml0eURpY3RbcGFyZW50SWRdID0gY2x1c3RlcjtcbiAgICBjbHVzdGVyLmZvckVhY2goKG5vZGUpID0+IHtcbiAgICAgIG5vZGVDbHVzdGVyRGljdFtub2RlXSA9IGluZGV4O1xuICAgIH0pO1xuICB9KTtcbiAgXG4gIHJldHVybiBbY2x1c3RlcnMsIGNvbW11bml0eURpY3QsIG5vZGVDbHVzdGVyRGljdF07XG59XG5cbmZ1bmN0aW9uIGZpbmRGYXJ0aGVzdE5vZGUoZW1iZWRkaW5ncywgY2VudHJvaWRzKSB7XG4gIGxldCBmYXJ0aGVzdE5vZGUgPSBudWxsO1xuICBsZXQgbWF4RGlzdGFuY2UgPSAtMTtcblxuICBmb3IgKGNvbnN0IG5vZGUgaW4gZW1iZWRkaW5ncykge1xuICAgIGNvbnN0IG5vZGVFbWJlZGRpbmcgPSBlbWJlZGRpbmdzW25vZGVdO1xuICAgIGxldCBtaW5EaXN0YW5jZSA9IEluZmluaXR5O1xuXG4gICAgZm9yIChjb25zdCBjZW50cm9pZCBvZiBjZW50cm9pZHMpIHtcbiAgICAgIGNvbnN0IGRpc3RhbmNlID0gZXVjbGlkZWFuRGlzdGFuY2Uobm9kZUVtYmVkZGluZywgY2VudHJvaWQpO1xuXG4gICAgICBpZiAoZGlzdGFuY2UgPCBtaW5EaXN0YW5jZSkge1xuICAgICAgICBtaW5EaXN0YW5jZSA9IGRpc3RhbmNlO1xuICAgICAgfVxuICAgIH1cblxuICAgIGlmIChtaW5EaXN0YW5jZSA+IG1heERpc3RhbmNlKSB7XG4gICAgICBtYXhEaXN0YW5jZSA9IG1pbkRpc3RhbmNlO1xuICAgICAgZmFydGhlc3ROb2RlID0gbm9kZTtcbiAgICB9XG4gIH1cblxuICByZXR1cm4gZmFydGhlc3ROb2RlO1xufVxuXG4vLyBIZWxwZXIgZnVuY3Rpb25zXG5mdW5jdGlvbiBldWNsaWRlYW5EaXN0YW5jZSh2ZWN0b3IxLCB2ZWN0b3IyKSB7XG4gIGxldCBzdW0gPSAwO1xuICBmb3IgKGxldCBpID0gMDsgaSA8IHZlY3RvcjEubGVuZ3RoOyBpKyspIHtcbiAgICBzdW0gKz0gTWF0aC5wb3codmVjdG9yMVtpXSAtIHZlY3RvcjJbaV0sIDIpO1xuICB9XG4gIHJldHVybiBNYXRoLnNxcnQoc3VtKTtcbn1cblxuZnVuY3Rpb24gY2FsY3VsYXRlQ2VudHJvaWQoZW1iZWRkaW5ncykge1xuICBjb25zdCBudW1EaW1lbnNpb25zID0gZW1iZWRkaW5nc1swXS5sZW5ndGg7XG4gIGNvbnN0IGNlbnRyb2lkID0gbmV3IEFycmF5KG51bURpbWVuc2lvbnMpLmZpbGwoXG4gICAgMFxuICApO1xuXG4gIGZvciAoY29uc3QgZW1iZWRkaW5nIG9mIGVtYmVkZGluZ3MpIHtcbiAgICBmb3IgKGxldCBpID0gMDsgaSA8IG51bURpbWVuc2lvbnM7IGkrKykge1xuICAgICAgY2VudHJvaWRbaV0gKz0gZW1iZWRkaW5nW2ldO1xuICAgIH1cbiAgfVxuXG4gIGZvciAobGV0IGkgPSAwOyBpIDwgbnVtRGltZW5zaW9uczsgaSsrKSB7XG4gICAgY2VudHJvaWRbaV0gLz0gZW1iZWRkaW5ncy5sZW5ndGg7XG4gIH1cblxuICByZXR1cm4gY2VudHJvaWQ7XG59XG5cbmZ1bmN0aW9uIGFyZUNlbnRyb2lkc0VxdWFsKFxuICBjZW50cm9pZHMxLFxuICBjZW50cm9pZHMyXG4pIHtcbiAgZm9yIChsZXQgaSA9IDA7IGkgPCBjZW50cm9pZHMxLmxlbmd0aDsgaSsrKSB7XG4gICAgY29uc3QgY2VudHJvaWQxID0gY2VudHJvaWRzMVtpXTtcbiAgICBjb25zdCBjZW50cm9pZDIgPSBjZW50cm9pZHMyW2ldO1xuXG4gICAgZm9yIChsZXQgaiA9IDA7IGogPCBjZW50cm9pZDEubGVuZ3RoOyBqKyspIHtcbiAgICAgIGlmIChjZW50cm9pZDFbal0gIT09IGNlbnRyb2lkMltqXSkge1xuICAgICAgICByZXR1cm4gZmFsc2U7XG4gICAgICB9XG4gICAgfVxuICB9XG5cbiAgcmV0dXJuIHRydWU7XG59XG5cblxuLy8gaGVscGVyIGZ1bmN0aW9uXG5mdW5jdGlvbiBjbnRCVFdDbHVzdChjbHVzdGVycywgZWRnZUxpc3QpIHtcbiAgY29uc3QgY2x1c3RlckVkZ2VzQ291bnQgPSB7fTtcblxuICBmb3IgKGxldCBpID0gMDsgaSA8IGNsdXN0ZXJzLmxlbmd0aDsgaSsrKSB7XG4gICAgZm9yIChsZXQgaiA9IGkgKyAxOyBqIDwgY2x1c3RlcnMubGVuZ3RoOyBqKyspIHtcbiAgICAgIGNvbnN0IGNsdXN0ZXJQYWlyID0gYCR7aX0tJHtqfWA7XG4gICAgICBjbHVzdGVyRWRnZXNDb3VudFtjbHVzdGVyUGFpcl0gPSAwO1xuICAgIH1cbiAgfVxuXG4gIGZvciAoY29uc3QgZWRnZSBvZiBlZGdlTGlzdCkge1xuICAgIGNvbnN0IHNvdXJjZUNsdXN0ZXJJbmRleCA9IGNsdXN0ZXJzLmZpbmRJbmRleChjbHVzdGVyID0+IGNsdXN0ZXIuaW5jbHVkZXMoZWRnZS5zb3VyY2UudG9TdHJpbmcoKSkpO1xuICAgIGNvbnN0IHRhcmdldENsdXN0ZXJJbmRleCA9IGNsdXN0ZXJzLmZpbmRJbmRleChjbHVzdGVyID0+IGNsdXN0ZXIuaW5jbHVkZXMoZWRnZS50YXJnZXQudG9TdHJpbmcoKSkpO1xuXG4gICAgaWYgKHNvdXJjZUNsdXN0ZXJJbmRleCAhPT0gdGFyZ2V0Q2x1c3RlckluZGV4KSB7XG4gICAgICBjb25zdCBjbHVzdGVyUGFpciA9IGAke01hdGgubWluKHNvdXJjZUNsdXN0ZXJJbmRleCwgdGFyZ2V0Q2x1c3RlckluZGV4KX0tJHtNYXRoLm1heChzb3VyY2VDbHVzdGVySW5kZXgsIHRhcmdldENsdXN0ZXJJbmRleCl9YDtcbiAgICAgIGNsdXN0ZXJFZGdlc0NvdW50W2NsdXN0ZXJQYWlyXSsrO1xuICAgIH1cbiAgfVxuICBcbiAgcmV0dXJuIGNsdXN0ZXJFZGdlc0NvdW50XG59O1xuXG5jbGFzcyBHcmFwaCB7XG4gIGNvbnN0cnVjdG9yKCkge1xuICAgIHRoaXMubm9kZXMgPSBuZXcgTWFwKCk7IC8vIE1hcCB0byBzdG9yZSBub2RlcyBhbmQgdGhlaXIgZWRnZXNcbiAgfVxuXG4gIGFkZE5vZGUobm9kZUlkKSB7XG4gICAgaWYgKCF0aGlzLm5vZGVzLmhhcyhub2RlSWQpKSB7XG4gICAgICB0aGlzLm5vZGVzLnNldChub2RlSWQsIG5ldyBTZXQoKSk7XG4gICAgfVxuICB9XG5cbiAgYWRkRWRnZShzb3VyY2VJZCwgdGFyZ2V0SWQpIHtcbiAgICAvLyBBZGQgc291cmNlIG5vZGUgaWYgbm90IGFscmVhZHkgcHJlc2VudFxuICAgIGlmICghdGhpcy5ub2Rlcy5oYXMoc291cmNlSWQpKSB7XG4gICAgICB0aGlzLmFkZE5vZGUoc291cmNlSWQpO1xuICAgIH1cblxuICAgIC8vIEFkZCB0YXJnZXQgbm9kZSBpZiBub3QgYWxyZWFkeSBwcmVzZW50XG4gICAgaWYgKCF0aGlzLm5vZGVzLmhhcyh0YXJnZXRJZCkpIHtcbiAgICAgIHRoaXMuYWRkTm9kZSh0YXJnZXRJZCk7XG4gICAgfVxuXG4gICAgLy8gQWRkIGVkZ2UgZnJvbSBzb3VyY2UgdG8gdGFyZ2V0XG4gICAgdGhpcy5ub2Rlcy5nZXQoc291cmNlSWQpLmFkZCh0YXJnZXRJZCk7XG4gIH1cblxuICBvdXROZWlnaGJvcnMobm9kZUlkKSB7XG4gICAgaWYgKHRoaXMubm9kZXMuaGFzKG5vZGVJZCkpIHtcbiAgICAgIHJldHVybiBbLi4udGhpcy5ub2Rlcy5nZXQobm9kZUlkKV07XG4gICAgfSBlbHNlIHtcbiAgICAgIHJldHVybiBbXTtcbiAgICB9XG4gIH1cblxuICBmb3JFYWNoTm9kZShjYWxsYmFjaykge1xuICAgIGZvciAobGV0IG5vZGVJZCBvZiB0aGlzLm5vZGVzLmtleXMoKSkge1xuICAgICAgY2FsbGJhY2sobm9kZUlkKTtcbiAgICB9XG4gIH1cbn1cblxuZXhwb3J0IGZ1bmN0aW9uIGNvbnN0cnVjdEdyYXBoKG5vZGVFZGdlcykge1xuICAvLyAvLyBFeGFtcGxlIHVzYWdlXG4gIGNvbnN0IGdyYXBoID0gbmV3IEdyYXBoKCk7XG4gIC8vIC4uLiBCdWlsZCBhbmQgcG9wdWxhdGUgdGhlIGdyYXBoIHdpdGggbm9kZXMgYW5kIGVkZ2VzXG4gIC8vIEFkZCBub2RlcyB0byB0aGUgZ3JhcGhcbiAgbm9kZUVkZ2VzLm5vZGVzLmZvckVhY2goKG5vZGUpID0+IHtcbiAgICBncmFwaC5hZGROb2RlKG5vZGUuaWQpO1xuICB9KTtcblxuICAvLyBBZGQgZWRnZXMgdG8gdGhlIGdyYXBoXG4gIG5vZGVFZGdlcy5saW5rcy5mb3JFYWNoKChsaW5rKSA9PiB7XG4gICAgZ3JhcGguYWRkRWRnZShsaW5rLnNvdXJjZSwgbGluay50YXJnZXQpO1xuICB9KTtcbiAgcmV0dXJuIGdyYXBoO1xufVxuXG5leHBvcnQgZnVuY3Rpb24gY29uc3RHcmFwaEZyb21DbHVzdChjbHVzdGVycywgZWRnZUxpc3QsIGNvbG9yUGFsZXR0ZSl7XG4gIGNvbnN0IGdyb3VwQ29sb3JzID0ge307XG4gIGxldCBjbHVzdGVyRWRnZXNDb3VudCA9IGNudEJUV0NsdXN0KGNsdXN0ZXJzLCBlZGdlTGlzdCk7XG4gIGxldCBjb21tTm9kZUVkZ2VzID0ge25vZGVzOiBbXSwgbGlua3M6IFtdfTtcblx0Zm9yIChsZXQgaSA9IDA7IGkgPCBjbHVzdGVycy5sZW5ndGg7IGkrKykge1xuICAgIGNvbnN0IGNvbG9ySW5kZXggPSBpICUgY29sb3JQYWxldHRlLmxlbmd0aDtcbiAgICBjb25zdCBjb2xvciA9IGNvbG9yUGFsZXR0ZVtjb2xvckluZGV4XTtcbiAgICBncm91cENvbG9yc1tpXSA9IGNvbG9yO1xuICAgIFxuICAgIC8vIENyZWF0ZSBub2RlcyBvZiB0aGUgZmlyc3QgY29tbXVuaXR5IGxheWVyXG4gIFx0bGV0IGN1cnJlbnRJZCA9ICdDb20uJyArIGk7XG4gICAgY29tbU5vZGVFZGdlcy5ub2Rlcy5wdXNoKHtpZDogY3VycmVudElkLCBncm91cDogaSwgZ3JvdXBDb2xvcjogZ3JvdXBDb2xvcnNbaV0sIG51bUVsZTpjbHVzdGVyc1tpXS5sZW5ndGh9KX07XG4gIFx0XG4gIFxuICAvLyBDcmVhdGUgZWRnZXMgb2YgdGhlIGZpcnN0IGNvbW11bml0eSBsYXllclxuICBPYmplY3QuZW50cmllcyhjbHVzdGVyRWRnZXNDb3VudCkuZm9yRWFjaCgoW2tleSwgdmFsdWVdKSA9PiB7XG4gICAgaWYodmFsdWUgPiAwKXtcbiAgICAgIGxldCBjb21tUGFpciA9IGtleS5zcGxpdCgnLScpO1xuICAgICAgbGV0IGN1cnJlbnRTcmNlID0gJ0NvbS4nICsgY29tbVBhaXJbMF07XG4gICAgICBsZXQgY3VycmVudFRHVCAgPSAnQ29tLicgKyBjb21tUGFpclsxXTtcbiAgICAgIGNvbW1Ob2RlRWRnZXMubGlua3MucHVzaCh7c291cmNlOiBjdXJyZW50U3JjZSwgdGFyZ2V0OiBjdXJyZW50VEdULCB2YWx1ZX0pO1xuICBcdH1cbiAgfSk7XG4gIHJldHVybiBbY29tbU5vZGVFZGdlcywgZ3JvdXBDb2xvcnNdO1xufTtcblxuXG5leHBvcnQgZnVuY3Rpb24gY29uc3RTdWJHRnJvbUNsdXN0KGNsdXN0LCBlZGdlTGlzdCwgbnVtQ2x1c3RlciwgZ3JvdXBDb2xvcnMsIG5vZGVDbHVzdGVyRGljdCwgY29uZmlnKXtcbiAgY29uc3QgbWF4TnVtTm9kZVBlckNsdXN0ZXIgPSBjb25maWdbMF07XG4gIGNvbnN0IG1heE51bUVkZ2VQZXJDbHVzdGVyID0gY29uZmlnWzFdO1xuICBjb25zdCBtYXhOdW1BZGpOb2RlUGVyQ2x1c3RlciA9IGNvbmZpZ1syXTtcbiAgY29uc3QgbWF4TnVtQWRqRWRnZVBlckNsdXN0ZXIgPSBjb25maWdbM107XG4gIFxuICAvLyBHZW5lcmF0ZSBlZGdlc1xuICBjb25zdCBjbHVzdGVyU2V0ID0gbmV3IFNldChjbHVzdCk7XG4gIGxldCBkYXRhID0ge307XG4gIFxuICAvLyBHZW5lcmF0ZSBub2RlcyBhbmQgbGlua3NcbiAgZGF0YS5ub2RlcyA9IEFycmF5LmZyb20oY2x1c3RlclNldCwgKGNsdXMpID0+ICh7IGlkOiBjbHVzLCBncm91cDogbnVtQ2x1c3RlciwgZ3JvdXBDb2xvcjogZ3JvdXBDb2xvcnNbbnVtQ2x1c3Rlcl0sIG51bUVsZTowfSkpO1xuICBkYXRhLmxpbmtzID0gZWRnZUxpc3QuZmlsdGVyKChlZGdlKSA9PiBjbHVzdGVyU2V0LmhhcyhlZGdlLnNvdXJjZSkgJiYgY2x1c3RlclNldC5oYXMoZWRnZS50YXJnZXQpKTtcbiAgaWYgKGRhdGEubm9kZXMubGVuZ3RoID4gbWF4TnVtTm9kZVBlckNsdXN0ZXIpIHtcbiAgXHRkYXRhID0gZmlsdGVyRGF0YUJ5VG9wTm9kZXMoZGF0YSwgbWF4TnVtTm9kZVBlckNsdXN0ZXIpO1xuICB9XG4gIFxuICBcbiAgZGF0YS5saW5rcyA9IGRhdGEubGlua3Muc2xpY2UoMCwgbWF4TnVtRWRnZVBlckNsdXN0ZXIpO1xuXG4gIGxldCBudW1BZGpFZGdlcyA9IDA7XG4gIGNvbnN0IGV4aXN0aW5nTm9kZXMgPSBuZXcgU2V0KGRhdGEubm9kZXMubWFwKChub2RlKSA9PiBub2RlLmlkKSk7IC8vIENyZWF0ZSBhIFNldCBvZiBleGlzdGluZyBub2RlIElEc1xuXG4gIGVkZ2VMaXN0LmZvckVhY2goKGVkZ2UpID0+IHtcbiAgICAvLyBFZGdlcyBvZiB0aGUgY3VycmVudCBjbHVzdGVyXG4gICAgY29uc3QgaXNFZGdlSW5EYXRhTGlua3MgPSBkYXRhLmxpbmtzLnNvbWUoKGV4aXN0aW5nRWRnZSkgPT4gKFxuICAgICAgKGV4aXN0aW5nRWRnZS5zb3VyY2UgPT09IGVkZ2Uuc291cmNlICYmIGV4aXN0aW5nRWRnZS50YXJnZXQgPT09IGVkZ2UudGFyZ2V0KSB8fFxuICAgICAgKGV4aXN0aW5nRWRnZS5zb3VyY2UgPT09IGVkZ2UudGFyZ2V0ICYmIGV4aXN0aW5nRWRnZS50YXJnZXQgPT09IGVkZ2Uuc291cmNlKVxuICAgICkpO1xuXG4gICAgY29uc3QgaXNTb3VyY2VJbkRhdGFOb2RlcyA9IGV4aXN0aW5nTm9kZXMuaGFzKGVkZ2Uuc291cmNlKTsgLy8gQ2hlY2sgaWYgc291cmNlIG5vZGUgYWxyZWFkeSBleGlzdHNcbiAgICBjb25zdCBpc1RhcmdldEluRGF0YU5vZGVzID0gZXhpc3RpbmdOb2Rlcy5oYXMoZWRnZS50YXJnZXQpOyAvLyBDaGVjayBpZiB0YXJnZXQgbm9kZSBhbHJlYWR5IGV4aXN0c1xuXG4gICAgLy8gRWRnZXMgb2YgYWRqYWNlbnQgY2x1c3RlcnNcbiAgICBpZiAoIWlzRWRnZUluRGF0YUxpbmtzICYmIG51bUFkakVkZ2VzIDw9IG1heE51bUFkakVkZ2VQZXJDbHVzdGVyICYmIChjbHVzdGVyU2V0LmhhcyhlZGdlLnNvdXJjZSkgfHwgY2x1c3RlclNldC5oYXMoZWRnZS50YXJnZXQpKSkge1xuICAgICAgZGF0YS5saW5rcy5wdXNoKGVkZ2UpO1xuICAgICAgbnVtQWRqRWRnZXMrKztcblxuICAgICAgaWYgKCFpc1NvdXJjZUluRGF0YU5vZGVzKSB7XG4gICAgICAgIGxldCBudW1DbHVzdGVyUyA9IG5vZGVDbHVzdGVyRGljdFtlZGdlLnNvdXJjZV07XG4gICAgICAgIGRhdGEubm9kZXMucHVzaCh7XG4gICAgICAgICAgaWQ6IGVkZ2Uuc291cmNlLFxuICAgICAgICAgIGdyb3VwOiBudW1DbHVzdGVyUyxcbiAgICAgICAgICBncm91cENvbG9yOiBncm91cENvbG9yc1tudW1DbHVzdGVyU10sXG4gICAgICAgICAgbnVtRWxlOjBcbiAgICAgICAgfSk7XG4gICAgICAgIGV4aXN0aW5nTm9kZXMuYWRkKGVkZ2Uuc291cmNlKTsgLy8gQWRkIHNvdXJjZSBub2RlIElEIHRvIHRoZSBleGlzdGluZ05vZGVzIHNldFxuICAgICAgfVxuICAgICAgaWYgKCFpc1RhcmdldEluRGF0YU5vZGVzKSB7XG4gICAgICAgIGxldCBudW1DbHVzdGVyVCA9IG5vZGVDbHVzdGVyRGljdFtlZGdlLnRhcmdldF07XG4gICAgICAgIGRhdGEubm9kZXMucHVzaCh7XG4gICAgICAgICAgaWQ6IGVkZ2UudGFyZ2V0LFxuICAgICAgICAgIGdyb3VwOiBudW1DbHVzdGVyVCxcbiAgICAgICAgICBncm91cENvbG9yOiBncm91cENvbG9yc1tudW1DbHVzdGVyVF0sXG4gICAgICAgICAgbnVtRWxlOjBcbiAgICAgICAgfSk7XG4gICAgICAgIGV4aXN0aW5nTm9kZXMuYWRkKGVkZ2UudGFyZ2V0KTsgLy8gQWRkIHRhcmdldCBub2RlIElEIHRvIHRoZSBleGlzdGluZ05vZGVzIHNldFxuICAgICAgfVxuICAgIH1cbiAgfSk7XG4gIHJldHVybiBkYXRhO1xufTtcblxuXG5cbmV4cG9ydCBmdW5jdGlvbiBmaWx0ZXJEYXRhQnlUb3BOb2RlcyhkYXRhLCBudW1Ub3BOb2Rlcykge1xuICBjb25zdCBlZGdlQ291bnQgPSB7fTtcbiAgZGF0YS5saW5rcy5mb3JFYWNoKGxpbmsgPT4ge1xuICAgIGNvbnN0IHNvdXJjZU5vZGVJRCA9IGxpbmsuc291cmNlO1xuICAgIGNvbnN0IHRhcmdldE5vZGVJRCA9IGxpbmsudGFyZ2V0O1xuXG4gICAgZWRnZUNvdW50W3NvdXJjZU5vZGVJRF0gPSAoZWRnZUNvdW50W3NvdXJjZU5vZGVJRF0gfHwgMCkgKyAxO1xuICAgIGVkZ2VDb3VudFt0YXJnZXROb2RlSURdID0gKGVkZ2VDb3VudFt0YXJnZXROb2RlSURdIHx8IDApICsgMTtcbiAgfSk7XG5cbiAgY29uc3QgdG9wTm9kZXMgPSBPYmplY3QuZW50cmllcyhlZGdlQ291bnQpXG4gICAgLnNvcnQoKGEsIGIpID0+IGJbMV0gLSBhWzFdKVxuICAgIC5zbGljZSgwLCBudW1Ub3BOb2RlcylcbiAgICAubWFwKChbbm9kZUlEXSkgPT4gbm9kZUlEKTtcblxuICBjb25zdCBmaWx0ZXJlZE5vZGVzID0gZGF0YS5ub2Rlcy5maWx0ZXIobm9kZSA9PiB0b3BOb2Rlcy5pbmNsdWRlcyhub2RlLmlkKSk7XG4gIGNvbnN0IGZpbHRlcmVkRWRnZXMgPSBkYXRhLmxpbmtzLmZpbHRlcihsaW5rID0+XG4gICAgdG9wTm9kZXMuaW5jbHVkZXMobGluay5zb3VyY2UpICYmIHRvcE5vZGVzLmluY2x1ZGVzKGxpbmsudGFyZ2V0KVxuICApO1xuICBcbiAgcmV0dXJuICB7bm9kZXM6IGZpbHRlcmVkTm9kZXMsIGxpbmtzOiBmaWx0ZXJlZEVkZ2VzLH1cbn1cblxuLy8gLS0tLS0tLS0tLS0tLSAgVXRpbGl0eSBmdW5jdGlvbnNcblxubGV0IGxldmVsc0RhdGEgPSB7fTtcbi8vIFVwZGF0ZSBkYXRhIGZvciBhIHNwZWNpZmljIGxldmVsXG5leHBvcnQgZnVuY3Rpb24gdXBkYXRlTGV2ZWxEYXRhKGxldmVsLCBkYXRhKSB7XG4gIGxldmVsc0RhdGFbbGV2ZWxdID0gZGF0YTtcbiAgY29uc29sZS5sb2cobGV2ZWxzRGF0YSk7XG59XG5cbi8vIEFjY2VzcyBkYXRhIGZvciBhIHNwZWNpZmljIGxldmVsXG5leHBvcnQgZnVuY3Rpb24gZ2V0TGV2ZWxEYXRhKGxldmVsKSB7XG4gIHJldHVybiBsZXZlbHNEYXRhW2xldmVsXTtcbn1cblxuLy8gQ3VzdG9tIGRldGVybWluaXN0aWMgcmFuZG9tIG51bWJlciBnZW5lcmF0b3JcbmZ1bmN0aW9uIHNlZWRlZFJhbmRvbShzZWVkKSB7XG4gIGxldCBtX3cgPSBzZWVkO1xuICBsZXQgbV96ID0gOTg3NjU0MzIxO1xuXG4gIC8vIFJldHVybnMgYSByYW5kb20gbnVtYmVyIGJldHdlZW4gMCBhbmQgMVxuICByZXR1cm4gZnVuY3Rpb24gKCkge1xuICAgIG1feiA9ICgzNjk2OSAqIChtX3ogJiA2NTUzNSkgKyAobV96ID4+IDE2KSkgJiAweGZmZmZmZmZmO1xuICAgIG1fdyA9ICgxODAwMCAqIChtX3cgJiA2NTUzNSkgKyAobV93ID4+IDE2KSkgJiAweGZmZmZmZmZmO1xuICAgIGxldCByZXN1bHQgPSAoKG1feiA8PCAxNikgKyBtX3cpICYgMHhmZmZmZmZmZjtcbiAgICByZXN1bHQgLz0gNDI5NDk2NzI5NjtcbiAgICByZXR1cm4gcmVzdWx0ICsgMC41O1xuICB9O1xufVxuIiwiaW1wb3J0IHtcbiAgc2VsZWN0LFxuICBzY2FsZUxpbmVhcixcbiAgZm9yY2VTaW11bGF0aW9uLFxuICBmb3JjZUxpbmssXG4gIGZvcmNlTWFueUJvZHksXG4gIGZvcmNlQ2VudGVyLFxuICBkcmFnLFxuICBtYXgsXG4gIHpvb20gYXMgZDNab29tLFxuICBzY2FsZU9yZGluYWwsXG4gIHNlbGVjdEFsbCxcbiAgZDNab29tSWRlbnRpdHksXG59IGZyb20gJ2QzJztcbmltcG9ydCB7XG4gIHBhcnNlR3JhcGhEYXRhc2V0LFxuICB1cGRhdGVMZXZlbERhdGEsXG4gIGdldExldmVsRGF0YSxcbiAgcmVhZFRleHQsXG4gIHBhcnNlR3JhcGhEYXRhc2V0U05BUCxcbiAgcGFyc2VHcmFwaERhdGFzZXRMYWJlbCxcbiAgcmVhZENsYXNzSWQsXG4gIG5vZGVFbWJlZGRpbmcsXG4gIGNvbnN0cnVjdEdyYXBoLFxuICBrTWVhbnNDbHVzdGVyaW5nLFxuICByZWFkVHh0LFxuICBjb25zdEdyYXBoRnJvbUNsdXN0LFxuICBjb25zdFN1YkdGcm9tQ2x1c3QsXG4gIGZpbHRlckRhdGFCeVRvcE5vZGVzLFxuICBzd2FwS2V5c0FuZFZhbHVlcyxcbiAgc2VsZWN0UGFyc2VGb3JtYXRcbn0gZnJvbSAnLi9ncmFwaFV0aWxzLmpzJztcbi8vIGltcG9ydCB7IHJlYWRGaWxlU3luYywgd3JpdGVGaWxlLCByZWFkRmlsZSB9IGZyb20gJ2ZzJ1xuXG5cbmNvbnN0IHdpZHRoID0gd2luZG93LmlubmVyV2lkdGg7XG5jb25zdCBoZWlnaHQgPSB3aW5kb3cuaW5uZXJIZWlnaHQ7XG5cbi8vIENyZWF0ZSBhbiBTVkcgY29udGFpbmVyIGZvciB0aGUgZ3JhcGhcbmNvbnN0IHN2ZyA9IHNlbGVjdCgnYm9keScpXG4gIC5hcHBlbmQoJ3N2ZycpXG4gIC5hdHRyKCd3aWR0aCcsIHdpZHRoKVxuICAuYXR0cignaGVpZ2h0JywgaGVpZ2h0KTtcblxuLy8gQ3JlYXRlIGEgPGc+IGVsZW1lbnQgYXMgdGhlIGNvbnRhaW5lciBmb3IgdGhlIGdyYXBoIGVsZW1lbnRzXG52YXIgY29udGFpbmVyID0gc3ZnLmFwcGVuZCgnZycpO1xuXG5mdW5jdGlvbiB6b29tZWQoZXZlbnQpIHtcbiAgY29udGFpbmVyLmF0dHIoJ3RyYW5zZm9ybScsIGV2ZW50LnRyYW5zZm9ybSk7XG59XG5cbi8vIENyZWF0ZSBhIHpvb20gYmVoYXZpb3JcbmNvbnN0IHpvb20gPSBkM1pvb20oKVxuICAuc2NhbGVFeHRlbnQoWzAuMSwgNF0pIC8vIFNldCB0aGUgem9vbSBzY2FsZSByYW5nZVxuICAub24oJ3pvb20nLCB6b29tZWQpOyAvLyBTcGVjaWZ5IHRoZSB6b29tIGV2ZW50IGhhbmRsZXJcblxuLy8gQXBwbHkgdGhlIHpvb20gYmVoYXZpb3IgdG8gdGhlIFNWRyBjb250YWluZXJcbnN2Zy5jYWxsKHpvb20pO1xuXG4vLy8vLy8vLy8vIFN0ZXAgMTogUmVhZCBkYXRhIGZyb20gYSBVUkwgXG5cbmFzeW5jIGZ1bmN0aW9uIGZldGNoRGF0YSgpIHtcbiAgdHJ5IHtcbiAgICBcbiAgICAvLyBGZXRjaCBjb25maWcuanNvblxuICAgIGNvbnN0IGNvbmZpZ1QgPSBhd2FpdCBmZXRjaCgnY29uZmlnLmpzb24nKTtcbiAgICBjb25zdCBjb25maWdUZXh0ID0gYXdhaXQgY29uZmlnVC50ZXh0KCk7XG4gICAgY29uc3QgY29uZmlnID0gSlNPTi5wYXJzZShjb25maWdUZXh0KTtcbiAgICBcbiAgICBjb25zdCBpZDJOYW1lVXJsID0gY29uZmlnLmlkMm5hbWVVUkw7XG4gICAgY29uc3Qgbm9kZXNFZGdlc1VybCA9IGNvbmZpZy5lZGdlc1VSTDtcbiAgICBjb25zdCBjb21tVXJsID0gY29uZmlnLmNvbW11bml0eUZpbGU7XG4gICAgLy8gaWQyTmFtZVVybCA9ICdodHRwczovL3Jhdy5naXRodWJ1c2VyY29udGVudC5jb20vS2hhbmluc2kvRDMtdmlzdWFsaXphdGlvbi9tYWluL2RhdGEvdXNlci5qc29uJztcbiAgICAvLyBsZXQgbm9kZXNFZGdlc1VybCA9ICdodHRwczovL3Jhdy5naXRodWJ1c2VyY29udGVudC5jb20vS2hhbmluc2kvRDMtdmlzdWFsaXphdGlvbi9tYWluL2RhdGEvZW1haWwtRXUtY29yZS50eHQnO1xuICAgIC8vIGxldCBub2Rlc0VkZ2VzVXJsID0gJ2h0dHBzOi8vcmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbS9LaGFuaW5zaS9EMy12aXN1YWxpemF0aW9uL21haW4vZGF0YS9lbnJvbl9lZGdlcy50eHQnO1xuXG4gICAgbGV0IG5vZGVMYWJlbE9iajtcbiAgICBsZXQgbm9kZUVkZ2VzO1xuICAgIGxldCBjb21wbGV4R3JhcGg7XG4gICAgbGV0IGlkMk5hbWVEaWN0O1xuICAgIGxldCBjb21tdW5pdHlEaWN0OztcbiAgICBcbiAgICAvLyBGZXRjaCBub2Rlc0VkZ2VzVXJsXG4gICAgY29uc3QgcmVzcG9uc2UxID0gYXdhaXQgZmV0Y2gobm9kZXNFZGdlc1VybCk7XG4gICAgY29uc3QgZmlsZUNvbnRlbnQxID0gYXdhaXQgcmVzcG9uc2UxLnRleHQoKTtcbiAgICBcbiAgICAvLyBGZXRjaCBpZDJOYW1lVXJsIGlmIGdpdmVuXG4gICAgaWYgKGlkMk5hbWVVcmwgJiYgaWQyTmFtZVVybC5sZW5ndGg+MCkge1xuICAgICAgY29uc3QgcmVzcG9uc2UyID0gYXdhaXQgZmV0Y2goaWQyTmFtZVVybCk7XG4gICAgICBjb25zdCByZXNwb25zZTJUZXh0ID0gYXdhaXQgcmVzcG9uc2UyLnRleHQoKTtcbiAgICAgIGNvbnN0IHJlc3BvbnNlMkpTID0gSlNPTi5wYXJzZShyZXNwb25zZTJUZXh0KTtcbiAgICAgIGlkMk5hbWVEaWN0ID0gIHN3YXBLZXlzQW5kVmFsdWVzKHJlc3BvbnNlMkpTKTtcbiAgICAgIG5vZGVFZGdlcyA9IHNlbGVjdFBhcnNlRm9ybWF0KG5vZGVzRWRnZXNVcmwsIGZpbGVDb250ZW50MSwgaWQyTmFtZURpY3QpO1xuICAgIH1cbiAgICBlbHNlIHtcblx0XHRcdG5vZGVFZGdlcyA9IHNlbGVjdFBhcnNlRm9ybWF0KG5vZGVzRWRnZXNVcmwsIGZpbGVDb250ZW50MSlcbiAgICB9XG5cbiAgICBpZiAoY29tbVVybC5sZW5ndGggPiAxKXtcbiAgICAgIGNvbnN0IHJlc3BvbnNlMyA9IGF3YWl0IGZldGNoKGNvbW1VcmwpO1xuICAgICAgY29uc3QgcmVzcG9uc2UzVGV4dCA9IGF3YWl0IHJlc3BvbnNlMy50ZXh0KCk7XG4gICAgICBjb21tdW5pdHlEaWN0ID0gSlNPTi5wYXJzZShyZXNwb25zZTNUZXh0KTtcbiAgICAgIC8vIGNvbnNvbGUubG9nKCdjb21tdW5pdHlEaWN0Jyk7XG4gICAgICAvLyBjb25zb2xlLmxvZyhjb21tdW5pdHlEaWN0KTtcbiAgICB9XG4gICAgXG4gICAgXG4gICAgXG4gICAgLy8gQ29udGludWUgd2l0aCB0aGUgbmV4dCBzZWN0aW9uIG9yIHBlcmZvcm0gYWRkaXRpb25hbCBvcGVyYXRpb25zXG4gICAgY29uc29sZS5sb2coJ0NvbXBsZXRlZCBmZXRjaGluZyBhbmQgcHJvY2Vzc2luZyBkYXRhLicpO1xuICAgIHJldHVybiB7XG4gICAgICBub2RlRWRnZXMsXG4gICAgICBjb25maWcsXG4gICAgICBjb21tdW5pdHlEaWN0LFxuICAgIH07XG4gIH0gY2F0Y2ggKGVycm9yKSB7XG4gICAgY29uc29sZS5sb2coJ0Vycm9yIG9jY3VycmVkIHdoaWxlIHJlYWRpbmcgdGhlIGZpbGU6JywgZXJyb3IpO1xuICB9XG59XG5cbi8vIENhbGwgdGhlIGZldGNoRGF0YSBmdW5jdGlvbiB0byBpbml0aWF0ZSB0aGUgZXhlY3V0aW9uXG5mZXRjaERhdGEoKVxuXHQudGhlbigoe25vZGVFZGdlcywgY29uZmlnLCBjb21tdW5pdHlEaWN0fSkgPT4ge1xuICAgIC8vIEV4YW1wbGU6IFBlcmZvcm0gb3BlcmF0aW9ucyB1c2luZyBub2RlRWRnZXNcbiAgICBjb25zb2xlLmxvZygnRmlyc3Qgbm9kZSBlZGdlOicsIG5vZGVFZGdlcy5ub2Rlc1swXSk7XG4gIFxuICBcdC8vIExpbWl0IG51bWJlciBvZiB0aGUgdG90YWwgbm9kZXMgYnkgdGhlIGNvbmZpZ3VyYXRpb24gaWYgZXhjZWVkc1xuICBcdGlmIChub2RlRWRnZXMubm9kZXMubGVuZ3RoID4gY29uZmlnLm1heE51bU5vZGVzKXtcbiAgICBcdG5vZGVFZGdlcyA9IGZpbHRlckRhdGFCeVRvcE5vZGVzKG5vZGVFZGdlcywgY29uZmlnLm1heE51bU5vZGVzKTtcbiAgICAgIGNvbnNvbGUubG9nKCdTZWxlY3QgdG9wbiBub2RlcycpO1xuICAgIH1cbiAgXHRjb25zdCBncmFwaE9iaiA9IGNvbnN0cnVjdEdyYXBoKG5vZGVFZGdlcyk7XG5cbiAgICAvLyBUT0RPOiBkZXRlcm1pbmUgdGhlc2UgMyBwYXJhbWV0ZXJzIGR5bmFtaWNhbGx5XG4gICAgY29uc3QgbnVtV2Fsa3MgPSBjb25maWcubnVtV2Fsa3M7IFxuICAgIGNvbnN0IHdhbGtMZW5ndGggPSBjb25maWcud2Fsa0xlbmd0aDtcbiAgICBjb25zdCBkaW1lbnNpb25zID0gY29uZmlnLmRpbWVuc2lvbnM7XG4gIFx0Y29uc3QgbnVtQ2x1c3RlcnMgPSBjb25maWcubnVtQ2x1c3RlcnM7XG4gIFx0Y29uc3QgbWF4SXRlcmF0aW9ucyA9IGNvbmZpZy5tYXhJdGVyYXRpb25zO1xuICAgIGNvbnN0IGNvbG9yUGFsZXR0ZSA9IFsnI2Y5NDE0NCcsICcjZjM3MjJjJywgJyNmOTg0NGEnLCBcbiAgICAnI2Y5Yzc0ZicsICcjOTBiZTZkJywgJyM0M2FhOGInLCAnIzRkOTA4ZScsICcjNTc3NTkwJywgJyMyNzdkYTEnLCBcbiAgICAnIzBmOGI4ZCcsICcjYTgyMDFhJywgJyMwNmQ2YTAnLCAnIzU4NTEyMycsICcjZWVjMTcwJywgXG4gICAgJyNmMmE2NWEnLCAnI2Y1ODU0OScsICcjNzcyZjFhJywgJyM0YzU3NjAnLCAnIzkzYThhYycsICcjYTU5ZThjJyxcbiAgICAnIzY2NjM1YicsICcjMDAyNTAwJywgJyMwZmEzYjEnLCAnI2VkZGVhNCcsICcjZjdhMDcyJywgJyNmZjliNDInLCAnI2ZlNWY1NScsICcjZjBiNjdmJ107XG4gIFxuXHRcdC8vIFRPRE86IEltcHJvdmUgdGhpcyBwZXJmb3JtYW5jZVxuICAgIGNvbnN0IGVtYmVkZGluZ3MgPSBub2RlRW1iZWRkaW5nKGdyYXBoT2JqLCBudW1XYWxrcywgd2Fsa0xlbmd0aCwgZGltZW5zaW9ucywgMCk7XG4gIFxuICBcdGxldCBjbHVzdGVycztcbiAgXHQvLyBsZXQgY29tbXVuaXR5RGljdDtcbiAgXHRsZXQgbm9kZUNsdXN0ZXJEaWN0O1xuICBcbiAgXHRpZiAoY29uZmlnLmNvbW11bml0eUZpbGUubGVuZ3RoID4gMSl7XG4gICAgICBcdG5vZGVDbHVzdGVyRGljdCA9IHt9O1xuICAgIFx0XHQvLyBSZWFkIGNsdXN0ZXJzIGFuZCBjb21tdW5pdHlEaWN0IGZyb20gdGhlIGZpbGUgaW5zdGVhZFxuICAgICAgXHRjbHVzdGVycyA9IE9iamVjdC52YWx1ZXMoY29tbXVuaXR5RGljdCk7XG4gICAgICAgIG5vZGVDbHVzdGVyRGljdCA9IHt9O1xuICAgICAgICBjbHVzdGVycy5mb3JFYWNoKChjbHVzdGVyLCBpbmRleCkgPT4ge1xuICAgICAgICAgIGNsdXN0ZXIuZm9yRWFjaCgobm9kZSkgPT4ge1xuICAgICAgICAgICAgbm9kZUNsdXN0ZXJEaWN0W25vZGVdID0gaW5kZXg7XG4gICAgICAgICAgfSk7XG4gIFx0XHRcdH0pO1xuICAgIH1cbiAgXHRlbHNlIHtcbiAgICAgIGNvbnN0IGNsdXN0ZXJzT3V0cHV0ID0ga01lYW5zQ2x1c3RlcmluZyhlbWJlZGRpbmdzLCBudW1DbHVzdGVycywgbWF4SXRlcmF0aW9ucywgMCk7XG4gICAgICBjbHVzdGVycyA9IGNsdXN0ZXJzT3V0cHV0WzBdO1xuICAgICAgY29tbXVuaXR5RGljdCA9IGNsdXN0ZXJzT3V0cHV0WzFdO1xuICAgICAgbm9kZUNsdXN0ZXJEaWN0ID0gY2x1c3RlcnNPdXRwdXRbMl07XG4gICAgfTtcbiAgICBcdFxuICBcdC8vIENyZWF0ZSBub2RlcyBhbmQgZWRnZXMgb2YgY29tbXVuaXR5IGZyb20gdGhlIGNsdXN0ZXJzXG4gIFx0Y29uc3Qgb3V0cHV0ID0gY29uc3RHcmFwaEZyb21DbHVzdChjbHVzdGVycywgbm9kZUVkZ2VzLmxpbmtzLCBjb2xvclBhbGV0dGUpO1xuICBcdGNvbnN0IGNvbW1Ob2RlRWRnZXMgPSBvdXRwdXRbMF07XG4gIFx0Y29uc3QgZ3JvdXBDb2xvcnMgICA9IG91dHB1dFsxXTtcbiAgXG4gIFx0Y29uc3QgZ3JhcGhSZWxhdGlvbiA9IFtdO1xuICBcblx0XHRjb25zdCBjb25maWdOdW1FZGdlc0NvbnN0ID0gW2NvbmZpZy5udW1Ob2RlUGVyQ2x1c3RlciwgY29uZmlnLm51bUVkZ2VQZXJDbHVzdGVyLCBjb25maWcubnVtQWRqTm9kZVBlckNsdXN0ZXIsIGNvbmZpZy5udW1BZGpFZGdlUGVyQ2x1c3Rlcl1cbiAgICBcbiAgICBsZXQgaSA9IDA7XG4gICAgZm9yIChjb25zdCBba2V5LCB2YWx1ZV0gb2YgT2JqZWN0LmVudHJpZXMoY29tbXVuaXR5RGljdCkpIHtcbiAgICAgICAgZ3JhcGhSZWxhdGlvbi5wdXNoKHtwYXJlbnRJZDoga2V5LCBjaGlsZHJlbjpjb25zdFN1YkdGcm9tQ2x1c3QodmFsdWUsIG5vZGVFZGdlcy5saW5rcywgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgaSwgZ3JvdXBDb2xvcnMsIG5vZGVDbHVzdGVyRGljdCwgXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgY29uZmlnTnVtRWRnZXNDb25zdCl9KVxuICAgICAgaSArPSAxXG5cdFx0fTtcbiAgICBcbiAgICBcbiAgXHQvLyBDcmVhdGUgbm9kZXMgYW5kIGVkZ2VzIG9mIHNlY29uZCBsYXllcnMgZnJvbSB0aGUgY2x1c3RlcnNcbiAgXHQvLyBmb3IobGV0IGkgPSAwOyBpIDwgY2x1c3RlcnMubGVuZ3RoOyBpKyspIHtcbiAgXHQvLyBsZXQgY2x1cyA9IGNsdXN0ZXJzW2ldO1xuICBcdC8vIGxldCBwYXJlbnRJZCA9ICdDb20uJytpO1xuICBcdC8vIGdyYXBoUmVsYXRpb24ucHVzaCh7cGFyZW50SWQ6IHBhcmVudElkLCBjaGlsZHJlbjpjb25zdFN1YkdGcm9tQ2x1c3QoY2x1cywgbm9kZUVkZ2VzLmxpbmtzLCBcbiAgXHQvLyBpLCBncm91cENvbG9ycywgbm9kZUNsdXN0ZXJEaWN0LCBcbiAgXHQvLyBjb25maWdOdW1FZGdlc0NvbnN0KX0pXG4gIFx0Ly8gfTtcbiAgXHRjb25zb2xlLmxvZyhncmFwaFJlbGF0aW9uKTtcblxuICAgIGNvbnN0IGdyYXBoID0gY29tbU5vZGVFZGdlcztcblxuICAgIC8vIEltcGxlbWVudCBkcmlsbFVwIGZ1bmN0aW9uXG4gICAgbGV0IGN1cnJlbnRMZXZlbCA9IDA7XG4gICAgbGV0IGN1cnJlbnROb2RlO1xuXG4gICAgZnVuY3Rpb24gZHJpbGxVcCgpIHtcbiAgICAgIGNvbnRhaW5lci5zZWxlY3RBbGwoJyonKS5yZW1vdmUoKTtcbiAgICAgIGlmIChjdXJyZW50TGV2ZWwgPiAwKSB7XG4gICAgICAgIC8vIEdldCB0aGUgcHJldmlvdXMgZ3JhcGggZGF0YSBiYXNlZCBvbiB0aGUgY3VycmVudCBsZXZlbFxuICAgICAgICBjdXJyZW50TGV2ZWwtLTtcbiAgICAgICAgY3VycmVudE5vZGUgPSBnZXRMZXZlbERhdGEoY3VycmVudExldmVsKTtcbiAgICAgICAgcmVuZGVyR3JhcGgoY3VycmVudE5vZGUpO1xuICAgICAgfVxuICAgIH1cblxuICAgIGZ1bmN0aW9uIHJlc2V0KCkge1xuICAgICAgY29udGFpbmVyLnNlbGVjdEFsbCgnKicpLnJlbW92ZSgpO1xuICAgICAgY3VycmVudExldmVsID0gMDtcbiAgICAgIGN1cnJlbnROb2RlID0gZ2V0TGV2ZWxEYXRhKGN1cnJlbnRMZXZlbCk7XG4gICAgICByZW5kZXJHcmFwaChjdXJyZW50Tm9kZSk7XG4gICAgfVxuXG4gICAgLy8gQWRkIGRyaWxsLXVwIGJ1dHRvblxuICAgIGNvbnN0IGJ1dHRvbiA9IHN2Z1xuICAgICAgLmFwcGVuZCgnZycpXG4gICAgICAuYXR0cignY2xhc3MnLCAnYnV0dG9uJylcbiAgICAgIC5hdHRyKFxuICAgICAgICAndHJhbnNmb3JtJyxcbiAgICAgICAgJ3RyYW5zbGF0ZSgnICsgKHdpZHRoIC0gMTAwKSArICcsIDIwKSdcbiAgICAgICkgLy8gQWRqdXN0IHRoZSBwb3NpdGlvbiBhcyBuZWVkZWRcbiAgICAgIC5vbignY2xpY2snLCBkcmlsbFVwKTtcblxuICAgIGJ1dHRvblxuICAgICAgLmFwcGVuZCgncmVjdCcpXG4gICAgICAuYXR0cignd2lkdGgnLCA4MClcbiAgICAgIC5hdHRyKCdoZWlnaHQnLCAzMClcbiAgICAgIC5hdHRyKCdyeCcsIDUpXG4gICAgICAuYXR0cigncnknLCA1KVxuICAgICAgLnN0eWxlKCdmaWxsJywgJ2xpZ2h0Ymx1ZScpO1xuXG4gICAgYnV0dG9uXG4gICAgICAuYXBwZW5kKCd0ZXh0JylcbiAgICAgIC5hdHRyKCd4JywgNDApXG4gICAgICAuYXR0cigneScsIDIwKVxuICAgICAgLmF0dHIoJ3RleHQtYW5jaG9yJywgJ21pZGRsZScpXG4gICAgICAuYXR0cignYWxpZ25tZW50LWJhc2VsaW5lJywgJ21pZGRsZScpXG4gICAgICAudGV4dCgnPC0gQmFjaycpO1xuXG4gICAgLy8gQWRkIGFub3RoZXIgYnV0dG9uIGJlbG93IHRoZSBmaXJzdCBidXR0b25cbiAgICBjb25zdCBidXR0b24yID0gc3ZnXG4gICAgICAuYXBwZW5kKCdnJylcbiAgICAgIC5hdHRyKCdjbGFzcycsICdidXR0b24nKVxuICAgICAgLmF0dHIoXG4gICAgICAgICd0cmFuc2Zvcm0nLFxuICAgICAgICAndHJhbnNsYXRlKCcgKyAod2lkdGggLSAxMDApICsgJywgNjApJ1xuICAgICAgKVxuICAgICAgLm9uKCdjbGljaycsIHJlc2V0KTtcblxuICAgIGJ1dHRvbjJcbiAgICAgIC5hcHBlbmQoJ3JlY3QnKVxuICAgICAgLmF0dHIoJ3dpZHRoJywgODApXG4gICAgICAuYXR0cignaGVpZ2h0JywgMzApXG4gICAgICAuYXR0cigncngnLCA1KVxuICAgICAgLmF0dHIoJ3J5JywgNSlcbiAgICAgIC5zdHlsZSgnZmlsbCcsICdsaWdodGdyZWVuJyk7XG5cbiAgICBidXR0b24yXG4gICAgICAuYXBwZW5kKCd0ZXh0JylcbiAgICAgIC5hdHRyKCd4JywgNDApXG4gICAgICAuYXR0cigneScsIDIwKVxuICAgICAgLmF0dHIoJ3RleHQtYW5jaG9yJywgJ21pZGRsZScpXG4gICAgICAuYXR0cignYWxpZ25tZW50LWJhc2VsaW5lJywgJ21pZGRsZScpXG4gICAgICAudGV4dCgnUmVzZXQnKTtcblxuICAgIHVwZGF0ZUxldmVsRGF0YShjdXJyZW50TGV2ZWwsIGdyYXBoKTtcbiAgICByZW5kZXJHcmFwaChncmFwaCk7XG5cbiAgICBmdW5jdGlvbiByZW5kZXJHcmFwaChncmFwaCkge1xuICAgICAgY29udGFpbmVyLnNlbGVjdEFsbCgnLmxpbmsnKS5yZW1vdmUoKTtcbiAgICAgIGNvbnRhaW5lci5zZWxlY3RBbGwoJy5ub2RlJykucmVtb3ZlKCk7XG4gICAgICBjb250YWluZXIuc2VsZWN0QWxsKCcqJykucmVtb3ZlKCk7XG5cbiAgICAgIGNvbnN0IG5vZGVzID0gZ3JhcGgubm9kZXM7XG4gICAgICBjb25zdCBsaW5rcyA9IGdyYXBoLmxpbmtzO1xuXG4gICAgICBjb25zdCB4ID0gc2NhbGVPcmRpbmFsKCkucmFuZ2UoW1xuICAgICAgICAyMCxcbiAgICAgICAgd2lkdGggLSAyMCxcbiAgICAgIF0pO1xuXG4gICAgICAvLyBDcmVhdGUgdGhlIGZvcmNlIHNpbXVsYXRpb25cbiAgICAgIGNvbnN0IHNpbXVsYXRpb24gPSBmb3JjZVNpbXVsYXRpb24obm9kZXMpXG4gICAgICAgIC5mb3JjZShcbiAgICAgICAgICAnbGluaycsXG4gICAgICAgICAgZm9yY2VMaW5rKGxpbmtzKS5pZCgoZCkgPT4gZC5pZClcbiAgICAgICAgKVxuICAgICAgICAuZm9yY2UoXG4gICAgICAgICAgJ2NoYXJnZScsXG4gICAgICAgICAgZm9yY2VNYW55Qm9keSgpLnN0cmVuZ3RoKC0yMDApXG4gICAgICAgIClcbiAgICAgICAgLmZvcmNlKFxuICAgICAgICAgICdjZW50ZXInLFxuICAgICAgICAgIGZvcmNlQ2VudGVyKHdpZHRoIC8gMiwgaGVpZ2h0IC8gMilcbiAgICAgICAgKTtcblxuICAgICAgLy8gSW5pdGlhbGl6ZSB0aGUgcG9zaXRpb25zIG9mIHRoZSBub2Rlc1xuICAgICAgbm9kZXMuZm9yRWFjaCgobm9kZSkgPT4ge1xuICAgICAgICBub2RlLnggPVxuICAgICAgICAgIE1hdGgucmFuZG9tKCkgKiAod2lkdGggLSAyMCkgKyAxMDtcbiAgICAgICAgbm9kZS55ID1cbiAgICAgICAgICBNYXRoLnJhbmRvbSgpICogKGhlaWdodCAtIDIwKSArIDEwO1xuICAgICAgfSk7XG5cbiAgICAgIGNvbnN0IHN0cm9rZVdpZHRoU2NhbGUgPSBzY2FsZUxpbmVhcigpXG4gICAgICAgIC5kb21haW4oWzAsIG1heChsaW5rcywgKGQpID0+IGQudmFsdWUpXSkgLy8gRGVmaW5lIHRoZSBpbnB1dCBkb21haW5cbiAgICAgICAgLnJhbmdlKFswLjEsIDNdKTsgLy8gRGVmaW5lIHRoZSBvdXRwdXQgcmFuZ2UgZm9yIHN0cm9rZS13aWR0aFxuXG4gICAgICBjb25zdCBzdHJva2VPcGFjaXR5ID0gc2NhbGVMaW5lYXIoKVxuICAgICAgICAuZG9tYWluKFswLCBtYXgobGlua3MsIChkKSA9PiBkLnZhbHVlKV0pIC8vIERlZmluZSB0aGUgaW5wdXQgZG9tYWluXG4gICAgICAgIC5yYW5nZShbMC4xLCAwLjRdKTsgLy8gRGVmaW5lIHRoZSBvdXRwdXQgcmFuZ2UgZm9yIHN0cm9rZS13aWR0aFxuXG4gICAgICAvLyBDcmVhdGUgdGhlIGxpbmtzXG4gICAgICBjb25zdCBsaW5rID0gY29udGFpbmVyXG4gICAgICAgIC5zZWxlY3RBbGwoJy5saW5rJylcbiAgICAgICAgLmRhdGEobGlua3MpXG4gICAgICAgIC5qb2luKCdsaW5lJylcbiAgICAgICAgLmF0dHIoJ2NsYXNzJywgJ2xpbmsnKVxuICAgICAgICAuc3R5bGUoJ29wYWNpdHknLCAwLjQpXG4gICAgICAgIC5zdHlsZSgnc3Ryb2tlJywgJyNhYWEnKVxuICAgICAgICAuc3R5bGUoJ3N0cm9rZS13aWR0aCcsIChkKSA9PlxuICAgICAgICAgIHN0cm9rZVdpZHRoU2NhbGUoZC52YWx1ZSB8fCAwLjEpXG4gICAgICAgICk7XG5cbiAgICAgIC8vIEZ1bmN0aW9uIHRvIGdlbmVyYXRlIGEgcmFuZG9tIGNvbG9yXG4gICAgICBjb25zdCBnZXRSYW5kb21Db2xvciA9ICgpID0+IHtcbiAgICAgICAgY29uc3QgbGV0dGVycyA9ICcwMTIzNDU2Nzg5QUJDREVGJztcbiAgICAgICAgbGV0IGNvbG9yID0gJyMnO1xuICAgICAgICBmb3IgKGxldCBpID0gMDsgaSA8IDY7IGkrKykge1xuICAgICAgICAgIGNvbG9yICs9XG4gICAgICAgICAgICBsZXR0ZXJzW1xuICAgICAgICAgICAgICBNYXRoLmZsb29yKE1hdGgucmFuZG9tKCkgKiAxNilcbiAgICAgICAgICAgIF07XG4gICAgICAgIH1cbiAgICAgICAgcmV0dXJuIGNvbG9yO1xuICAgICAgfTtcblxuICAgICAgY29uc3QgZ3JvdXBDb2xvcnMgPSB7fTtcbiAgICAgIC8vIE1ha2UgYSBoYXNoIHRhYmxlIGZvciBlYWNoIGdyb3VwXG4gICAgICBjb25zdCBjb2xvclBhbGV0dGUgPSBbJyNmOTQxNDQnLCAnI2YzNzIyYycsICcjZjg5NjFlJywgJyNmOTg0NGEnLCBcbiAgICAgICcjZjljNzRmJywgJyM5MGJlNmQnLCAnIzQzYWE4YicsICcjNGQ5MDhlJywgJyM1Nzc1OTAnLCAnIzBmOGI4ZCcsICcjYTgyMDFhJywgXG4gICAgICAnI2VmNDc2ZicsICcjZmZkMTY2JywgJyMwNmQ2YTAnLCAnIzA3M2I0YycsICcjNTg1MTIzJywgJyNlZWMxNzAnLCAnI2YyYTY1YScsIFxuICAgICAgJyNmNTg1NDknLCAnIzc3MmYxYScsICcjNGM1NzYwJywgJyM5M2E4YWMnLCAnI2E1OWU4YycsXG4gICAgICAnIzY2NjM1YicsICcjMDAyNTAwJywgJyMwZmEzYjEnLCAnI2VkZGVhNCcsICcjZjdhMDcyJywgJyNmZjliNDInLCAnI2ZlNWY1NScsICcjZjBiNjdmJ107XG4gICAgICBub2Rlcy5mb3JFYWNoKGZ1bmN0aW9uIChub2RlLCBpbmRleCkge1xuICAgICAgICBpZiAoIShub2RlLmdyb3VwIGluIGdyb3VwQ29sb3JzKSkge1xuICAgICAgICAgIGNvbnN0IGNvbG9ySW5kZXggPSBpbmRleCAlIGNvbG9yUGFsZXR0ZS5sZW5ndGg7XG4gICAgICAgICAgY29uc3QgY29sb3IgPSBjb2xvclBhbGV0dGVbY29sb3JJbmRleF07XG4gICAgICAgICAgZ3JvdXBDb2xvcnNbXG4gICAgICAgICAgICBub2RlLmdyb3VwXG4gICAgICAgICAgXSA9IGNvbG9yO1xuICAgICAgICB9XG4gICAgICB9KTtcbiAgICAgIFxuICAgICAgY29uc3QgcmFkaXVzU2NhbGUgPSBzY2FsZUxpbmVhcigpXG4gICAgICAuZG9tYWluKFswLCBtYXgobm9kZXMsIChkKSA9PiBkLm51bUVsZSldKSAvLyBEZWZpbmUgdGhlIGlucHV0IGRvbWFpblxuICAgICAgLnJhbmdlKFs3LCAyNV0pOyAvLyBEZWZpbmUgdGhlIG91dHB1dCByYW5nZSBmb3Igc3Ryb2tlLXdpZHRoXG5cbiAgICAgIGNvbnN0IG5vZGUgPSBjb250YWluZXJcbiAgICAgICAgLnNlbGVjdEFsbCgnLm5vZGVHcm91cCcpXG4gICAgICAgIC5kYXRhKG5vZGVzKVxuICAgICAgICAuam9pbignZycpXG4gICAgICAgIC5hdHRyKCdjbGFzcycsICdub2RlJylcbiAgICAgICAgLm9uKCdjbGljaycsIGhhbmRsZUNsaWNrKTtcblxuICAgICAgbm9kZVxuICAgICAgICAuYXBwZW5kKCdjaXJjbGUnKVxuICAgICAgICAuYXR0cihcInJcIiwgZCA9PiByYWRpdXNTY2FsZShkLm51bUVsZSkpXG4gICAgICAgIC5hdHRyKCdjbGFzcycsICdub2RlQ2lyY2xlJylcbiAgICAgICAgLmF0dHIoXG4gICAgICAgICAgJ2ZpbGwnLFxuICAgICAgICAgIChkKSA9PiBkLmdyb3VwQ29sb3JcbiAgICAgICAgKTtcblxuICAgICAgbm9kZVxuICAgICAgICAuYXBwZW5kKCd0ZXh0JylcbiAgICAgICAgLmF0dHIoJ2NsYXNzJywgJ25vZGVOYW1lJylcbiAgICAgICAgLmF0dHIoJ2R4JywgMTApXG4gICAgICAgIC5hdHRyKCdkeScsIDUpXG4gICAgICAgIC5hdHRyKCdmb250LXNpemUnLCAxMClcbiAgICAgICAgLmF0dHIoJ2ZvbnQtZmFtaWx5JywgJ1NoYXJlIFRlY2gnKVxuICAgICAgICAudGV4dCgoZCkgPT4gZC5pZClcbiAgICAgICAgLnN0eWxlKCd2aXNpYmlsaXR5JywgKGQpID0+IHtcbiAgICAgICAgICBpZiAobm9kZXMubGVuZ3RoIDwgMjApIHtcbiAgICAgICAgICAgIHJldHVybiAndmlzaWJsZSc7XG4gICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHJldHVybiAnaGlkZGVuJztcbiAgICAgICAgICB9XG4gICAgICAgIH0pO1xuXG4gICAgICBsZXQgc2VsZWN0ZWROb2RlID0gbnVsbDtcblxuICAgICAgZnVuY3Rpb24gaGFuZGxlQ2xpY2soZXZlbnQsIGQpIHtcbiAgICAgICAgaWYgKHNlbGVjdGVkTm9kZSA9PT0gZCkge1xuICAgICAgICAgIC8vIFNlY29uZCBjbGljayBvbiB0aGUgc2FtZSBub2RlLCByZXNldCBzZWxlY3Rpb25cbiAgICAgICAgICBzZWxlY3RlZE5vZGUgPSBudWxsO1xuICAgICAgICAgIHNlbGVjdEFsbCgnLm5vZGVOYW1lJykuc3R5bGUoXG4gICAgICAgICAgICAndmlzaWJpbGl0eScsXG4gICAgICAgICAgICAoZCkgPT4ge1xuICAgICAgICAgICAgICBpZiAobm9kZXMubGVuZ3RoIDwgMTApIHtcbiAgICAgICAgICAgICAgICByZXR1cm4gJ3Zpc2libGUnO1xuICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIHJldHVybiAnaGlkZGVuJztcbiAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfVxuICAgICAgICAgICk7XG4gICAgICAgICAgcmVzZXREYXRhKCk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgLy8gRmlyc3QgY2xpY2sgb3IgY2xpY2sgb24gYSBkaWZmZXJlbnQgbm9kZVxuICAgICAgICAgIHNlbGVjdEFsbCgnLm5vZGVOYW1lJykuc3R5bGUoXG4gICAgICAgICAgICAndmlzaWJpbGl0eScsXG4gICAgICAgICAgICAnaGlkZGVuJ1xuICAgICAgICAgICk7XG4gICAgICAgICAgc2VsZWN0KHRoaXMpXG4gICAgICAgICAgICAuc2VsZWN0KCcubm9kZU5hbWUnKVxuICAgICAgICAgICAgLnN0eWxlKCd2aXNpYmlsaXR5JywgJ3Zpc2libGUnKTtcbiAgICAgICAgICBzZWxlY3RlZE5vZGUgPSBkO1xuXG4gICAgICAgICAgLy8gUmVuZGVyIGNoaWxkcmVuIG5vZGVzXG4gICAgICAgICAgZ3JhcGhSZWxhdGlvbi5mb3JFYWNoKChyZWwpID0+IHtcbiAgICAgICAgICAgIGxldCBjdXJyZW50Q2hyID0gcmVsLmNoaWxkcmVuO1xuICAgICAgICAgICAgaWYgKHJlbC5wYXJlbnRJZCA9PT0gZC5pZCkge1xuICAgICAgICAgICAgICBjdXJyZW50TGV2ZWwrKztcbiAgICAgICAgICAgICAgdXBkYXRlTGV2ZWxEYXRhKFxuICAgICAgICAgICAgICAgIGN1cnJlbnRMZXZlbCxcbiAgICAgICAgICAgICAgICBjdXJyZW50Q2hyXG4gICAgICAgICAgICAgICk7XG4gICAgICAgICAgICAgIGxldCBjdXJyZW50RGF0YSA9IGdldExldmVsRGF0YShcbiAgICAgICAgICAgICAgICBjdXJyZW50TGV2ZWxcbiAgICAgICAgICAgICAgKTtcblxuICAgICAgICAgICAgICByZW5kZXJHcmFwaChjdXJyZW50RGF0YSk7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIC8vIElmIHRoZXJlJ3Mgbm8gY2hpbGRyZW4gbm9kZXMsIGNsaWNrIHRvIGhpZ2hsaWdodCB0aGUgbm9kZSBhbmQgaXRzIGFkamFjZW50XG4gICAgICAgICAgICB1cGRhdGVEYXRhKGQpO1xuICAgICAgICAgICAgdXBkYXRlR3JhcGgoKTtcbiAgICAgICAgICB9KTtcbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgICBmdW5jdGlvbiByZXNldERhdGEoKSB7XG4gICAgICAgIC8vIFJlc2V0IHNlbGVjdGVkIG5vZGUgYW5kIHVwZGF0ZSBkYXRhXG4gICAgICAgIC8vIGN1cnJlbnRMZXZlbCA9IDA7XG4gICAgICAgIHNlbGVjdGVkTm9kZSA9IG51bGw7XG4gICAgICAgIHVwZGF0ZURhdGEoKTtcbiAgICAgICAgdXBkYXRlR3JhcGgoKTtcbiAgICAgIH1cblxuICAgICAgZnVuY3Rpb24gdXBkYXRlRGF0YShzZWxlY3RlZE5vZGUpIHtcbiAgICAgICAgLy8gUmVzZXQgb3BhY2l0eSBmb3IgYWxsIG5vZGVzIGFuZCBsaW5rc1xuICAgICAgICBub2RlLnN0eWxlKCdvcGFjaXR5JywgMSk7XG4gICAgICAgIGxpbmsuc3R5bGUoJ29wYWNpdHknLCAwLjQpO1xuXG4gICAgICAgIGlmIChzZWxlY3RlZE5vZGUpIHtcbiAgICAgICAgICBjb25zdCBjbGlja2VkTm9kZUlkID0gc2VsZWN0ZWROb2RlLmlkO1xuICAgICAgICAgIGNvbnN0IG5laWdoYm9yTm9kZXMgPSBuZXcgU2V0KCk7XG5cbiAgICAgICAgICAvLyBGaW5kIG5laWdoYm9yIG5vZGVzXG4gICAgICAgICAgbGlua3MuZm9yRWFjaCgobGluaykgPT4ge1xuICAgICAgICAgICAgaWYgKFxuICAgICAgICAgICAgICBsaW5rLnNvdXJjZS5pZCA9PT0gY2xpY2tlZE5vZGVJZFxuICAgICAgICAgICAgKSB7XG4gICAgICAgICAgICAgIG5laWdoYm9yTm9kZXMuYWRkKGxpbmsudGFyZ2V0LmlkKTtcbiAgICAgICAgICAgIH0gZWxzZSBpZiAoXG4gICAgICAgICAgICAgIGxpbmsudGFyZ2V0LmlkID09PSBjbGlja2VkTm9kZUlkXG4gICAgICAgICAgICApIHtcbiAgICAgICAgICAgICAgbmVpZ2hib3JOb2Rlcy5hZGQobGluay5zb3VyY2UuaWQpO1xuICAgICAgICAgICAgfVxuICAgICAgICAgIH0pO1xuXG4gICAgICAgICAgLy8gRmlsdGVyIG5vZGVzIGFuZCBsaW5rc1xuICAgICAgICAgIG5vZGUuc3R5bGUoJ29wYWNpdHknLCAobm9kZSkgPT5cbiAgICAgICAgICAgIG5vZGUuaWQgPT09IGNsaWNrZWROb2RlSWQgfHxcbiAgICAgICAgICAgIG5laWdoYm9yTm9kZXMuaGFzKG5vZGUuaWQpXG4gICAgICAgICAgICAgID8gMVxuICAgICAgICAgICAgICA6IDAuMVxuICAgICAgICAgICk7XG4gICAgICAgICAgbGluay5zdHlsZSgnb3BhY2l0eScsIChsaW5rKSA9PlxuICAgICAgICAgICAgbGluay5zb3VyY2UuaWQgPT09IGNsaWNrZWROb2RlSWQgfHxcbiAgICAgICAgICAgIGxpbmsudGFyZ2V0LmlkID09PSBjbGlja2VkTm9kZUlkXG4gICAgICAgICAgICAgID8gMVxuICAgICAgICAgICAgICA6IDAuMVxuICAgICAgICAgICk7XG4gICAgICAgICAgLy8gc2hvdyB0aGUgbGFiZWxzIG9mIG5laWdoYm9yc1xuICAgICAgICAgIG5vZGVcbiAgICAgICAgICAgIC5zZWxlY3QoJy5ub2RlTmFtZScpXG4gICAgICAgICAgICAuc3R5bGUoJ3Zpc2liaWxpdHknLCAobm9kZSkgPT5cbiAgICAgICAgICAgICAgbm9kZS5pZCA9PT0gY2xpY2tlZE5vZGVJZCB8fFxuICAgICAgICAgICAgICBuZWlnaGJvck5vZGVzLmhhcyhub2RlLmlkKVxuICAgICAgICAgICAgICAgID8gJ3Zpc2libGUnXG4gICAgICAgICAgICAgICAgOiAnaGlkZGVuJ1xuICAgICAgICAgICAgKTtcbiAgICAgICAgfVxuICAgICAgfVxuXG4gICAgICBmdW5jdGlvbiB1cGRhdGVHcmFwaCgpIHtcbiAgICAgICAgLy8gVXBkYXRlIHRoZSB2aXN1YWwgZWxlbWVudHMgb2YgdGhlIGdyYXBoXG4gICAgICAgIGxpbmtcbiAgICAgICAgICAuYXR0cigneDEnLCAoZCkgPT4gZC5zb3VyY2UueClcbiAgICAgICAgICAuYXR0cigneTEnLCAoZCkgPT4gZC5zb3VyY2UueSlcbiAgICAgICAgICAuYXR0cigneDInLCAoZCkgPT4gZC50YXJnZXQueClcbiAgICAgICAgICAuYXR0cigneTInLCAoZCkgPT4gZC50YXJnZXQueSlcbiAgICAgICAgICAuYXR0cignc3Ryb2tlJywgJ2JsYWNrJyk7XG5cbiAgICAgICAgbm9kZS5hdHRyKFxuICAgICAgICAgICd0cmFuc2Zvcm0nLFxuICAgICAgICAgIChkKSA9PiBgdHJhbnNsYXRlKCR7ZC54fSwgJHtkLnl9KWBcbiAgICAgICAgKTtcbiAgICAgIH1cblxuICAgICAgZnVuY3Rpb24gaGFuZGxlTW91c2VPdmVyKGV2ZW50LCBkKSB7XG4gICAgICAgIGlmIChzZWxlY3RlZE5vZGUgPT09IG51bGwpIHtcbiAgICAgICAgICBzZWxlY3QodGhpcylcbiAgICAgICAgICAgIC5zZWxlY3QoJy5ub2RlTmFtZScpXG4gICAgICAgICAgICAuc3R5bGUoJ3Zpc2liaWxpdHknLCAndmlzaWJsZScpO1xuICAgICAgICB9XG4gICAgICB9XG5cbiAgICAgIGZ1bmN0aW9uIGhhbmRsZU1vdXNlT3V0KGV2ZW50LCBkKSB7XG4gICAgICAgIGlmIChzZWxlY3RlZE5vZGUgPT09IG51bGwpIHtcbiAgICAgICAgICBzZWxlY3QodGhpcylcbiAgICAgICAgICAgIC5zZWxlY3QoJy5ub2RlTmFtZScpXG4gICAgICAgICAgICAuc3R5bGUoJ3Zpc2liaWxpdHknLCAnaGlkZGVuJyk7XG4gICAgICAgIH1cbiAgICAgIH1cblxuICAgICAgLy8gQWRkIGRyYWcgYmVoYXZpb3IgdG8gdGhlIG5vZGVzXG4gICAgICBub2RlLmNhbGwoXG4gICAgICAgIGRyYWcoKVxuICAgICAgICAgIC5vbignc3RhcnQnLCBkcmFnU3RhcnQpXG4gICAgICAgICAgLm9uKCdkcmFnJywgZHJhZ2dpbmcpXG4gICAgICAgICAgLm9uKCdlbmQnLCBkcmFnRW5kKVxuICAgICAgKTtcblxuICAgICAgLy8gRGVmaW5lIHRoZSBkcmFnIGZ1bmN0aW9uc1xuICAgICAgZnVuY3Rpb24gZHJhZ1N0YXJ0KGV2ZW50LCBkKSB7XG4gICAgICAgIGlmICghZXZlbnQuYWN0aXZlKVxuICAgICAgICAgIHNpbXVsYXRpb24uYWxwaGFUYXJnZXQoMC4zKS5yZXN0YXJ0KCk7XG4gICAgICAgIGQuZnggPSBkLng7XG4gICAgICAgIGQuZnkgPSBkLnk7XG4gICAgICB9XG5cbiAgICAgIGZ1bmN0aW9uIGRyYWdnaW5nKGV2ZW50LCBkKSB7XG4gICAgICAgIGQuZnggPSBldmVudC54O1xuICAgICAgICBkLmZ5ID0gZXZlbnQueTtcbiAgICAgIH1cblxuICAgICAgZnVuY3Rpb24gZHJhZ0VuZChldmVudCwgZCkge1xuICAgICAgICBpZiAoIWV2ZW50LmFjdGl2ZSlcbiAgICAgICAgICBzaW11bGF0aW9uLmFscGhhVGFyZ2V0KDApO1xuICAgICAgICBkLmZ4ID0gbnVsbDtcbiAgICAgICAgZC5meSA9IG51bGw7XG4gICAgICB9XG5cbiAgICAgIC8vIFVwZGF0ZSB0aGUgbm9kZSBhbmQgbGluayBwb3NpdGlvbnMgb24gZWFjaCB0aWNrIG9mIHRoZSBzaW11bGF0aW9uXG4gICAgICBzaW11bGF0aW9uLm9uKCd0aWNrJywgdXBkYXRlR3JhcGgpO1xuICAgIH1cbiAgfSlcbiAgLmNhdGNoKChlcnJvcikgPT4ge1xuICAgIGNvbnNvbGUubG9nKCdFcnJvcjonLCBlcnJvcik7XG4gIH0pO1xuIl0sIm5hbWVzIjpbInNlbGVjdCIsImQzWm9vbSIsInNjYWxlT3JkaW5hbCIsImZvcmNlU2ltdWxhdGlvbiIsImZvcmNlTGluayIsImZvcmNlTWFueUJvZHkiLCJmb3JjZUNlbnRlciIsInNjYWxlTGluZWFyIiwibWF4Iiwic2VsZWN0QWxsIiwiZHJhZyJdLCJtYXBwaW5ncyI6Ijs7O0VBQUE7QUFvRkE7RUFDTyxTQUFTLHFCQUFxQjtFQUNyQyxFQUFFLE9BQU87RUFDVCxFQUFFLE9BQU8sQ0FBQyxJQUFJO0VBQ2QsRUFBRSxZQUFZLENBQUMsSUFBSTtFQUNuQixFQUFFO0VBQ0YsRUFBRSxNQUFNLEtBQUssR0FBRztFQUNoQixJQUFJLEtBQUssRUFBRSxFQUFFO0VBQ2IsSUFBSSxLQUFLLEVBQUUsRUFBRTtFQUNiLEdBQUcsQ0FBQztFQUNKLEVBQUUsTUFBTSxjQUFjLEdBQUcsSUFBSSxHQUFHLEVBQUUsQ0FBQztBQUNuQztFQUNBLEVBQUUsTUFBTSxLQUFLLEdBQUcsT0FBTyxDQUFDLEtBQUssQ0FBQyxJQUFJLENBQUMsQ0FBQztBQUNwQztFQUNBLEVBQUUsS0FBSyxNQUFNLElBQUksSUFBSSxLQUFLLEVBQUU7RUFDNUIsSUFBSSxNQUFNLENBQUMsS0FBSyxFQUFFLEtBQUssQ0FBQyxHQUFHLElBQUksQ0FBQyxLQUFLLENBQUMsUUFBUSxDQUFDLENBQUM7RUFDaEQsSUFBSSxJQUFJLGFBQWEsQ0FBQyxJQUFJLENBQUMsSUFBSSxDQUFDLElBQUksYUFBYSxDQUFDLElBQUksQ0FBQyxJQUFJLENBQUMsRUFBRTtBQUM5RDtFQUNBLE1BQU0sSUFBSSxLQUFLLEtBQUssRUFBRSxJQUFJLEtBQUssS0FBSyxFQUFFLEVBQUU7RUFDeEMsUUFBUSxNQUFNLE9BQU8sR0FBRyxPQUFPLEdBQUcsT0FBTyxDQUFDLEtBQUssQ0FBQyxRQUFRLEVBQUUsQ0FBQyxHQUFHLEtBQUssQ0FBQyxRQUFRLEVBQUUsQ0FBQztFQUMvRSxRQUFRLE1BQU0sT0FBTyxHQUFHLE9BQU8sR0FBRyxPQUFPLENBQUMsS0FBSyxDQUFDLFFBQVEsRUFBRSxDQUFDLEdBQUcsS0FBSyxDQUFDLFFBQVEsRUFBRSxDQUFDO0FBQy9FO0VBQ0EsUUFBUSxJQUFJLFVBQVUsR0FBRyxLQUFLLENBQUMsS0FBSyxDQUFDLFNBQVMsQ0FBQyxDQUFDLElBQUksS0FBSyxJQUFJLENBQUMsRUFBRSxLQUFLLE9BQU8sQ0FBQyxDQUFDO0VBQzlFLFFBQVEsSUFBSSxVQUFVLEtBQUssQ0FBQyxDQUFDLEVBQUU7RUFDL0IsVUFBVSxNQUFNLFNBQVMsR0FBRyxFQUFFLEVBQUUsRUFBRSxPQUFPLEVBQUUsQ0FBQztFQUM1QyxVQUFVLElBQUksWUFBWSxJQUFJLFlBQVksQ0FBQyxPQUFPLENBQUMsRUFBRTtFQUNyRCxZQUFZLFNBQVMsQ0FBQyxLQUFLLEdBQUcsWUFBWSxDQUFDLE9BQU8sQ0FBQyxDQUFDO0VBQ3BELFdBQVc7RUFDWCxVQUFVLEtBQUssQ0FBQyxLQUFLLENBQUMsSUFBSSxDQUFDLFNBQVMsQ0FBQyxDQUFDO0VBQ3RDLFNBQVM7QUFDVDtFQUNBLFFBQVEsSUFBSSxVQUFVLEdBQUcsS0FBSyxDQUFDLEtBQUssQ0FBQyxTQUFTLENBQUMsQ0FBQyxJQUFJLEtBQUssSUFBSSxDQUFDLEVBQUUsS0FBSyxPQUFPLENBQUMsQ0FBQztFQUM5RSxRQUFRLElBQUksVUFBVSxLQUFLLENBQUMsQ0FBQyxFQUFFO0VBQy9CLFVBQVUsTUFBTSxTQUFTLEdBQUcsRUFBRSxFQUFFLEVBQUUsT0FBTyxFQUFFLENBQUM7RUFDNUMsVUFBVSxJQUFJLFlBQVksSUFBSSxZQUFZLENBQUMsT0FBTyxDQUFDLEVBQUU7RUFDckQsWUFBWSxTQUFTLENBQUMsS0FBSyxHQUFHLFlBQVksQ0FBQyxPQUFPLENBQUMsQ0FB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