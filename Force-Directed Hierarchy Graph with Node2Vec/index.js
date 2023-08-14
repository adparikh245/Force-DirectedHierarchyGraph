import {
  select,
  scaleLinear,
  forceSimulation,
  forceLink,
  forceManyBody,
  forceCenter,
  drag,
  max,
  zoom as d3Zoom,
  scaleOrdinal,
  selectAll,
  d3ZoomIdentity,
} from 'd3';
import {
  parseGraphDataset,
  updateLevelData,
  getLevelData,
  readText,
  parseGraphDatasetSNAP,
  parseGraphDatasetLabel,
  readClassId,
  nodeEmbedding,
  constructGraph,
  kMeansClustering,
  readTxt,
  constGraphFromClust,
  constSubGFromClust,
  filterDataByTopNodes,
  swapKeysAndValues,
  selectParseFormat
} from './graphUtils.js';
// import { readFileSync, writeFile, readFile } from 'fs'


const width = window.innerWidth;
const height = window.innerHeight;

// Create an SVG container for the graph
const svg = select('body')
  .append('svg')
  .attr('width', width)
  .attr('height', height);

// Create a <g> element as the container for the graph elements
var container = svg.append('g');

function zoomed(event) {
  container.attr('transform', event.transform);
}

// Create a zoom behavior
const zoom = d3Zoom()
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
			nodeEdges = selectParseFormat(nodesEdgesUrl, fileContent1)
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
    };
    	
  	// Create nodes and edges of community from the clusters
  	const output = constGraphFromClust(clusters, nodeEdges.links, colorPalette);
  	const commNodeEdges = output[0];
  	const groupColors   = output[1];
  
  	const graphRelation = [];
  
		const configNumEdgesConst = [config.numNodePerCluster, config.numEdgePerCluster, config.numAdjNodePerCluster, config.numAdjEdgePerCluster]
    
    let i = 0;
    for (const [key, value] of Object.entries(communityDict)) {
        graphRelation.push({parentId: key, children:constSubGFromClust(value, nodeEdges.links, 
                                                                      i, groupColors, nodeClusterDict, 
                                                                      configNumEdgesConst)})
      i += 1
		};
    
    
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

      const x = scaleOrdinal().range([
        20,
        width - 20,
      ]);

      // Create the force simulation
      const simulation = forceSimulation(nodes)
        .force(
          'link',
          forceLink(links).id((d) => d.id)
        )
        .force(
          'charge',
          forceManyBody().strength(-200)
        )
        .force(
          'center',
          forceCenter(width / 2, height / 2)
        );

      // Initialize the positions of the nodes
      nodes.forEach((node) => {
        node.x =
          Math.random() * (width - 20) + 10;
        node.y =
          Math.random() * (height - 20) + 10;
      });

      const strokeWidthScale = scaleLinear()
        .domain([0, max(links, (d) => d.value)]) // Define the input domain
        .range([0.1, 3]); // Define the output range for stroke-width

      const strokeOpacity = scaleLinear()
        .domain([0, max(links, (d) => d.value)]) // Define the input domain
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

      // Function to generate a random color
      const getRandomColor = () => {
        const letters = '0123456789ABCDEF';
        let color = '#';
        for (let i = 0; i < 6; i++) {
          color +=
            letters[
              Math.floor(Math.random() * 16)
            ];
        }
        return color;
      };

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
      
      const radiusScale = scaleLinear()
      .domain([0, max(nodes, (d) => d.numEle)]) // Define the input domain
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
          selectAll('.nodeName').style(
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
          selectAll('.nodeName').style(
            'visibility',
            'hidden'
          );
          select(this)
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

      function handleMouseOver(event, d) {
        if (selectedNode === null) {
          select(this)
            .select('.nodeName')
            .style('visibility', 'visible');
        }
      }

      function handleMouseOut(event, d) {
        if (selectedNode === null) {
          select(this)
            .select('.nodeName')
            .style('visibility', 'hidden');
        }
      }

      // Add drag behavior to the nodes
      node.call(
        drag()
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
