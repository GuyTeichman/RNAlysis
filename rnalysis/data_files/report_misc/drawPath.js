function drawPath(params) {
    // if something is selected:
    if (params.nodes.length > 0) {
        highlightActive = true;
        var selectedNode = params.nodes[0];
        var nodesToSelect = [selectedNode];
        var edgesToSelect = [];

        // mark all nodes as hard to read.
        for (var nodeId in allNodes) {
            allNodes[nodeId].color = "rgba(200,200,200,0.5)";
            allNodes[nodeId].width = 1;  // reset width
            if (allNodes[nodeId].hiddenLabel === undefined) {
                allNodes[nodeId].hiddenLabel = allNodes[nodeId].label;
                allNodes[nodeId].label = undefined;
            }
        }

        // Function to get all ancestors of a node
        function getAncestors(nodeId) {
            var ancestors = [];
            var queue = [nodeId];

            while (queue.length > 0) {
                var currentNode = queue.shift();
                var connectedEdges = network.getConnectedEdges(currentNode);
                var connectedNodes = network.getConnectedNodes(currentNode, 'from');

                for (var i = 0; i < connectedNodes.length; i++) {
                    if (!ancestors.includes(connectedNodes[i])) {
                        ancestors.push(connectedNodes[i]);
                        queue.push(connectedNodes[i]);
                        // Find the edge connecting these nodes and add to edgesToSelect
                        var edge = connectedEdges.find(edgeId => network.getConnectedNodes(edgeId).includes(connectedNodes[i]));
                        if (edge) edgesToSelect.push(edge);
                    }
                }
            }

            return ancestors;
        }

        // Get all ancestors of the selected node
        var ancestors = getAncestors(selectedNode);
        nodesToSelect = nodesToSelect.concat(ancestors);

        // Highlight all ancestors and the selected node
        for (var i = 0; i < nodesToSelect.length; i++) {
            var nodeId = nodesToSelect[i];
            allNodes[nodeId].color = nodeColors[nodeId];
            allNodes[nodeId].width = 2;  // make nodes bolder
            if (allNodes[nodeId].hiddenLabel !== undefined) {
                allNodes[nodeId].label = allNodes[nodeId].hiddenLabel;
                allNodes[nodeId].hiddenLabel = undefined;
            }
        }

    } else if (highlightActive === true) {
        // reset all nodes
        for (var nodeId in allNodes) {
            allNodes[nodeId].color = nodeColors[nodeId];
            allNodes[nodeId].width = 1;  // reset width
            if (allNodes[nodeId].hiddenLabel !== undefined) {
                allNodes[nodeId].label = allNodes[nodeId].hiddenLabel;
                allNodes[nodeId].hiddenLabel = undefined;
            }
        }
        highlightActive = false;
        nodesToSelect = [];
        edgesToSelect = [];
    }

    // Update the nodes in the network
    network.setData({nodes: new vis.DataSet(Object.values(allNodes)), edges: data.edges});

    // Set selection
    network.setSelection({nodes: nodesToSelect, edges: edgesToSelect}, {highlightEdges: false});
}
