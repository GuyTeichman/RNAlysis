// Add a one-time event listener for stabilization
network.once("stabilizationIterationsDone", function () {
    // Fit the network to view
    network.fit();

    // Disable further fitting
    network.setOptions({physics: {stabilization: {fit: false}}});
});
// use the "drawPath" event to redraw the graph with path to root highlighted
network.on("click", drawPath);
