# 4.0.0 (2024-06-29)

Version 4.0.0 of RNAlysis introduces two major improvements:
First, the software has switched from Pandas to Polars, which should substantially improve overall performance, particularly when generating automatic analysis reports.

Second, significant enhancements have been made to the automatic analysis reports.
These include new interactive features such as the ability to highlight analysis paths and filter the report by node types,
improved hierarchical layout for better readability, performance optimizations, and additional customization options.
The reports now offer a more comprehensive and user-friendly way to navigate and understand analysis workflows.

This update also includes various bug fixes and minor improvements to other functions within the software.

Happy analysis!

## Added

- Added additional optional parameters to automatic report generation.
- Added additional parameters to 'Hierarchical clustergram plot' (CountFilter.clustergram).

## Changed

- RNAlysis now uses Polars instead of Pandas. This change should improve the performance of RNAlysis by an order of magnitude, especially when generating automatic analysis reports.
- Automatic analysis reports are now laid out in a hierarchical structure, making it easier to navigate through the report.
- When clicking on a node in an automatic analysis report, the analysis path leading to that node will be highlighted.
- When clicking on a legend node in an automatic analysis report, all nodes of that type will be highlighted.
- The 'Hierarchical clustergram plot' function should now run faster on large datasets.

## Fixed

- Fixed bug where the visualization functions 'Plot histogram of p-values' and 'Histogram of a column' would not display a graph immediately when using the stand-alone app.
- Fixed bug re-loading a saved session report would sometimes lead to missing connections in the analysis flow graph.
- Fixed bug where applying not-inplace operations to differential expression tables would change the expected name of the p-value column.
- Fixed bug where table previews in automatic analysis reports would display an extra empty line.

## Backwards-incompatible changes

- Since RNAlysis now uses Polars instead of Pandas, all functions that previously returned Pandas DataFrames will now return Polars DataFrames. To keep your code compatible with the new version, you may need to adjust the way you interact with the returned DataFrames, or convert them back into Pandas DataFrames.
- Sessions saved with previous versions of RNAlysis are mostly compatible with RNAlysis 4.0.0, but some of the new features of automatic analysis reports (such as highlighting the report by node type) may not work in older sessions.
