# 4.1.0 (2024-09-16)

Version 4.1.0 of RNAlysis features improved performance and stability across the board, thanks to a switch to Qt6 and
Numpy 2.

## Changed

- RNAlysis now runs on Qt6 instead of Qt5. This change should improve both performance and stability across RNAlysis.
- When importing parameters for differential expression analysis, RNAlysis will now automatically load the design matrix
  and chosen comparisons. This should work even with parameters that were exported in previous versions of RNAlysis.
- Made small improvements to the RNAlysis graphical interface.
- RNAlysis and its dependencies now run on Numpy 2 instead of Numpy 1.
- RNAlysis now uses a different implementation of the K-Medoids clustering algorithm, which should be more stable and
  faster than the previous implementation. However, note that the two implementations may give slightly different
  results.
- When running differential expression or feature counting, RNAlysis session reports will automatically include a
  logfile with R session info.
- Added optional parameters to all differential expression functions, allowing users to return a path to a logfile with
  R session info.

## Fixed

- Fixed bug where table previews in the graphical interface would sometimes fail to update when applying functions
  inplace.
