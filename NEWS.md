# NetworkLayout Release Notes

## v0.4.11 Changelog
- Fixes two bugs in the iteration scheme:
  - The layout iterator returned the second-to-last layout; now it returns the final positions after full iteration.
  - Iterative layouts computed `iterations` layouts, including the initial guess. Now `iterations` describes the number of actual *iterations*, and the layout iterator returns at most `iterations + 1` layouts.
