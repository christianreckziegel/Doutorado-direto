# Merging TTrees in ROOT Files

This repository contains macros written in ROOT to merge TTrees from different directories into single TTrees in single directories, and instructions on how to use them.

## Files

- **mergeTTrees.C**: This macro merges TTrees from different directories into single TTrees in single directories.
- **Merged_a.root**: Output file generated by `mergeTTrees.C`, containing merged TTrees from input files.
- **Merged_b.root**: Another output file generated by `mergeTTrees.C`, containing merged TTrees from input files.

## Usage

### Merging TTrees

1. Place your input ROOT files in the same directory as `mergeTTrees.C`.
2. Open a ROOT session.
3. Load the macro `mergeTTrees.C`.
4. Run the `mergeTTrees()` function.

```cpp
root [0] .L mergeTTrees.C
root [1] mergeTTrees()
```

This will merge TTrees from different ROOT directories in the input files and create `Merged_a.root` and `Merged_b.root` in the same directory as their respective ROOT input file.

### Merging Output Files

To merge the output files `Merged_a.root` and `Merged_b.root` into a single file called `AO2D_merged.root`, use the `hadd` command in the terminal or execute the macro `mergeFiles.C`.

```bash
hadd AO2D_merged.root Merged_a.root Merged_b.root
```
Or use the `mergeFiles.C` macro providing the correct merged files addresses.

## Dependencies

- ROOT: The macros are written in ROOT, a data analysis framework developed by CERN.

