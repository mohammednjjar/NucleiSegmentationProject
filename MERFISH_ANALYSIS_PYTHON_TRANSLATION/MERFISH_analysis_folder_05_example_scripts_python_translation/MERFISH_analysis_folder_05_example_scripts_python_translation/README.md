# MERFISH_analysis Folder 05: `example_scripts/`

This package translates the MATLAB scripts under `example_scripts/`, including:

- `analysis_script.m`
- `code_construction_script.m`
- `library_design_example.m`
- `decoding/runMERFISH.m`
- `decoding/startup.m`
- `deprecated/ExampleLibraryConstruction_140genes_script.m`

The translated Python files preserve the MATLAB scripts' procedural flow and call the translated MERFISH modules by their corresponding Python names. Scripts that depend on later folders, such as `analysis/` and `deprecated/`, call those translated modules when available.

## Files

See `translation_manifest.csv` for every MATLAB file, Python file, purpose, and use.

## Checks

- Python compilation is checked for every `.py` file.
- Basic unit tests cover the standalone barcode demo and FASTA utilities.
- The scan report checks for unresolved scaffold markers.
