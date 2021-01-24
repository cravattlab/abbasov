Chemoproteomics data processing and analysis for "A Proteome-Wide Atlas of Lysine-Reactive Chemistry" by Abbasov et. al

# Requirements:

- Python 3.7+ with pandas and numpy packages
- Rust programming language, with cargo package manager, which can be installed from <https://rustup.rs/> or your operating system's package manager, where applicable.

# Data structure:

Experiments (DTASelect, and CIMAGE combined output files) should be divided into 4 folders in `./data/`:
- `processed`: empty folder
- `abbasov_all`: all data sets
- `abbasov_noscouts`: partition the data sets into the next two folders, this one a copy of all data sets, except scouts
- `scouts`: the other partition of the data sets, these being only scout fragment treated experiments

Data files for downstream processing can then be generated by running `cargo run --release`.

Individual analyses (python scripts) are in the `figures` folder, and can be run after running the initial processing script
