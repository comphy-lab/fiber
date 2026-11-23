# FIBER: Flow-Induced Breakup of Elastic Rheology

[![Website Status](https://img.shields.io/website?url=https%3A%2F%2Fcomphy-lab.github.io%2Ffiber%2F&style=flat-square&logo=github&label=Website)](https://comphy-lab.github.io/fiber/)
[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue?style=flat-square&logo=github)](https://github.com/comphy-lab/fiber)

---

🚀 **Fiber** is a state-of-the-art, open-source framework for simulating viscoelastic fluid flows in 2D and 3D, built as an extension to the [Basilisk C](http://basilisk.fr/) CFD library. It implements the log-conformation method for robust, high-Weissenberg number simulations, with a focus on clarity, extensibility, and scientific rigor.

## ✨ Key Features

- **Full 3D Log-Conformation Method**: Complete scalar implementation for 3D viscoelastic fluids ([log-conform-viscoelastic-scalar-3D.h](src-local/log-conform-viscoelastic-scalar-3D.h))
- **Robust 2D/Axi Support**: Scalar and tensor-based log-conformation for 2D and axisymmetric cases ([log-conform-viscoelastic-scalar-2D.h](src-local/log-conform-viscoelastic-scalar-2D.h), [log-conform-viscoelastic.h](src-local/log-conform-viscoelastic.h))
- **Optimized Matrix Algebra**: Efficient, stable eigenvalue and tensor operations ([eigen_decomposition.h](src-local/eigen_decomposition.h))
- **Advanced Initialization**: Functions for pseudo vectors/tensors in 2D/3D
- **Error Handling**: Negative eigenvalue checks, eigenvalue clamping, and detailed diagnostics
- **Performance**: Simplified acceleration term calculations and optimized tensor operations
- **Documentation**: Extensive inline documentation, mathematical background, and verification notes
- **Compatibility**: GPLv3 license, fully compatible with Basilisk and previous versions

## 🐛 Bug Fixes (v2.5/v2.6)
- Corrected matrix algebra in 3D
- Fixed rotation tensor and eigenvalue edge cases
- Improved error reporting and diagnostics
- Enhanced axisymmetric and 2D/3D compatibility

## 🗂️ Repository Structure

```
├── basilisk/src/               # Core Basilisk CFD library (reference only, do not modify)
├── src-local/                  # Custom viscoelastic solvers and tensor utilities
│   ├── log-conform-viscoelastic-scalar-2D.h   # 2D/axi log-conformation (scalar)
│   ├── two-phaseVE.h                          # Two-phase viscoelastic extension
├── simulationCases/            # Example/test cases and post-processing scripts
│   ├── filamentExt.c                            # Pinch-off of viscoelastic jet (2D/axi)
├── postProcess/                # Project-specific post-processing tools and utilities
│   ├── getData-elastic-scalar2D.c            # Data extraction utility
│   ├── getFacet2D.c                          # Facet extraction utility
│   └── VideoAxi.py                           # Python visualization script
```

## 📚 Documentation
- Inline documentation in all major headers (see `src-local/`)
- Example simulation and post-processing scripts in `simulationCases/`

## 🚀 Quick Start

### 1. Prerequisites
- [Basilisk C](http://basilisk.fr/) (use a trusted local installation; `basilisk/` is reference source)
- C compiler (e.g., gcc)
- Python 3 (for post-processing)
- Optional: Jupyter for notebooks

### 2. Compiling & Running Simulations

Configure a trusted, separately provisioned local Basilisk source directory. Do not use a legacy checked-in or pre-installed network viewer.

```bash
bash reset_install_requirements.sh /absolute/path/to/basilisk/src
source .project_config
```

Compile without `DISPLAY`:
```sh
qcc -O2 -Wall -I./src-local -disable-dimensions simulationCases/{CaseName}.c -o {CaseName} -lm
./{CaseName}
```

The network viewer is disabled in this checkout. Render saved simulation
outputs offline. Previously installed binaries or running processes are not changed.

Maintainers: update the website search index only through a reviewed pull request.
Verify the exact source revision, public content, and destination URLs before publishing;
do not mirror an unreviewed remote index or use branch-protection bypass credentials.

### 3. Post-Processing & Analysis
- Python scripts and Jupyter notebooks for data extraction and visualization are in `simulationCases/` (e.g., `VideoAxi.py`, `verifyWtihPlots.ipynb`).
- Example utilities: `getData-elastic-scalar2D.c`, `getFacet2D.c`.

## 📝 Example: Running a 3D Drop Atomisation Simulation

```sh
qcc -O2 -Wall -I./src-local -disable-dimensions simulationCases/dropAtomisation.c -o dropAtomisation -lm
./dropAtomisation
```

## 🔍 Technical Details
- **Log-Conformation Method**: See [src-local/log-conform-viscoelastic-scalar-3D.h](src-local/log-conform-viscoelastic-scalar-3D.h) and [src-local/log-conform-viscoelastic-scalar-2D.h](src-local/log-conform-viscoelastic-scalar-2D.h) for mathematical background and implementation notes.
- **Eigenvalue Solver**: [src-local/eigen_decomposition.h](src-local/eigen_decomposition.h) provides robust 3x3 symmetric eigensystem routines.
- **Two-Phase Flows**: [src-local/two-phaseVE.h](src-local/two-phaseVE.h) extends Basilisk's two-phase solver for viscoelasticity.
- **Axisymmetric/2D/3D**: Use the appropriate header for your geometry (see comments in each header for guidance).

## 🧑‍💻 Contributing
- See [CLAUDE.md](CLAUDE.md) for code style and development guidelines.
- Issue templates and feature requests: [GitHub Issue Templates](.github/ISSUE_TEMPLATE/)
- Pull requests are welcome! Please document your changes and update relevant tests/examples.

## 📋 License
This project is licensed under the [GNU GPLv3](LICENSE), in line with the Basilisk codebase.

## 🙏 Acknowledgments
- Thanks to all contributors and the Basilisk community

## 🔗 References
- Fattal & Kupferman (2004, 2005): Log-conformation method
- Comminal et al. (2015): Constitutive model functions
- Hao & Pan (2007): Split scheme implementation
- [Basilisk C](http://basilisk.fr/)

---
For detailed documentation, see the [docs/](docs/) folder or open `docs/index.html` in your browser.
