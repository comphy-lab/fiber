# FIBER Development Guidelines

## project structure 
├-- `basilisk/src/`: Core Basilisk CFD library (reference only, do not modify) <- This contains the codebase for Basilisk C--a wrapper around C to do computational fluid dynamics. 
## Some selected files inside basilisk/src:
│   ├── `navier-stokes/centered.h`: Main centered Navier-Stokes solver
│   ├── `navier-stokes/conserving.h`: Conservative form solver with VoF momentum advection
│   ├── `two-phase*.h`: Two-phase flow implementations (VoF/Level-set/CLSVOF)
│   ├── `reduced.h`: this implements the reduced gravity approach. 
│   ├── `curvature.h`: this implements the curvature module for computing interface properties.
│   ├── `tension.h`: this implements the surface tension force using the Brackbill method.
│   ├── `integral.h`: Integral formulation for surface tension.
│   ├── `tracer.h`: this implements the advection equation for passive tracers.
│   ├── `diffusion.h`: this implements the diffusion equation for passive tracers.
│   ├── `vof.h`: this implements the volume of fluid (VOF) method.
│   ├── `viscosity.h`: this implements the implicit viscous stress solver in Basilisk.
│   ├── `axi.h`: this implements axisymmetric metric so that 2D equations can be extended to axi. 
├── `postProcess/`: Project-specific post-processing tools
├── `src-local/`: Custom header files extending Basilisk functionality
├── `testCases/`: Test cases with their own Makefile

## Code Style
- Indentation: 2 spaces (no tabs)
- Line length: 80 characters max
- Documentation: Use markdown in comments starting with `/**`
- Spacing: Space after commas, spaces around operators (+, -)
- Files: Core functionality in `.h` headers, tests in `.c` files
- Naming: Snake_case for variables, camelCase for functions
- Error handling: Return values with stderr messages
- **Documentation**: Use detailed header comments with title, features, author, update history
- **Functions**: Document mathematical equations in comments before implementation
- **Comments**: Include TODO notes for future improvements

## Build & Test Commands
- Compile: `qcc -autolink file.c -o executable -lm`
- Compile with specific headers: `qcc -I$PWD/src-local -O2 -Wall -disable-dimensions file.c -o executable -lm`


## Best Practices
- Keep simulations modular and reusable
- Document physical assumptions and numerical methods
- Perform mesh refinement studies to ensure solution convergence
- Include visualization scripts in the postProcess directory


## Documentation Generation

- Read `.github/Website-generator-readme.md` for the website generation process.
- Do not auto-deploy the website; generating HTML is permitted using `.github/scripts/build.sh`.
- Avoid editing HTML files directly; they are generated using `.github/scripts/build.sh`, which utilizes `.github/scripts/generate_docs.py`.
- The website is deployed at `https://comphy-lab.org/repositoryName`; refer to the `CNAME` file for configuration. Update if not done already. 

### Directory Tree Formatting
- When creating or modifying directory trees in Markdown files:
  - Use proper tree characters: `├──`, `└──`, and `│   `
  - Maintain consistent 4-space indentation
  - Place trees inside code blocks (triple backticks)
  - Follow the pattern: `[tree chars]filename/    # description`
  - Use `└──` for the last item in each group
- Always refer to `.github/Website-generator-readme.md` for detailed formatting requirements
- Test the tree rendering by running the documentation generation script

## Purpose

This rule provides guidance for maintaining and generating documentation for code repositories in the CoMPhy Lab, ensuring consistency and proper workflow for website generation.

## Process Details

The documentation generation process utilizes Python scripts to convert source code files into HTML documentation. The process handles C/C++, Python, Shell, and Markdown files, generating a complete documentation website with navigation, search functionality, and code highlighting.

## Best Practices

- Always use the build script for generating documentation rather than manually editing HTML files
- Customize styling through CSS files in `.github/assets/css/`
- Modify functionality through JavaScript files in `.github/assets/js/`
- For template changes, edit `.github/assets/custom_template.html`
- Troubleshoot generation failures by checking error messages and verifying paths and dependencies