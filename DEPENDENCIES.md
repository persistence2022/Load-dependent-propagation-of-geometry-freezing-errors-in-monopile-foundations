# Dependencies and scientific sources

MATLAB R2024b is the primary computation environment. MATLAB is a separate commercial product and is not included in this repository.

The optional Python comparison scripts import installed SciPy and OpenPile. Install them in a separate Python environment with `python -m pip install -r research_base/research/pisa_model/requirements_reference.txt`. This file records the dependency versions of the original comparison environment. Those optional Python comparisons are distinct from MATLAB execution.

OpenPile 1.0.2 is available at https://pypi.org/project/openpile/1.0.2/. Its own licence and notices apply to the separately installed package. This repository does not redistribute its source tree or wheel. The PISA curve comparison calls the installed OpenPile functions. The conventional sand modulus fit and capacity expressions are documented in OpenPile `utils.py_curves.api_sand` and `utils.hooks.InitialSubgradeReaction.api_sand`.

SciPy supplies the optional adaptive boundary-value solution through `scipy.integrate.solve_bvp`. That solver uses a separate discretization of the same constitutive equations. It verifies numerical implementation without supplying a physical ground truth. The associated SciPy reference is https://doi.org/10.1038/s41592-019-0686-2.

The PISA constitutive sources are Byrne et al. https://doi.org/10.1680/jgeot.18.P.255 and Burd et al. https://doi.org/10.1680/jgeot.18.P.277. Published numerical curves and interpreted soil profiles retain their source distinctions in DATA_SOURCES.csv. The field comparisons use published Cowden and Dunkirk data and are related to the PISA development lineage.

The conventional Matlock and API sand formulas are implemented as analytical expressions with explicit units. The Matlock core regularization is an implementation choice with separate numerical checks. The auxiliary layered soil profile is transcribed from the project manuscript and defines a controlled numerical scenario. Its underlying site identity and geotechnical-report provenance have not been independently verified.
