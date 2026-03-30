# README

## Overview

A simple tool to calculate Onsager coefficients from FCC AKMC simulations defined with pairwise nearest-neighbor interactions.

For overview of KMC simulations of RIS, see F. Soisson, "Kinetic Monte Carlo simulations of radiation induced segregation and precipitation." *Journal of Nuclear Materials*, 349 (3), pp. 235-250, 2006, https://doi.org/10.1016/j.jnucmat.2005.11.003. 

For definitions of KMC simulation inputs used in `settings.ini` and how they can be converted into pairwise interactions, see S. Shu, P. Bellon, and R. S. Averback, "Role of point-defect sinks on irradiation-induced compositional patterning in model binary alloys." *Physical Review B*, 91 (21), p. 214107, 2015, https://doi.org/10.1103/PhysRevB.91.214107

For the calculation methods used to obtain dilute limit Onsager coefficients, see J. D. Tucker, R. Najafabadi, T. R. Allen, and D. Morgan, "*Ab initio*-based diffusion theory and tracer diffusion in Ni-Cr and Ni-Fe alloys." *Journal of Nuclear Materials*, 405, pp. 216-234, 2010, https://doi.org/10.1016/j.jnucmat.2010.08.003

## Installation

#### Method 1) Installing via `uv` python package manager:

First, install necessary python packages:
```bash
uv sync
```

Run the script with:
```bash
uv run main.py
```

#### Method 2: Installing via `pip` python package manager:
First, create a virtual environment:
```bash
python -m venv onsager-calc
```

Activate the virtual environment with `source onsager-calc/Scripts/activate` or `source onsager-calc/bin/activate` depending on your OS. (You can deactivate the virtual environment with `deactive`.)


Then install necessary python packages:
```bash
pip install -r requirements.txt 
```

Run the script with 
```bash
python main.py
```

## Usage 
This script assumes one is working with an FCC binary alloy defined via nearest-neighbor pairwise interactions. 

Thermodynamic interactions, temperatures, compositions, geometry, etc. can be edited in `settings.ini`.

To calculate Onsager coefficients as a function of concentration or temperature, the code has an example overriding values in `settings.ini` with those taken from `example.csv` file. Edit the function `main()` in `main.py` and `read_csv()` in `io_handler.py` for your specific use case.

## Contributing
Pull requests and issues are welcome. When describing an issue, please include a minimal reproducible example, the expected behavior, and the actual behavior.

To make a contribution, please fork this repository an open a new branch. Make commits to the branch and open a pull request for review. In the pull request description, please explain your contribution and any issues you may have addressed.

## License
Copyright (c) 2026 nickplaysjazz [MIT License](https://mit-license.org/)

## TODO
- Add point defect profiles for grain boundary and bulk geometries. 

