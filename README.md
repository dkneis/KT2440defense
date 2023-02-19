# Content

This repository contains observed data and model source code for the publication
**"Shifts from cooperative to individual-based predation defense determine microbial predator-prey dynamics"** by
Magali de la Cruz Barron, Ellen van Velzen, Uli Klümper, Markus Weitere, Thomas U. Berendonk, and David Kneis,
published in the *ISME Journal*, 2023, [DOI: 10.1038/s41396-023-01381-5](htpps://doi.org/10.1038/s41396-023-01381-5).

Here, we provide the original observation data together with the source code of
the ODE-based mathematical model presented in the publication. For any additional
information, please contact the corresponding author, David Kneis (firstname.lastname@tu-dresden.de).

## Observation data

We monitored the abundance of *Pseudomonas putida* (strain KT2440)
and *Poteriospumella lacustris* (strain JBM10) in semicontinuous co-cultures.
Generally, we observed an initial, toxin-based bacterial defense which resulted
in a marked decline of predator numbers. However, after a period of depression,
predators always recovered while a filamentous bacterial phenotype became
dominant in all replicates.

Data are provided as a tab-delimited text file ```data/observed.txt```. This file has
a typical database-like layout as shown below

```
day_zero    day  experiment  replicate_system  variable  value
2020-09-17  1    F0_10       a                 Bs        10300000
2020-09-17  1    F4_10       a                 Bs        6550000
2020-09-17  1    F5_10       a                 Bs        980000
[...]
```

with the columns

- *day_zero*: The day when the experiment started (YYYY-MM-DD)
- *day*: Number of the day after 'day_zero' when samples were taken (integer)
- *experiment*: Identifier for the experiment (string); see details below
- *replicate_system*: Identifier for the physical replicate (string) 
- *variable*: Identifier for the variable (string); see details below
- *value*: Observed value

The experiment identifiers are of the form "FX_YY" where X and YY are integers. The first number, X, specifies the initial density of flagellates (individuals/mL) in log10 units. The value of zero indicates a control experiment where predators were absent throughout. The second number, YY, indicates the dilution regime; specifically, it is the percentage of the culture which has been replaced by fresh medium in intervals of 24 hours. Example: The experiment identifier "F3_50" indicates that 50% of the culture was transferred to a new flask every 24 hours and that the initial flagellate density was 1000 ind/mL.

The variable identifiers are as follows:

- *F*: Flagellates (ind/mL)
- *Bs*: Single-cell bacteria (ind/mL)
- *Bf_num*: Number of bacterial filaments (objects/mL)
- *Bf_len*: Total length of filaments, summed over individuals (µm/mL)

Based on the measured volume/length of single cells, the data for the variable "Bf_len" can roughly be interpreted as single-cell equivalents (1 cell eq. = 1 µm of filament).

## ODE-based simulation model

### Prerequisites

The model was implemented in the [R language](https://www.r-project.org) and environment
for statistical computing. For solving simultaneous ODE numerically, we use the
[deSolve package](https://CRAN.R-project.org/package=deSolve) package. In addition, we make use of the [rodeo package](https://CRAN.R-project.org/package=rodeo)
written by the corresponding author. The latter serves two purposes:

1. It builds the system of differential equations from a description of the model
   in tabular form. See [this paper](https://doi.org/10.1016/j.envsoft.2017.06.036)
   for rationals and examples.
   
2. It generates [deSolve](https://CRAN.R-project.org/package=deSolve) compliant
   source code in modern Fortran and triggers compilation of the latter. Using
   compiled code speeds up numerical integration by orders of magnitude.

Thus, to execute the model source code provided here, you need the following free, open-source software:

- A recent installation of [R](https://www.r-project.org)

- The add-on package [deSolve](https://CRAN.R-project.org/package=deSolve)

- The add-on package [rodeo](https://CRAN.R-project.org/package=rodeo)

- A Fortran compiler recognized by your R installation. On Linux-based systems,
  you may want to install [gfortran](https://gcc.gnu.org/wiki/GFortran) from
  the GNU compiler collection. On Windows, you probably want to use the
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to have the compiler
  installed and set up.

### Files provided

#### Directory ```model/rodeo```

Contains the model description in tabular form. The extension ".tsv" indicates tab-separated text. The contents of these files is processed into Fortran source code by the rodeo package (see [this paper](https://doi.org/10.1016/j.envsoft.2017.06.036) for details).

- File ```variables.tsv```: Declares the state variables.

- File ```processesAndVariables.tsv```: Defines the ODEs right hand sides. They are
  constructed by multiplication of a vector of process rate expressions with the so-called
  stoichiometry matrix. The latter links processes and state variables.
  
- File ```parameters.tsv```: Declares biological parameters appearing in the ODE. Also contains the set of values to generate the published results.

- File ```functions.tsv```: Declares all functions appearing in the ODE.
  
- File ```functions.f95```: Fortran source code implementing the declared functions
  unless these are native Fortran function.

#### Files in the top-level directory ```model```

- File ```experiments.tsv```: Defines initial and boundary conditions for the set of experiments conducted. See above for an explanation of experiment identifiers. The example source code (see below) runs the model for all experimental setups defined in this file.

- File ```functions.R```: Implements two functions.

  1. The Function ```buildSingleStepModel``` constructs the model's equation system,
  generates the respective Fortran code, and subsequently triggers compilation into a
  shared library compliant with solvers from the ```deSolve``` package.
  
  2. The function ```multiStepSimulation``` makes it convenient to apply the model to a series of culture transfers where the values of state variables change instantaneously in response to dilution events. On every dilution event, integration of the ODE is paused and continued with updated initial conditions.

- File ```main.R```: This is the all-in-one script to build, compile, and run the model for all experiments. It also plots the dynamics of all state variables in a basic manner. You could run this script, for example, from within a dedicated IDE (if any is installed) or by typing ```source("main.R")``` at the R prompt or by running the command ```Rscript --vanilla main.R``` within a terminal. Since all used file paths are relative, your working directory should generally be the one where the file ```main.R``` resides (i.e. the ```model``` folder).
