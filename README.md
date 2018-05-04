# Project4
Ising Model of Ferromagnetism

### Nicolas Dronchi

## What is Ising Model; What is Monte Carlo Metrapolis

Ising Model: A lattice in 2d, 3d or 1d statement. For each grid point there is one particle and there are two spinning statement(up or down).
For the interactions between particles we only consider the nearest four grid points: up, low, left and right ones. (2D)

Monte Carlo is a stochastic simulation which can get an approximate value for a question through computer modeling.

Metrapolis rules a probablity for particles to have a balance spinning statement which is given by e^(-delta(T)/(kT).
k is the Boltzzman Constant

### Prerequisites

If you want to run the ising_visualization.py
```
Python 3.0
!pip install imageio --user
```

### Other Packages Installed

```
from math import exp
from random import choice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from copy import deepcopy
```

## This Repository
Pics -- Under pictures you will find pictures that were either created by matplotlib.pyplot or by myself. Most pictures are in the report but there are some extra including a gif visualizing the Ising model.

Report -- In this you will find both the pdf of the final report as well as the LaTeX file.

SavedVar -- These are numpy files for saved arrays from different runs of the ising model. These were created from the ising_main.py program. They are later used by ShowGraph.py

SourceFiles --
  * energyDistribution.py -- plots the distribution of energy for different tempeartures as well as showing the thermalization time
  * ising_main.py -- main program for running the ising model. This will save the arrays which is useful for long runs.
  * ising_visualization.py -- creates gifs of the evolution of the lattice after spin flips
  * ShawGraph.py -- Plots the variables saved from longer runs
  * small.py -- This is the 2x2 case for the ising model. It prints out the analytical and numerical values to compare and confirm the validity of the algorithm.


## Acknowledgments
*Thank you to my professor Morten Hjorth-Jensen. https://github.com/mhjensen

*Special thanks to: Joseph Slivka, Casey Chartier, and Zhiyang Yu my CMSE 202 final project team

*prtkm from Github(Dec. 2014) Monte-Carlo Simulations of the 2-D Ising Model.

              https://github.com/prtkm/ising-monte-carlo/blob/master/ising-monte-carlo.org
