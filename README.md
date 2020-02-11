# surF - Surrogate Fourier modeling 

surF is a novel surrogate modeling technique that leverages the discrete Fourier transform to generate a smoother, and possibly easier to explore, fitness landscape.

## Usage

First of all, import surF as follows (please mind the upper case F):

```
from surfer import surF
```

Assume now that you have a fitness function ```f()``` defined over a search space ```hypercube```.

In order to build a surrogate model with surF, considering ```gamma``` Fourier coefficients, built with ```sigma``` samples of the fitness landscape and interpolated with a grid with ```rho``` steps, use the following code:

```
S = surF()
S.specify_fitness(fitness)
S.specify_search_space(hypercube)
S.build_model(coefficients=gamma, numpoints=sigma, resolution=rho)
```

Now, it is possible to exploit surF's  ```approximate(x)``` method to calculate the fitness value of a candidate solution ```x``` using the Fourier surrogate model.

## Citing surF

If you find surF useful for your research, please cite our work as follows:

Manzoni L., Papetti D.M., Cazzaniga P., Spolaor S., Mauri G., Besozzi D., and Nobile M.S.: Surfing on Fitness Landscapes: FST-PSO Powered by Fourier Surrogate Modeling (under revision)

## Additional information

For any information please contact: 
Luca Manzoni (luca.manzoni@units.it)
Marco S. Nobile (m.s.nobile@tue.nl)
