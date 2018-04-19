# Yin

Fast Python implementation of the Yin algorithm: a fundamental frequency estimator.

Based on the article:

\[1\] De Cheveigné, A., & Kawahara, H. (2002). YIN, a fundamental frequency estimator for speech and music. The Journal of the Acoustical Society of America, 111(4), 1917-1930.

All the functions in the code correspond to steps in the article \[1\]. Meanwhile, the difference function has been modify substantially in order to improve speed. Finally, speed has been improved by more than 1000x.


## Prerequisites

 * [Numpy](http://www.numpy.org/)
 * [Scipy](http://www.scipy.org/)
 * [Matlplotlib](http://matplotlib.org/) (for graphing)

## Usage

$python yin.py

All parameters (i.e frequence min, frequence max, harmonic threshold) in the yin.py function should be adapted to obtain good results. See the article \[1\] for more details. 


## Authors

Patrice Guyot

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1220947.svg)](https://doi.org/10.5281/zenodo.1220947)

If you use this code, please cite:
Patrice Guyot. (2018, April 19). Fast Python implementation of the Yin algorithm (Version v1.1.1). Zenodo. http://doi.org/10.5281/zenodo.1220947


Previous works on the implementation of the YIN algorithm have been made thanks to Robin Larvor, Maxime Le Coz and Lionel Koenig.
