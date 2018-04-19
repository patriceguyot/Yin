# Yin

Fast Python Implementation of the Yin algorithm: a fundamental frequency estimator.

Based on the article:

\[1\] De Cheveigné, A., & Kawahara, H. (2002). YIN, a fundamental frequency estimator for speech and music. The Journal of the Acoustical Society of America, 111(4), 1917-1930.

All the functions in the code correspond to step in the article \[1\]. Meanwhile, the difference function has been modify substantially in order to improve speed. Finally, speed have been improved by more than 1000x.


## Prerequisites

 * [Numpy](http://www.numpy.org/)
 * [Scipy](http://www.scipy.org/)
 * [Matlplotlib](http://matplotlib.org/) (for graphing)

## Usage

$python yin.py

All parameters (i.e frequence min, frequence max, harmonic threshold) in the yin.py function should be adapted to obtain good results.


## Authors

Patrice Guyot

Previous works on the implementation of the YIN algorithm have been made thanks to Robin Larvor, Maxime Le Coz and Lionel Koenig.
