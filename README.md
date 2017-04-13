# mbbp
Implementation of reduced restart tabu search for mbbp problem.
Support instances include [random instances](instances/GraphU_250_0.05_1.clq)  and [KONECT instances](http://konect.uni-koblenz.de/networks/)

## Installation
```
git clone  https://github.com/joey001/mbbp.git
make
```
## Usage
```
mbbp -f <filename> -t <maximum seconds> [-s <seed>] [-o <best know result>]
```
## Branches
* TSGR-MBBP: An implementation of TSGR-MBBP, which combines an Constraint-Balanced Tabu Search (CBTS) and two graph reduction techniques.
* Greedy: An implementation of A. Al-Yamani et al's algorithm 	https://doi.org/10.1109/TCSI.2007.907875

## License
GNU General Public License
