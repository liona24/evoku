# evoku
Solve sudokus using evolutionary algorithms

This is a simple demo project.

To compile simply run
```
$ g++ sudoku.cpp -o sudoku
```

Usage:
```
$ ./sudoku
Usage: ./sudoku <path_to_sudoku> [<MAX_ITER> [<POP_SIZE> [<MUT_RATE>]]]
```

Some demo puzzles are available in `puzzles/`. The easy ones are solved quite fast. Struggles become visible at medium difficulty and the program eventually fails at hard ones.

`sudoku2.cpp` is a bit more complex approach which tries to mimic a kind of stochastic backtracking solution. The optimization technique is quite naive though. It will generate a set of solutions and is then more likely to take solutions which are better (i.e. more cells filled).
Since the implementation logic got quite hard at times I chose a simple algorithm for exploring solutions which is very close to backtracking.
If one chooses a more sophisticated approach for recombination and selection I believe this greatly outperforms the other solution presented. However this way it is only able to solve the easy ones (more robust though).


Example output:
```
$ ./sudoku puzzles/easy1.txt
INPUT:
  ---------------------
|     1 | 7     | 4     |

|   5   |   9   |   8 6 |

| 4     | 6 2   | 3 7   |
  ---------------------
| 5     |     6 | 9     |

| 9 3   | 1   5 |   2 4 |

|     2 | 3     |     1 |
  ---------------------
|   9 6 |   1 2 |     3 |

| 1 4   |   3   |   6   |

|     3 |     4 | 7     |
  ---------------------
=========================
POP_SIZE: 1000
MUT_RATE: 0.2
MAX_ITER: 500
=========================
20/500		Best fitness: 16
40/500		Best fitness: 20
SOLUTION FOUND!
  ---------------------
| 3 6 1 | 7 5 8 | 4 9 2 |

| 2 5 7 | 4 9 3 | 1 8 6 |

| 4 8 9 | 6 2 1 | 3 7 5 |
  ---------------------
| 5 1 4 | 2 8 6 | 9 3 7 |

| 9 3 8 | 1 7 5 | 6 2 4 |

| 6 7 2 | 3 4 9 | 8 5 1 |
  ---------------------
| 7 9 6 | 8 1 2 | 5 4 3 |

| 1 4 5 | 9 3 7 | 2 6 8 |

| 8 2 3 | 5 6 4 | 7 1 9 |
  ---------------------
  ```
