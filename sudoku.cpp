#include "sudoku.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>

// #define NDEBUG
#include <cassert>

// #define _WEIGHTED_FITNESS

std::random_device rd;
std::mt19937 rng(rd());

/** constant size of the currently loaded puzzle */
int SIZE;
/** constant size of each sub grid */
int SUB_SIZE;
/** constant expression each subset of cells evaluates to when accumulated, multiplied respectivley */
int SUM_SHOULD, PROD_SHOULD;

inline int best_fitness() {
#ifdef _WEIGHTED_FITNESS
    return 500;  // not true but is higher than any error possible
#else
    return 3 * SIZE;
#endif
}
int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
int dist2(int x1, int y1, int x2, int y2) {
    return (int)std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}
void print_usage(char *exec_name) {
    std::cout << "Usage: " << exec_name << " <path_to_sudoku> [<MAX_ITER> [<POP_SIZE> [<MUT_RATE>]]]" << std::endl;
}

int main(int argc, char **argv) {

    if (argc < 2) {
        print_usage(argv[0]);
        std::exit(1);
    }

    const int MAX_ITER = argc > 2 ? atoi(argv[2]) : 500;
    const int POP_SIZE = argc > 3 ? atoi(argv[3]) : 1000;
    const double MUT_RATE = argc > 4 ? atof(argv[4]) : 0.2;

    INT2 grid = read_grid(argv[1]);

    std::cout << "INPUT:" << std::endl;
    print_grid(grid);

    std::cout << "=========================" << std::endl;
    std::cout << "POP_SIZE: " << POP_SIZE << std::endl;
    std::cout << "MUT_RATE: " << MUT_RATE << std::endl;
    std::cout << "MAX_ITER: " << MAX_ITER << std::endl;
    std::cout << "=========================" << std::endl;

    INT3 candidates = preprocess(grid);

    // use two alternating buffers representing children and parents respectively
    INT3 pops[2] = { INT3(POP_SIZE), INT3(POP_SIZE) };
    INT3 *pop = &pops[0];
    INT3 *pop_buffer = &pops[1];

    INT1 fitness_values(POP_SIZE);
    std::vector<int>::iterator best;

    // Initialize population
    for (int i = 0; i < POP_SIZE; ++i) {
        (*pop)[i] = generate_individual(candidates);
    }

    for (int generation = 0; generation < MAX_ITER; ++generation) {
        // Calculate fitness for each individual
        for (int i = 0; i < POP_SIZE; ++i) {
            fitness_values[i] = fitness((*pop)[i]);
        }

        // Check if the best fitness is reached (i.e. the solution is found)
        best = std::max_element(fitness_values.begin(), fitness_values.end());
        if (*best == best_fitness()) {
            int idx = std::distance(fitness_values.begin(), best);
            std::cout << "SOLUTION FOUND!" << std::endl;
            print_grid((*pop)[idx]);
            std::exit(0);
        } else if ((generation+1) % 20 == 0) {
            std::cout << generation + 1 << "/" << MAX_ITER << "\t\tBest fitness: " << *best << std::endl;
        }

        std::discrete_distribution<int> d(fitness_values.begin(), fitness_values.end());
        std::bernoulli_distribution b(MUT_RATE);

        for (int i = 0; i < POP_SIZE; ++i) {
            const INT2 *parents[2] = { &(*pop)[d(rng)], &(*pop)[d(rng)] };
            (*pop_buffer)[i] = crossover(parents, grid);
            if (b(rng)) {
                mutate((*pop_buffer)[i], candidates);
            }
        }

        std::swap(pop, pop_buffer);
    }

    best = std::max_element(fitness_values.begin(), fitness_values.end());
    int idx = std::distance(fitness_values.begin(), best);
    std::cout << " -- Cancel at maximum iterations. --" << std::endl;
    std::cout << "Best fitness: " << *best << std::endl;
    print_grid((*pop)[idx]);

    return 0;
}


void print_grid(const INT2& grid) {
    std::cout << "  ---------------------" << std::endl;
    for (int li = 0; li < grid.size(); ++li) {
        auto line = grid[li];

        std::cout << "| ";
        for (int i = 0; i < line.size(); ++i) {
            if (line[i] == 0) {
                std::cout << " ";
            } else {
                std::cout << line[i];
            }

            if ((i + 1) % SUB_SIZE == 0) {
                std::cout << " | ";
            } else {
                std::cout << " ";
            }
        }

        std::cout << std::endl;

        if ((li + 1) % SUB_SIZE == 0) {
            std::cout << "  ---------------------";
        }
        std::cout << std::endl;
    }
}

INT2 read_grid(char *path) {
    std::ifstream stream(path);

    if (stream.is_open()) {
        stream >> SIZE;
        SUB_SIZE = (int)std::sqrt(SIZE);

        assert(SIZE == SUB_SIZE * SUB_SIZE);

        SUM_SHOULD = (SIZE + 1) * SIZE / 2;
        PROD_SHOULD = factorial(SIZE);

        INT2 rv(SIZE);
        for (int i = 0; i < SIZE; ++i) {
            rv[i] = INT1(SIZE);
            for (int j = 0; j < SIZE; ++j) {
                stream >> rv[i][j];
            }
        }

        stream.close();

        return rv;
    }

    std::cerr << "Could not open file!" << std::endl;
    std::exit(1);
}

INT3 preprocess(INT2& grid) {

    const int ALL = (1<<(SIZE+1)) - 1;

    INT2 bits(SIZE);
    for (int i = 0; i < SIZE; ++i) {
        bits[i] = INT1(SIZE, ALL);
    }

    // process all rows
    for (int r = 0; r < SIZE; ++r) {
        int mask = ALL;

        for (int i = 0; i < SIZE; ++i) {
            mask &= ~(1<<grid[r][i]);
        }

        for (int i = 0; i < SIZE; ++i) {
            bits[r][i] &= mask;
        }
    }

    // process all columns
    for (int c = 0; c < SIZE; ++c) {
        int mask = ALL;

        for (int i = 0; i < SIZE; ++i) {
            mask &= ~(1<<grid[i][c]);
        }

        for (int i = 0; i < SIZE; ++i) {
            bits[i][c] &= mask;
        }
    }

    // process all sub squares
    for (int r = 0; r < SIZE; r += SUB_SIZE) {
        for (int c = 0; c < SIZE; c += SUB_SIZE) {
            int mask = ALL;

            for (int i = 0; i < SUB_SIZE; ++i) {
                for (int j = 0; j < SUB_SIZE; ++j) {
                    mask &= ~(1<<grid[r+i][c+j]);
                }
            }

            for (int i = 0; i < SUB_SIZE; ++i) {
                for (int j = 0; j < SUB_SIZE; ++j) {
                    bits[r+i][c+j] &= mask;
                }
            }

        } // for c
    } // for r

    INT3 rv(SIZE);
    for (int i = 0; i < SIZE; ++i) {
        rv[i] = INT2(SIZE);
        for (int j = 0; j < SIZE; ++j) {
            if (grid[i][j] != 0) {
                rv[i][j] = INT1 { grid[i][j] };
                continue;
            }

            INT1 candidates;
            for (int n = 1; n <= SIZE; ++n) {
                if (bits[i][j] & (1<<n)) {
                    candidates.push_back(n);
                }
            }

            assert(candidates.size() > 0);

            if (candidates.size() == 1) {
                grid[i][j] = candidates[0];
            }

            rv[i][j] = candidates;
        }
    }

    return rv;
}

INT2 generate_individual(const INT3& candidates) {

    INT2 rv(SIZE);

    for (int i = 0; i < SIZE; ++i) {
        rv[i] = INT1(SIZE);

        for (int j = 0; j < SIZE; ++j) {
            std::uniform_int_distribution<int> random(0, candidates[i][j].size() - 1);
            rv[i][j] = candidates[i][j][random(rng)];
        }
    }

    return rv;
}

int fitness(const INT2& grid) {
    // interpret each constraint violation as binary error (i.e. 1 if constraint was violated)
    int value = 0;

    // check rows
    for (int r = 0; r < SIZE; ++r) {
        int prod = 1;
        int sum = 0;
        for (int i = 0; i < SIZE; ++i) {
            prod *= grid[r][i];
            sum += grid[r][i];
        }

        #ifdef _WEIGHTED_FITNESS
        value += (r + 1) * (int)(prod != PROD_SHOULD || sum != SUM_SHOULD);
        #else
        value += (int)(prod != PROD_SHOULD || sum != SUM_SHOULD);
        #endif
    }

    // check columns
    for (int c = 0; c < SIZE; ++c) {
        int prod = 1;
        int sum = 0;
        for (int i = 0; i < SIZE; ++i) {
            prod *= grid[i][c];
            sum += grid[i][c];
        }

        #ifdef _WEIGHTED_FITNESS
        value += (c + 1) * (int)(prod != PROD_SHOULD || sum != SUM_SHOULD);
        #else
        value += (int)(prod != PROD_SHOULD || sum != SUM_SHOULD);
        #endif
    }

    // check sub squares
    for (int r = 0; r < SIZE; r += SUB_SIZE) {
        for (int c = 0; c < SIZE; c += SUB_SIZE) {
            int prod = 1;
            int sum = 0;
            for (int i = 0; i < SUB_SIZE; ++i) {
                for (int j = 0; j < SUB_SIZE; ++j) {
                    prod *= grid[r+i][c+j];
                    sum += grid[r+i][c+j];
                }
            }

            #ifdef _WEIGHTED_FITNESS
            value += (dist2(0, 0, r, c) + 1) * (int)(prod != PROD_SHOULD || sum != SUM_SHOULD);
            #else
            value += (int)(prod != PROD_SHOULD || sum != SUM_SHOULD);
            #endif
        } // for c
    } // for r

    return best_fitness() - value;
}

INT2 crossover(const INT2* parents[2], const INT2& ref_grid) {
    // uniform crossover over two parents

    INT2 rv(SIZE);
    std::bernoulli_distribution b(0.5);

    for (int i = 0; i < SIZE; ++i) {
        rv[i] = INT1(SIZE);
        for (int j = 0; j < SIZE; ++j) {
            if (ref_grid[i][j] != 0) {
                rv[i][j] = ref_grid[i][j];
                continue;
            }

            rv[i][j] = (*parents[b(rng)])[i][j];
        }
    }

    return rv;
}

void mutate(INT2& grid, const INT3& candidates) {
    std::bernoulli_distribution b(0.5);

    for (int i = 0; i < SIZE; ++i) {
        for (int j = 0; j < SIZE; ++j) {
            if (b(rng)) {
                std::uniform_int_distribution<int> random(0, candidates[i][j].size() - 1);
                grid[i][j] = candidates[i][j][random(rng)];
            }
        }
    }
}