#include "sudoku2.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>

// #define NDEBUG
#include <cassert>

std::random_device rd;
std::mt19937 rng(rd());

inline int best_fitness() {
    return SIZE * SIZE;
}
void print_usage(char *exec_name) {
    std::cout << "Usage: " << exec_name << " <path_to_sudoku> [<MAX_ITER> [<POP_SIZE>]]" << std::endl;
}

int main(int argc, char **argv) {

    if (argc < 2) {
        print_usage(argv[0]);
        std::exit(1);
    }

    const int MAX_ITER = argc > 2 ? atoi(argv[2]) : 500;
    const int POP_SIZE = argc > 3 ? atoi(argv[3]) : 1000;

    A_CELL2 grid = read_grid(argv[1]);

    std::cout << "INPUT:" << std::endl;
    print_grid(grid);

    std::cout << "=========================" << std::endl;
    std::cout << "POP_SIZE: " << POP_SIZE << std::endl;
    std::cout << "MAX_ITER: " << MAX_ITER << std::endl;
    std::cout << "=========================" << std::endl;

    preprocess(grid);

    // Initialize population
    for (int i = 0; i < POP_SIZE; ++i) {
        generate_individual(grid);
    }
    // initialize a buffer grid which is used to hold the new generation
    A_CELL2 next_gen = deep_copy_grid(grid);
    A_CELL2 grid_buffer = deep_copy_grid(grid);

    A_INT1 fitness_values(POP_SIZE*2);
    A_INT1::iterator best;

    for (int generation = 0; generation < MAX_ITER; ++generation) {
        // Calculate fitness for each individual
        for (int i = 0; i < POP_SIZE; ++i) {
            fitness_values[i] = fitness(grid, i);
        }

        // Check if the best fitness is reached (i.e. the solution is found)
        best = std::max_element(fitness_values.begin(), fitness_values.begin()+POP_SIZE);
        if (*best == best_fitness()) {
            int idx = std::distance(fitness_values.begin(), best);
            std::cout << "SOLUTION FOUND!" << std::endl;
            print_grid(grid, idx);
            std::exit(0);
        } else if ((generation+1) % 20 == 0) {
            std::cout << generation + 1 << "/" << MAX_ITER << "\t\tBest fitness: " << *best << std::endl;
        }

        for (int i = 0; i < POP_SIZE; ++i) {
            mutate(grid, next_gen, i);
            fitness_values[i+POP_SIZE] = fitness(next_gen, i);
        }
        std::discrete_distribution<int> d(fitness_values.begin(), fitness_values.end());
        Cell *tmp;
        int k = 0;
        for (int i = 0; i < POP_SIZE; ++i) {
            int p = d(rng);
            A_CELL2 *ref = p >= POP_SIZE ? &next_gen : &grid;
            if (p >= POP_SIZE) {
                p -= POP_SIZE;
            }
            for (int r = 0; r < SIZE; ++r) {
                for (int c = 0; c < SIZE; ++c) {
                    tmp = (*ref)[r][c];
                    grid_buffer[r][c]->p_candidates[k].assign(tmp->p_candidates[p].begin(), tmp->p_candidates[p].end());
                    grid_buffer[r][c]->p_values[k] = tmp->p_values[p];
                }
            }
            ++k;
        }

        std::swap(grid, grid_buffer);
    }

    best = std::max_element(fitness_values.begin(), fitness_values.end());
    int idx = std::distance(fitness_values.begin(), best);
    std::cout << " -- Cancel at maximum iterations. --" << std::endl;
    std::cout << "Best fitness: " << *best << std::endl;
    print_grid(grid, idx);

    return 0;
}

Cell::Cell(int value) {

    assert(value >= 0);

    this->init_value = value;
    if (value != 0) {
        this->raw_candidates = A_INT1 { value };
    }

    this->p_candidates = A_INT2();
    this->p_values = A_INT1();

    this->next = A_CELL1();
}

int Cell::value(int i) {
    if (i < 0) {
        return this->init_value;
    }

    assert(i < this->p_values.size());

    return this->p_values[i];
}

void print_grid(const A_CELL2& grid, int individual) {
    std::cout << "  ---------------------" << std::endl;
    for (int li = 0; li < grid.size(); ++li) {
        auto line = grid[li];

        std::cout << "| ";
        for (int i = 0; i < line.size(); ++i) {
            if (line[i]->value(individual) == 0) {
                std::cout << " ";
            } else {
                std::cout << line[i]->value(individual);
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

A_CELL2 read_grid(char *path) {
    std::ifstream stream(path);

    if (stream.is_open()) {
        stream >> SIZE;
        SUB_SIZE = (int)std::sqrt(SIZE);

        assert(SIZE == SUB_SIZE * SUB_SIZE);

        A_CELL2 rv(SIZE);
        for (int i = 0; i < SIZE; ++i) {
            rv[i] = A_CELL1(SIZE);
            for (int j = 0; j < SIZE; ++j) {
                int v; stream >> v;
                rv[i][j] = new Cell(v);
            }
        }

        stream.close();

        return rv;
    }

    std::cerr << "Could not open file!" << std::endl;
    std::exit(1);
}

void init_grid_successors(A_CELL2& grid) {
    // set connections between cells
    for (int i = 0; i < SIZE; ++i) {
        for (int j = 0; j < SIZE; ++j) {
            for (int r = i+1; r < SIZE; ++r) {
                grid[i][j]->next.push_back(grid[r][j]);
            }
            for (int c = j+1; c < SIZE; ++c) {
                grid[i][j]->next.push_back(grid[i][c]);
            }

            int block_r = i - (i % SUB_SIZE);
            int block_c = j - (j % SUB_SIZE);
            for (int r = i+1; r < block_r+SUB_SIZE; ++r) {
                for (int c = block_c; c < block_c+SUB_SIZE; ++c) {
                    if (c == j) {
                        continue;
                    }
                    grid[i][j]->next.push_back(grid[r][c]);
                }
            }
        } // for j
    } // for i
}

A_CELL2 deep_copy_grid(const A_CELL2& grid) {
    A_CELL2 rv(SIZE);
    Cell* tmp;
    for (int i = 0; i < SIZE; ++i) {
        rv[i] = A_CELL1(SIZE);
        for (int j = 0; j < SIZE; ++j) {
            tmp = grid[i][j];
            Cell* c = new Cell(tmp->init_value);
            c->raw_candidates.assign(tmp->raw_candidates.begin(), tmp->raw_candidates.end());
            c->p_candidates = A_INT2(tmp->p_candidates.size());
            c->p_values = A_INT1(tmp->p_values.size());
            for (int k = 0; k < c->p_values.size(); ++k) {
                c->p_candidates[k].assign(tmp->p_candidates[k].begin(), tmp->p_candidates[k].end());
                c->p_values[k] = tmp->p_values[k];
            }
            rv[i][j] = c;
        }
    }
    init_grid_successors(rv);

    return rv;
}

void preprocess(A_CELL2& grid) {

    const int ALL = (1<<(SIZE+1)) - 1;

    A_INT2 bits(SIZE);
    for (int i = 0; i < SIZE; ++i) {
        bits[i] = A_INT1(SIZE, ALL);
    }

    // process all rows
    for (int r = 0; r < SIZE; ++r) {
        int mask = ALL;

        for (int i = 0; i < SIZE; ++i) {
            mask &= ~(1<<grid[r][i]->value(-1));
        }

        for (int i = 0; i < SIZE; ++i) {
            bits[r][i] &= mask;
        }
    }

    // process all columns
    for (int c = 0; c < SIZE; ++c) {
        int mask = ALL;

        for (int i = 0; i < SIZE; ++i) {
            mask &= ~(1<<grid[i][c]->value(-1));
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
                    mask &= ~(1<<grid[r+i][c+j]->value(-1));
                }
            }

            for (int i = 0; i < SUB_SIZE; ++i) {
                for (int j = 0; j < SUB_SIZE; ++j) {
                    bits[r+i][c+j] &= mask;
                }
            }

        } // for c
    } // for r

    for (int i = 0; i < SIZE; ++i) {
        for (int j = 0; j < SIZE; ++j) {
            if (grid[i][j]->init_value != 0) {
                continue;
            }

            for (int n = 1; n <= SIZE; ++n) {
                if (bits[i][j] & (1<<n)) {
                    grid[i][j]->raw_candidates.push_back(n);
                }
            }
            assert(grid[i][j]->raw_candidates.size() > 0);
        }
    }

    init_grid_successors(grid);
}

void propagate_choice(Cell* cell, int p, int value) {
    for (int i = 0; i < cell->next.size(); ++i) {
        // check if this cell was visited yet
        if (cell->next[i]->p_values[p] != 0) {
            assert(false);
            continue;
        }
        A_INT1* n_c = &cell->next[i]->p_candidates[p];
        auto it = std::find(n_c->begin(), n_c->end(), value);

        if (it != n_c->end()) {
            n_c->erase(it);
        }
    }
}
int place_random_value_into(Cell* cell, int p) {
    // we want our sub solutions to be correct, meaning that we stop setting values if
    // we cannot set them correctly (i.e. we leave cells empty)
    int value = 0;
    if (cell->p_candidates[p].size() > 0) {
        // select random value at current position
        std::uniform_int_distribution<int> random(0, cell->p_candidates[p].size() - 1);
        value = cell->p_candidates[p][random(rng)];

        // remove this value from all other cells influenced by this choice
        propagate_choice(cell, p, value);
    }

    cell->p_values[p] = value;

    return value;
}

int generate_individual(A_CELL2& grid) {

    Cell* tmp;

    // initialize new candidate lists for the new individual
    for (int i = 0; i < SIZE; ++i) {
        for (int j = 0; j < SIZE; ++j) {
            tmp = grid[i][j];
            A_INT1 candidates(tmp->raw_candidates.size());
            std::copy(tmp->raw_candidates.begin(), tmp->raw_candidates.end(), candidates.begin());
            tmp->p_candidates.push_back(candidates);
            tmp->p_values.push_back(0);
        }
    }

    // index of our individual
    int p = grid[0][0]->p_candidates.size() - 1;

    for (int r = 0; r < SIZE; ++r) {
        for (int c = 0; c < SIZE; ++c) {
            place_random_value_into(grid[r][c], p);
        }
    }

    return p;
}

int fitness(const A_CELL2& grid, int p) {
    // since we enforce correctness all the time (by stopping to fill in values
    // which would violate the constrains) we only have to count the empty
    // cells to know if the solution was correct.

    int empty_count = 0;

    for (int i = 0; i < SIZE; ++i) {
        for (int j = 0; j < SIZE; ++j) {
            if (grid[i][j]->p_values[p] == 0) {
                ++empty_count;
            }
        }
    }

    return best_fitness() - empty_count;
}

void mutate(A_CELL2& grid, A_CELL2& buffer, int p) {
    std::uniform_int_distribution<int> random_pos(0, SIZE - 1);
    int init_r = random_pos(rng);
    int init_c = random_pos(rng);

    // first we have to find a cell which has a 'significant state'. This means
    // we have to find a cell which has more than one candidate available
    // We search backwards from our random starting point, if all fails we will
    // mutate the top-left most cell
    while (init_r != 0 || init_c != 0) {
        if (grid[init_r][init_c]->p_candidates[p].size() > 1) {
            break;
        }

        if (--init_c < 0) {
            --init_r;
            init_c = SIZE - 1;
        }
    }

    // in order to reconstruct our individual we just re-do the
    // construction process up to the position where we mutate it.
    // this simplifies the state managment A LOT and performance is still okay
    // since sudokus are not very large
    A_INT1 reroll_values(init_r*SIZE+init_c);
    int i, j;
    for (i = 0; i <= init_r; ++i) {
        int limit = i == init_r ? init_c : SIZE;
        for (j = 0; j < limit; ++j) {
            reroll_values[i*SIZE+j] = grid[i][j]->value(p);
        }
    }

    Cell* tmp;
    // Reset the individual
    for (i = 0; i < SIZE; ++i) {
        for (j = 0; j < SIZE; ++j) {
            tmp = buffer[i][j];
            tmp->p_candidates[p] = A_INT1(tmp->raw_candidates.size());
            std::copy(tmp->raw_candidates.begin(), tmp->raw_candidates.end(), tmp->p_candidates[p].begin());
            tmp->p_values[p] = 0;
        }
    }
    // Re-build it
    for (i = 0; i*SIZE < reroll_values.size(); ++i) {
        for (j = 0; j < SIZE && i*SIZE+j < reroll_values.size(); ++j) {
            tmp = buffer[i][j];
            int value = reroll_values[i*SIZE+j];
            propagate_choice(tmp, p, value);
            tmp->p_values[p] = value;
        }
    }
    --i;
    // Perform actual mutation
    for (; j < SIZE; ++j) {
        place_random_value_into(buffer[i][j], p);
    }
    ++i;
    for (; i < SIZE; ++i) {
        for (j = 0; j < SIZE; ++j) {
            place_random_value_into(buffer[i][j], p);
        }
    }
}