#include <vector>

typedef std::vector<int> A_INT1;
typedef std::vector<A_INT1> A_INT2;
typedef std::vector<A_INT2> A_INT3;

/** constant size of the currently loaded puzzle */
int SIZE;
/** constant size of each sub grid */
int SUB_SIZE;

struct Cell {
    A_INT1 raw_candidates;

    A_INT2 p_candidates;
    A_INT1 p_values;

    std::vector<Cell*> next;

    int init_value;

    Cell(int value);

    int value(int i);
};

typedef std::vector<Cell*> A_CELL1;
typedef std::vector<A_CELL1> A_CELL2;
typedef std::vector<A_CELL2> A_CELL3;


/**
 * Sets the next values of each cell in the grid accordingly
 */
void init_grid_successors(A_CELL2& grid);

/**
 * Propagate the choice of the given value to the given cell's successors
 * (for the given individual p)
 */
void propagate_choice(Cell* cell, int p, int value);

/**
 * Choose a random value from the cell's candidate list for the given individual
 *
 * Note: Value may be 0 if candidate list was empty
 *
 * @return the chosen value
 */
int place_random_value_into(Cell* cell, int p);

void print_grid(const A_CELL2& grid, int individual = -1);

/**
 * Reads a puzzle from a file located at path
 *
 * It also sets SIZE, SUM_SHOULD, PROD_SHOULD accordingly
 *
 * The file should be structured as follows:
 * n
 * i11 .. i1n
 * i21 .. i2n
 *  .  ..  .
 * in1 .. inn
 *
 * @return the grid
 */
A_CELL2 read_grid(char *path);

A_CELL2 deep_copy_grid(const A_CELL2& grid);


/**
 *  Fills obvious cells of the given grid and
 *  computes a candidate list for each cell
 */
void preprocess(A_CELL2& grid);

/**
 * Generates a pseudo-random grid based on the given reference. The individual is added to the population of the given grid
 *
 * @return the index of the generated individual
 */
int generate_individual(A_CELL2& grid);

/**
 * Alter an individual at position p randomly and place the result in buffer
 */
void mutate(A_CELL2& grid, A_CELL2& buffer, int p);

/**
 * Calculate the fitness of a given grid
 *
 * @return the fitness value
 */
int fitness(const A_CELL2& grid, int p);
