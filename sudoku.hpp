#include <vector>

typedef std::vector<int> INT1;
typedef std::vector<INT1> INT2;
typedef std::vector<INT2> INT3;

void print_grid(const INT2& grid);
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
INT2 read_grid(char *path);

/**
 *  Fills obvious cells of the given grid and
 *  computes a candidate list for each cell
 *
 *  @return the candidate list
 */
INT3 preprocess(INT2& grid);

/**
 * Generates a pseudo-random grid based on the given canidates.
 *
 * @return the generated grid
 */
INT2 generate_individual(const INT3& candidates);

/**
 * Recombine two individuals to form a derived one.
 *
 * @return the derived individual
 */
INT2 crossover(const INT2* parents[2], const INT2& ref_grid);

/**
 * Alter an individual randomly
 */
void mutate(INT2& grid, const INT3& candidates);

/**
 * Calculate the fitness of a given grid
 *
 * @return the fitness value
 */
int fitness(const INT2& grid);