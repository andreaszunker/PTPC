/*
 * Enumeration of Minimum Weight Codewords of Pre-Transformed Polar Codes by Tree Intersection
 * -------------------------------------------------------------------------------------------
 * Compile with: gcc -march=native -Ofast -o enumeration enumeration.c
 * 
 * Please cite the following reference if you want to use this algorithm in your 
 * research:
 * 
 * @INPROCEEDINGS{10480163,
 *   author={Zunker, Andreas and Geiselhart, Marvin and Ten Brink, Stephan},
 *   booktitle={2024 58th Annual Conference on Information Sciences and Systems (CISS)}, 
 *   title={Enumeration of Minimum Weight Codewords of Pre-Transformed Polar Codes by Tree Intersection}, 
 *   year={2024},
 *   doi={10.1109/CISS59072.2024.10480163}}
 */


// Includes
// --------
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


/*
 * Macro: MIN
 * ----------
 * Returns the minimum of 'a' and 'b'.
 */
#define MIN(a,b) ((a) < (b) ? (a) : (b))


/*
 * Macro: PRETRANSFORM
 * -------------------
 * Allows to use 2D array indexing of the 'pretransform' member of 'args'.
 */
#define PRETRANSFORM(args) \
    ((uint64_t(*)[(args).message_size])(args).pretransform)


/*
 * Struct: args
 * ------------
 */
struct args {
    size_t   coset_index, stop_level, message_size;
    uint8_t  *rate_profile, *sibling_levels;
    uint64_t *pretransform;
};


/*
 * Struct: result
 * --------------
 * Struct representing the result of "wmin"-weight codeword enumeration. It 
 * contains:
 * - wmin: Minimum weight "wmin"
 * - A_wmin: Number of "wmin"-weight codewords
 */
struct result {
  uint64_t wmin, A_wmin;
};


/*
 * Function: reed_muller_rate_profile
 * ----------------------------------
 * Constructs the binary representation of the rate-profile corresponding to the 
 * RM(r,m) code.
 *
 * r: Order of the Reed-Muller code
 * m: m = log2(N) of the code length N
 * rate_profile: Buffer of length N for the binary representation of the 
 *     rate-profile
 */
void reed_muller_rate_profile(int r, int m, uint8_t rate_profile[]) {
    for (size_t index = 0; index < (1 << m); ++index) {
        rate_profile[index] = __builtin_popcountl(index) >= m-r;
    }
}


/*
 * Function: convolutional_pretransform
 * ------------------------------------
 * Constructs the pre-transformation matrix corresponding to the given 
 * rate-profile and polynomial.
 *
 * K: Code dimension
 * N: Code length
 * pretransform: Buffer for the K×N pretransfom
 * rate_profile: Binary representation of the rate-profile
 * degree: Degree of the polynomial
 * polynomial: Polynomial in binary representation with least significant 
 *     coefficient first
 */
void convolutional_pretransform(
    size_t  K, 
    size_t  N, 
    uint8_t pretransform[K][N], 
    uint8_t rate_profile[N], 
    size_t  degree, 
    uint8_t polynomial[degree+1]
) {
    memset(pretransform, 0, sizeof(uint8_t[K][N]));
    size_t row = 0;
    for (size_t i = 0; i < N; ++i) {
        if (!rate_profile[i]) {
            continue;
        }
        for (size_t pos = 0, col = i; pos <= degree && col < N; ++pos, ++col) {
            pretransform[row][col] = polynomial[pos];
        }
        ++row;
    }
}


/*
 * Function: polar_transform
 * -------------------------
 * Applies the polar transform to the given matrix along the second dimension.
 * 
 * rows: Number of rows of the matrix
 * columns: Number of columns of the matrix
 * matrix: Matrix on which the polar transform is to be applied
 */
void polar_transform(
    size_t  rows, 
    size_t  columns, 
    uint8_t matrix[rows][columns]
) {
    for (size_t row = 0; row < rows; ++row) {
        for (size_t distance = 1; distance < columns; distance *= 2) { // Separation of the two inputs to be XORed
            for (size_t group = 0; group < columns; group += 2*distance) { // Group iterator
                for (size_t butterfly = 0; butterfly < distance; ++butterfly) { // Butterfly iterator
                     matrix[row][group+butterfly] ^= 
                         matrix[row][group+butterfly+distance];
                }
            }
        }
    } 
}


/*
 * Function: reduced_row_echelon_form
 * ----------------------------------
 * Brings the given matrix into reduced row echelon form (RREF).
 * 
 * rows: Number of rows of the matrix
 * columns: Number of columns of the matrix
 * matrix: Matrix that is to be brought into RREF
 */
void reduced_row_echelon_form(
    size_t  rows, 
    size_t  columns, 
    uint8_t matrix[rows][columns]
) {
    size_t current_row = 0, pivot_column = 0;
        
    while (current_row < rows && pivot_column < columns) {
        // Find the pivot element with a non-zero value in the current column
        size_t pivot_row = current_row;
        while (pivot_row < rows && matrix[pivot_row][pivot_column] == 0) {
            ++pivot_row;
        }
        if (pivot_row >= rows) {
            // No non-zero pivot element found in the current column
            ++pivot_column;
            continue;
        }
        if (current_row != pivot_row) {
            // Swap the pivot row with the current row
            for (size_t column = 0; column < columns; column++) {
                uint8_t temp = matrix[current_row][column];
                matrix[current_row][column] = matrix[pivot_row][column];
                matrix[pivot_row][column] = temp;
            }
        }
        // Eliminate all elements in the pivot column except the pivot itself
        for (size_t row = 0; row < rows; ++row) {
            if (row != current_row && matrix[row][pivot_column] != 0) {
                for (size_t column = pivot_column; column < columns; ++column) {
                    matrix[row][column] ^= matrix[current_row][column];
                }
            }
        }
        ++current_row; ++pivot_column;
    }
}


/*
 * Function: update_message
 * ------------------------
 * Updates the message to form a "wmin"-weight codeword of a universal polar 
 * coset.
 */ 
inline void update_message(int coset_index, int level, uint64_t message[]) {
    // Update the message according to the "M"-set formulation (using int instead of size_t makes a noticeable speed difference here).
    for (size_t i = coset_index+1; i < level; ++i) {
        if (message[i >> 6] & ((uint64_t)1 << (i & 63))
            && (~coset_index & level & i) == 0
        ) {
            int update_index = (~coset_index & (level | i)) | (level & i);
            message[update_index >> 6] ^= (uint64_t)1 << (update_index & 63);
        }
    }
    message[level >> 6] |= (uint64_t)1 << (level & 63); // set the "level"-th bit of the message to one. 
}


/*
 * Function: enumerate_subtree
 * ---------------------------
 * Counts the "wmin"-weight codewords contained in the given sub-tree.
 */ 
uint64_t enumerate_subtree(
    struct args args, 
    size_t      level, 
    uint64_t    start_message[]
) {
    // Copy the message so the original message can be used later to enumerate the other sub-tree
    uint64_t *message = malloc(sizeof(uint64_t[args.message_size])); 
    if (message == NULL) {
        exit(EXIT_FAILURE);
    }
    memcpy(message, start_message, sizeof(uint64_t[args.message_size]));
    
    uint64_t A_wmin = 0UL;
    update_message(args.coset_index, level, message);
    for (++level; level <= args.stop_level; ++level) {
        if (args.rate_profile[level]) { // Sibling level of the tree correspoding to the PTPC coset
            // Sibling level of the "wmin"-weight codeword tree of a universal polar coset -> sibling level of the intersection tree
            A_wmin += args.sibling_levels[level] ? enumerate_subtree(args, level, message) : 0;
            continue;
        }
        // Einzelchild level of the intersection tree
        uint64_t bit = 0L, message_bit = (message[level >> 6] >> (level & 63)) & 1; 
        for (size_t i = args.coset_index >> 6; i < ((level-1) >> 6)+1; ++i) {
            bit ^= message[i] & PRETRANSFORM(args)[level][i];
        }
        if ((__builtin_popcountl(bit) & 1) == message_bit) {
            continue;
        }
        // The pre-transformation does not match with the current message
        if (args.sibling_levels[level]) {
            // Sibling level of the "wmin"-weight codeword tree of a universal polar coset -> message can be updated
            update_message(args.coset_index, level, message);
            continue;
        }
        // Both trees have einzelchild levels -> the message cannot be adjusted -> the message path does not form a "wmin"-weight codeword
        free(message); 
        return A_wmin;
    }
    // The message path does form a "wmin"-weight codeword
    free(message); 
    return A_wmin + 1;
}


/*
 * Function: enumerate_minimum_weight_codewords
 * --------------------------------------------
 * Counts the "wmin"-weight codewords of the given pre-transformed polar code 
 * (PTPC).
 *
 * K: Code dimension
 * N: Code length
 * generator_matrix: K×N generator matrix
 *
 * Returns: A struct of type enumeration_result containing:
 * - wmin: Minimum weight "wmin"
 * - A_wmin: Number of "wmin"-weight codewords
 */ 
struct result enumerate_minimum_weight_codewords(
    size_t  K, 
    size_t  N, 
    uint8_t generator_matrix[K][N]
) {
    // Compute the pre-transformation matrix and bring it into RREF
    uint8_t (*pretransform)[N] = malloc(sizeof(uint8_t[K][N]));
    if (pretransform == NULL) {
        exit(EXIT_FAILURE);
    }
    memcpy(pretransform, generator_matrix, sizeof(uint8_t[K][N])); // copy so that the given generator matrix remains unchanged
    polar_transform(K, N, pretransform); // get the pre-transformation matrix
    reduced_row_echelon_form(K, N, pretransform); // bring the pre-transformation matrix into RREF 
    
    struct args args = {
        .message_size = ((N-1) >> 6)+1, // the message is stored bitwise -> faster checking of the dynamic frozen bits
        .rate_profile = calloc(N, sizeof(uint8_t)), // indicator of the sibling level of the PTPC 
        .pretransform = calloc(N*(((N-1) >> 6)+1), sizeof(uint64_t)),
        .sibling_levels = malloc(sizeof(uint8_t[N])), // indicator of the sibling level of the "wmin"-weight codeword tree of a universal polar coset
    };
    if (args.rate_profile == NULL || args.pretransform == NULL || args.sibling_levels == NULL) {
        exit(EXIT_FAILURE);
    } 
    
    // Expand the pre-transform into a N×N matrix and store it with bitwise columns and in Fortran order
    for (size_t row = 0, pivot_column = 0; row < K; ++row) {
        // Find the pivot point columns of the pre-transformation matrix -> information bits
        for (; !pretransform[row][pivot_column]; ++pivot_column);
        for (size_t column = 0; column < N; ++column) {
            PRETRANSFORM(args)[column][pivot_column >> 6] |= 
                (uint64_t)pretransform[row][column] << (pivot_column & 63);
        }
        args.rate_profile[pivot_column] = 1;  
    } 
    uint64_t *message = malloc(sizeof(uint64_t[args.message_size])); 
    if (message == NULL) {
        exit(EXIT_FAILURE);
    }      
    struct result result = {UINT64_MAX, 0UL};
    
    // Find minimum Hamming weight "wmin" of a coset leader
    for (size_t index = 0; index < N; ++index) {
        if (args.rate_profile[index]) {
            result.wmin = MIN(result.wmin, 1U << __builtin_popcount(index));
        }
    }
    
    // Count the "wmin"-weight codewords in each coset
    for (args.coset_index = 0; args.coset_index < N; ++args.coset_index) {
        uint64_t row_weight = 1UL << __builtin_popcount(args.coset_index);
        if (!args.rate_profile[args.coset_index] || row_weight > result.wmin) {
            continue; 
        }
        // The coset is led by a "wmin"-weight row
        memset(message, 0, sizeof(uint64_t[args.message_size])); 
        memset(args.sibling_levels, 0, sizeof(uint8_t[N])); 
        
        // Find the level after which the pre-transformation cannot prevent the formation of "wmin"-weight codewords and compute their number
        args.stop_level = args.coset_index; // f*(I) level
        uint64_t shifts = 0; // the PDBTs after the f*(I) level have "1 << shifts" codewords, where "shifts" = |K°(I) ∩ I|
        for (size_t level = args.coset_index+1; level < N; ++level) {
            if (__builtin_popcount(~args.coset_index & level) == 1) {
                args.sibling_levels[level] = 1;
                shifts += args.rate_profile[level];
            } else if (!args.rate_profile[level]) {
                args.stop_level = level; 
                shifts = 0;
            }
        }
        result.A_wmin += enumerate_subtree(args, args.coset_index, message) << shifts;
    }
    free(pretransform); free(message);
    free(args.rate_profile); free(args.pretransform); free(args.sibling_levels); 
    return result;
}


/*
 * Function: enumerate_minimum_weight_codewords_wrapper
 * ----------------------------------------------------
 * Wraps "enumerate_minimum_weight_codewords" so that it is easily callable from
 * a Numba accelerated Python function.
 *
 * K: Code dimension
 * N: Code length
 * generator_matrix: K×N generator matrix
 * wmin: A pointer to the minimum weight "wmin"
 * A_wmin: A pointer to the number of "wmin"-weight codewords
 */ 
void enumerate_minimum_weight_codewords_wrapper(
    size_t   K, 
    size_t   N, 
    uint8_t  generator_matrix[K][N],
    uint64_t *wmin, 
    uint64_t *A_wmin
) {
    struct result result = enumerate_minimum_weight_codewords(K, N, generator_matrix);
    *wmin = result.wmin; *A_wmin = result.A_wmin; // dealing with C structs inside of a Numba JIT function is difficult -> just use pointers 
}


/*
 * Function: main
 * --------------
 * Enumeration of the "wmin"-weight codeowords of a PAC code with polynomial
 * 0o155 and RM(3,7) rate-profile.
 */
int main() {
    // Code parameters
    int r = 3, m = 7;
    size_t N = 1 << m; 
    uint8_t polynomial[] = {1,0,1,1,0,1,1};
    size_t degree = sizeof(polynomial) - 1;
    
    // Construct the rate-profile
    uint8_t *rate_profile = malloc(sizeof(uint8_t[N])); 
    if (rate_profile == NULL) {
        exit(EXIT_FAILURE);
    }    
    reed_muller_rate_profile(r, m, rate_profile);
    
    size_t K = 0;
    for (size_t index = 0; index < N; ++index) {
        if (rate_profile[index]) {
            ++K;
        }
    }
    
    // Construct the generator matrix
    uint8_t (*generator_matrix)[N] = malloc(sizeof(uint8_t[K][N]));
    if (generator_matrix == NULL) {
        exit(EXIT_FAILURE);
    } 
    convolutional_pretransform(K, N, generator_matrix, rate_profile, degree, polynomial);
    polar_transform(K, N, generator_matrix);
    
    // Evaluate
    struct result result;
    
    int runs = 1000;
    clock_t start = clock();
    for (int i = 0; i < runs; ++i) {
        result = enumerate_minimum_weight_codewords(K, N, generator_matrix);
    }
    clock_t end = clock();
    
    // Print result
    printf("PAC RM(%d,%d):\n", r, m);
    printf("wmin: %lu, A_wmin: %lu\n", result.wmin, result.A_wmin);
    printf("Average elapsed time of %d runs: %.3e s\n", runs, (double)(end - start) / (CLOCKS_PER_SEC * runs));
    
    free(rate_profile); free(generator_matrix);
    return 0;
} 
