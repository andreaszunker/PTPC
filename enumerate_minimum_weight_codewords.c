#include <math.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>

#define min(a,b) ((a) < (b) ? (a) : (b))


/*
 * Enumeration of Minimum Weight Codewords of Pre-Transformed Polar Codes by Tree Intersection
 * -------------------------------------------------------------------------------------------
 * Compile with: gcc -march=native -Ofast -o enumerate_minimum_weight_codewords enumerate_minimum_weight_codewords.c -lm
 */


/*
 * Function: reed_muller_rate_profile
 * ----------------------------------
 * Constructs the binary representation of the rate-profile corresponding to the RM(r,m) code.
 *
 * r: Order of the Reed-Muller code
 * m: m = log2(N) of the code length N
 * rate_profile: Buffer of length N for the binary representation of the rate-profile
 */
void reed_muller_rate_profile(int r, int m, uint8_t rate_profile[]) {
    for (size_t index; index < (1 << m); ++index) {
        rate_profile[index] = __builtin_popcountl(index) >= m-r;
    }
}


/*
 * Function: convolutional_pretransform
 * ------------------------------------
 * Constructs the pre-transformation matrix corresponding to the given rate-profile and polynomial.
 *
 * K: Code dimension
 * N: Code length
 * pretransform: Buffer for the K×N pretransfom
 * rate_profile: Binary representation of the rate-profile
 * degree: Degree of the polynomial
 * polynomial: Polynomial in binary representation with least significant coefficient first
 */
void convolutional_pretransform(size_t K, size_t N, uint8_t pretransform[K][N], uint8_t rate_profile[N], size_t degree, uint8_t polynomial[degree+1]) {
    memset(pretransform, 0, sizeof(uint8_t[K][N]));
    size_t row = 0;
    for (size_t index = 0; index < N; ++index) {
        if (rate_profile[index]) {
            for (size_t position = 0, column = index; position <= degree && column < N; ++position, ++column)
                pretransform[row][column] = polynomial[position];
            ++row;
        }
    }
}


/*
 * Function: fast_transform2
 * -------------------------
 * Applies the polar transform to the given matrix al-ong the second dimension.
 * 
 * rows: Number of rows of the matrix
 * columns: Number of columns of the matrix
 * matrix: Matrix on which the polar transform is to be applied
 */
void fast_transform2(size_t rows, size_t columns, uint8_t matrix[rows][columns]) {
    for (int stage = 0; stage < (int)log2(columns); ++stage) {
        size_t distance = 1 << stage; // Separation of the two inputs to be XORed
        for (size_t group = 0; group < columns; group += 2*distance) // Group iterator
            for (size_t butterfly = 0; butterfly < distance; ++butterfly) // Butterfly iterator
                for (size_t row = 0; row < rows; ++row)
                    matrix[row][group+butterfly] ^= matrix[row][group+butterfly+distance];
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
void reduced_row_echelon_form(size_t rows, size_t columns, uint8_t matrix[rows][columns]) {
    size_t current_row = 0, pivot_column = 0;
        
    while (current_row < rows && pivot_column < columns) {
        // Find the pivot element with a non-zero value in the current column
        size_t pivot_row = current_row;
        while (pivot_row < rows && matrix[pivot_row][pivot_column] == 0)
            ++pivot_row;

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
        // Eliminate all elements of in pivot current column except the pivot itself
        for (size_t row = 0; row < rows; ++row)
            if (row != current_row && matrix[row][pivot_column] != 0)
                for (size_t column = pivot_column; column < columns; ++column)
                    matrix[row][column] ^= matrix[current_row][column];
        
        ++current_row;
        ++pivot_column;
    }
}


/*
 * Function: update_message
 * ------------------------
 Updates the message to form a "wmin"-weight codeword of a universal polar coset.
 */ 
void update_message(size_t coset_index, size_t level, int64_t message[]) {
    // Update the message according to the "M"-set formulation
    for (size_t message_index = coset_index+1; message_index < level; ++message_index) {
        if (message[message_index >> 6] & ((int64_t)1 << (message_index & 63)) && (~coset_index & level & message_index) == 0) {
            size_t update_index = (~coset_index & (level | message_index)) | (level & message_index);
            message[update_index >> 6] ^= (int64_t)1 << (update_index & 63);
        }
    }
    message[level >> 6] |= (int64_t)1 << (level & 63); // set the "level"-th bit of the message to one. 
}


/*
 * Function: enumerate_subtree
 * ---------------------------
 * Counts the "wmin"-weight codewords contained in the given sub-tree.
 */ 
unsigned long enumerate_subtree(
    size_t coset_index, size_t start_level, size_t stop_level, size_t message_size, int64_t start_message[], 
    uint8_t rate_profile[], uint8_t sibling_levels[], int64_t pretransform[][message_size]
    ) {
    unsigned long A_wmin = 0UL;
    // Copy the message so the original message can be used later to enumerate the other sub-tree
    int64_t *message = malloc(sizeof(int64_t[message_size])); memcpy(message, start_message, sizeof(int64_t[message_size]));
    update_message(coset_index, start_level, message);
    for (size_t level = start_level + 1; level <= stop_level; ++level) {
        if (rate_profile[level]) { 
             // Sibling level of the tree correspoding to the PTPC coset
            if (sibling_levels[level]) {
                // Sibling level of the "wmin"-weight codeword tree of a universal polar coset -> sibling level of the intersection tree
                A_wmin += enumerate_subtree(coset_index, level, stop_level, message_size, message, rate_profile, sibling_levels, pretransform);
            }
        } else {
            // Einzelchild level of the intersection tree
            int64_t bit = 0L; 
            for (size_t index = coset_index >> 6; index < ((level-1) >> 6)+1; ++index)
                bit ^= message[index] & pretransform[level][index];
            if (__builtin_popcountl(bit) % 2 != ((message[level >> 6] >> (level & 63)) & 1)) {
                // The pre-transformation does not match with the current message
                if (sibling_levels[level]) {
                    // Sibling level of the "wmin"-weight codeword tree of a universal polar coset -> message can be updated
                    update_message(coset_index, level, message);
                } else {
                    // Both trees have einzelchild levels -> the message cannot be adjusted -> the message path does not form a "wmin"-weight codeword
                    free(message); 
                    return A_wmin;
                }
            }
        }
    }
    // The message path does form a "wmin"-weight codeword
    free(message); 
    return A_wmin + 1;
}


/*
 * Struct: enumeration_result
 * ---------------
 * Struct representing the result of "wmin"-weight codeword enumeration. It contains:
 * - wmin: Minimum weight "wmin"
 * - A_wmin: Number of "wmin"-weight codewords
 */
struct enumeration_result {
  unsigned int wmin;
  unsigned long A_wmin;
};


/*
 * Function: enumerate_minimum_weight_codewords
 * --------------------------------------------
 * Counts the "wmin"-weight codewords of the given pre-transformed polar code (PTPC).
 *
 * K: Code dimension
 * N: Code length
 * generator_matrix: K×N generator matrix
 *
 * Returns: A struct of type enumeration_result containing:
 * - wmin: Minimum weight "wmin"
 * - A_wmin: Number of "wmin"-weight codewords
 */ 
struct enumeration_result enumerate_minimum_weight_codewords(size_t K, size_t N, uint8_t generator_matrix[K][N]) {
    // Compute the pre-transformation matrix and bring it into RREF
    uint8_t (*pretransform)[N] = malloc(sizeof(uint8_t[K][N]));
    memcpy(pretransform, generator_matrix, sizeof(uint8_t[K][N])); // copy so that the given generator matrix remains unchanged
    fast_transform2(K, N, pretransform); // get the pre-transformation matrix
    reduced_row_echelon_form(K, N, pretransform); // bring the pre-transformation matrix into RREF 
    
    // Expand the pre-transform into a square matrix and store it with bitwise columns and in Fortran order
    uint8_t *rate_profile = calloc(N, sizeof(uint8_t)); // indicator of the sibling level of the PTPC
    size_t message_size = ((N-1) >> 6)+1; // the message is stored bitwise
    int64_t (*expanded_pretransform)[message_size] = calloc(N*message_size, sizeof(int64_t)); 
                
    size_t pivot_column = 0;
    for (size_t row = 0; row < K; ++row) {
        // find the pivot point columns of the pre-transformation matrix -> information bits
        while (pretransform[row][pivot_column] == 0)
            ++pivot_column;
        for (size_t column = 0; column < N; ++column)
            // store the rows in an expanded N×N pre-transformation with bitwise columns -> faster checking of the dynamic frozen bits
            expanded_pretransform[column][pivot_column >> 6] |= (int64_t)pretransform[row][column] << (pivot_column & 63);
        rate_profile[pivot_column] = 1;  
    } 
    int64_t *message = malloc(sizeof(int64_t[message_size])); 
    uint8_t *sibling_levels = malloc(sizeof(uint8_t[N])); // indicator of the sibling level of the "wmin"-weight codeword tree of a universal polar coset
    struct enumeration_result result = {UINT_MAX, 0UL};
    
    // Find minimum Hamming weight "wmin" of a coset leader
    for (size_t index = 0; index < N; ++index)
        if (rate_profile[index])
            result.wmin = min(result.wmin, 1U << __builtin_popcount(index));
    
    // Count the "wmin"-weight codewords in each coset
    for (size_t coset_index = 0; coset_index < N; ++coset_index) {
        if (rate_profile[coset_index] == 0 || (1U << __builtin_popcount(coset_index)) > result.wmin) continue; 
        
        // The coset is led by a "wmin"-weight row
        memset(message, 0, sizeof(int64_t[message_size])); memset(sibling_levels, 0, sizeof(uint8_t[N])); 
        
        // Find the level after which the pre-transformation cannot prevent the formation of "wmin"-weight codewords and compute their number
        size_t stop_level = coset_index; int coefficient = 1;
        for (size_t level = coset_index+1; level < N; ++level) {
            if (__builtin_popcount(~coset_index & level) == 1) {
                sibling_levels[level] = 1;
                if (rate_profile[level]) coefficient *= 2;
            } else if (rate_profile[level] == 0) {
                stop_level = level; coefficient = 1;
            }
        }
        result.A_wmin += coefficient * enumerate_subtree(
            coset_index, coset_index, stop_level, message_size, message, rate_profile, sibling_levels, expanded_pretransform
        );
    }
    free(rate_profile); free(pretransform); free(expanded_pretransform); free(message); free(sibling_levels); 
    return result;
}


/*
 * Function: main
 * --------------
 * Enumeration of the "wmin"-weight codeowords of a PAC code with polynomial 0o155 and RM(3,7) rate-pofile
 */
int main() {
    // Code parameters
    int r = 3, m = 7;
    size_t N = 1 << m; 
    uint8_t polynomial[] = {1,0,1,1,0,1,1};
    size_t degree = sizeof(polynomial) - 1;
    
    // Construct the rate-profile
    uint8_t *rate_profile = malloc(sizeof(uint8_t[N])); 
    reed_muller_rate_profile(r, m, rate_profile);
    
    size_t K = 0;
    for (size_t index = 0; index < N; ++index) 
        if (rate_profile[index]) ++K;
    
    // Construct the generator matrix
    uint8_t (*generator_matrix)[N] = malloc(sizeof(uint8_t[K][N]));
    convolutional_pretransform(K, N, generator_matrix, rate_profile, degree, polynomial);
    fast_transform2(K, N, generator_matrix);
    
    // Evaluate
    struct enumeration_result result;
    
    int runs = 1000;
    clock_t start = clock();
    for (int i = 0; i < runs; ++i)
        result = enumerate_minimum_weight_codewords(K, N, generator_matrix);
    clock_t end = clock();
    
    printf("PAC RM(%d,%d):\n", r, m);
    printf("wmin: %u, A_wmin: %lu\n", result.wmin, result.A_wmin);
    printf("Average elapsed time of %d runs: %.3e s\n", runs, (double)(end - start) / (CLOCKS_PER_SEC * runs));
    
    free(rate_profile); free(generator_matrix);
    return 0;
} 
