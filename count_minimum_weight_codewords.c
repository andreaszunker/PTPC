#include <math.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define min(a,b) ((a) < (b) ? (a) : (b))


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
            for (size_t position = 0, column = index; position <= degree && column < N; ++position, ++column) {
                pretransform[row][column] = polynomial[position];
            }
            ++row;
        }
    }
}


/*
 * Function: fast_transform2
 * -------------------------
 * Applies the polar transform to the given matrix along the second dimension.
 * 
 * rows: Number of rows of the matrix
 * columns: Number of columns of the matrix
 * matrix: Matrix on which the polar transform is to be applied
 */
void fast_transform2(size_t rows, size_t columns, uint8_t matrix[rows][columns]) {
    for (int stage = 0; stage < (int)log2(columns); ++stage) {
        int distance = 1 << stage; // Separation of the two inputs to be XORed
        for (int group = 0; group < columns; group += 2*distance) // Group iterator
            for (int butterfly = 0; butterfly < distance; ++butterfly) // Butterfly iterator
                for (int row = 0; row < rows; ++row)
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
    int current_row = 0; 
    int pivot_column = 0;

    while (current_row < rows && pivot_column < columns) {
        // Find the pivot element with a non-zero value in the current column
        int pivot_row = current_row;
        while (pivot_row < rows && matrix[pivot_row][pivot_column] == 0)
            ++pivot_row;

        if (pivot_row >= rows) {
            // No non-zero pivot element found in the current column
            ++pivot_column;
            continue;
        }
        if (current_row != pivot_row) {
            // Swap the pivot row with the current row
            for (int column = 0; column < columns; column++) {
                int temp = matrix[current_row][column];
                matrix[current_row][column] = matrix[pivot_row][column];
                matrix[pivot_row][column] = temp;
            }
        }
        // Eliminate all elements of in pivot current column except the pivot itself
        for (int row = 0; row < rows; ++row)
            if (row != current_row && matrix[row][pivot_column] != 0)
                for (int column = pivot_column; column < columns; ++column)
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
void update_message(int coset_index, int level, long message[]) {
    for (int message_index = coset_index+1; message_index < level; ++message_index) {
        if (message[message_index >> 6] & (1L << (message_index & 63)) && (~coset_index & level & message_index) == 0) {
            int update_index = (~coset_index & (level | message_index)) | (level & message_index);
            message[update_index >> 6] ^= 1L << (update_index & 63);
        }
    }
    message[level >> 6] |= 1L << (level & 63);
}


/*
 * Function: enumerate_subtree
 * ---------------------------
 * Counts the "wmin"-weight codewords contained in the given sub-tree.
 */ 
unsigned long enumerate_subtree(
    int coset_index, int start_level, int stop_level, size_t message_size, long start_message[], 
    uint8_t rate_profile[], uint8_t sibling_levels[], long pretransform[][message_size]
    ) {
    unsigned long A_wmin = 0UL;
    long *message = malloc(sizeof(long[message_size])); memcpy(message, start_message, sizeof(long[message_size]));
    update_message(coset_index, start_level, message);
    for (int level = start_level + 1; level <= stop_level; ++level) {
        if (rate_profile[level]) {
            if (sibling_levels[level]) {
                A_wmin += enumerate_subtree(coset_index, level, stop_level, message_size, message, rate_profile, sibling_levels, pretransform);
            }
        } else {
            long bit = 0L;
            for (int i = coset_index >> 6; i < ((level-1) >> 6)+1; ++i)
                bit ^= message[i] & pretransform[level][i];
            if (__builtin_popcountl(bit) % 2 != ((message[level >> 6] >> (level & 63)) & 1)) {
                if (sibling_levels[level]) {
                    update_message(coset_index, level, message);
                } else {
                    free(message); 
                    return A_wmin;
                }
            }
        }
    }
    free(message); 
    return A_wmin + 1UL;
}

/*
 * Function: count_minimum_weight_codewords
 * ----------------------------------------
 * Counts the "wmin"-weight codewords of the given pre-transformed polar code (PTPC).
 *
 * K: Code dimension
 * N: Code length
 * generator_matrix: K×N generator matrix
 * wmin_pointer: Pointer for the minimum weight "wmin"
 *
 * Returns: Number of "wmin"-weight codewords
 */ 
unsigned long count_minimum_weight_codewords(size_t K, size_t N, uint8_t generator_matrix[K][N], int* wmin_pointer) {
    // Compute the pre-transformation matrix and bring it into RREF
    uint8_t (*pretransform)[N] = malloc(sizeof(uint8_t[K][N]));
    memcpy(pretransform, generator_matrix, sizeof(uint8_t[K][N]));
    fast_transform2(K, N, pretransform);
    reduced_row_echelon_form(K, N, pretransform);    
    
    // Expand the pre-transform into a square matrix and store it with bitwise columns and in Fortran order
    uint8_t *rate_profile = calloc(N, sizeof(uint8_t)); 
    int message_size = ((N-1) >> 6)+1;
    long (*expanded_pretransform)[message_size] = calloc(N*message_size, sizeof(long)); 
                
    int pivot_column = 0;
    for (int row = 0; row < K; ++row) {
        while (pretransform[row][pivot_column] == 0)
            ++pivot_column;
        for (int column = 0; column < N; ++column)
            expanded_pretransform[column][pivot_column >> 6] |= (long)pretransform[row][column] << (pivot_column & 63);
        rate_profile[pivot_column] = 1;  
    } 
    long *message = malloc(sizeof(long[message_size])); 
    uint8_t *sibling_levels = malloc(sizeof(uint8_t[N])); 
        
    // Find minimum Hamming weight "wmin" of a coset leader
    *wmin_pointer = INT_MAX;
    for (int index = 0; index < N; ++index)
        if (rate_profile[index])
            *wmin_pointer = min(*wmin_pointer, 1 << __builtin_popcount(index));
    
    // Count the "wmin"-weight codewords in each coset
    unsigned long A_wmin = 0UL;   
    for (int coset_index = 0; coset_index < N; ++coset_index) {
        if (rate_profile[coset_index] == 0 || (1 << __builtin_popcount(coset_index)) > *wmin_pointer) continue; 
        
        // The coset is led by a "wmin"-weight row
        memset(message, 0, sizeof(long[message_size])); memset(sibling_levels, 0, sizeof(uint8_t[N])); 
        int stop_level = coset_index, coefficient = 1;
        for (int level = coset_index+1; level < N; ++level) {
            if (__builtin_popcount(~coset_index & level) == 1) {
                sibling_levels[level] = 1;
                if (rate_profile[level]) coefficient *= 2;
            } else if (rate_profile[level] == 0) {
                stop_level = level; coefficient = 1;
            }
        }
        A_wmin += coefficient * enumerate_subtree(
            coset_index, coset_index, stop_level, message_size, message, rate_profile, sibling_levels, expanded_pretransform
        );
    }
    free(rate_profile); free(expanded_pretransform); free(message); free(sibling_levels); 
    return A_wmin;
}


int main() {
    // PAC code with polynomial 0o155 and RM(3,7) rate-pofile
    int r = 3, m = 7;
    size_t N = 1 << m; 
    size_t degree = 6;
    uint8_t polynomial[] = {1,0,1,1,0,1,1};
    
    uint8_t *rate_profile = malloc(sizeof(uint8_t[N])); 
    reed_muller_rate_profile(r, m, rate_profile);
    
    size_t K = 0;
    for (size_t index = 0; index < N; ++index) 
        if (rate_profile[index]) ++K;
    
    uint8_t (*generator_matrix)[N] = malloc(sizeof(uint8_t[K][N]));
    convolutional_pretransform(K, N, generator_matrix, rate_profile, degree, polynomial);
    fast_transform2(K, N, generator_matrix);
    
    int w_min;
    unsigned long A_wmin;
    int runs = 1000;
    
    clock_t start = clock();
    for (int i = 0; i < runs; ++i)
        A_wmin = count_minimum_weight_codewords(K, N, generator_matrix, &w_min);
    clock_t end = clock();
    printf("PAC RM(%d,%d):\n", r, m);
    printf("wmin: %d, A_wmin: %lu\n", w_min, A_wmin);
    printf("Average elapsed time of %d runs: %.3e s\n", runs, (double)(end - start) / (CLOCKS_PER_SEC * runs));

    free(rate_profile); free(generator_matrix);
    return 0;
} 
