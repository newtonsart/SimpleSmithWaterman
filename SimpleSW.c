/*
 *   Arturo Navarro Muñoz
 *   First simple implementation of the Smith-Waterman algorithm with no optimizations
 *
 *            ⎧  H(i−1,j−1)+score(A[i],B[j])     (diagonal)
 * H(i,j)=max    H(i−1,j)+gap penalty            (up)
 *            ⎨  H(i,j−1)+gap penalty            (left)
 *            ⎩  0                               (no alignment)
 *
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <getopt.h>
 #include <omp.h>
 #include <time.h>
 #include "fasta_parser.h"
 
 int GAP = -4;
 char QUERY[1024];
 char DATABASE[1024];
 static int n_procs = 0;
 
 int iBlosum62[] = {
 //  A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z    
     4, -2,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -1, -2, -4, -1, -1, -1,  1,  0, -4,  0, -3, -1, -2, -1, 
    -2,  4, -3,  4,  1, -3, -1,  0, -3, -3,  0, -4, -3,  4, -4, -2,  0, -1,  0, -1, -4, -3, -4, -1, -3,  0, 
     0, -3,  9, -3, -4, -2, -3, -3, -1, -1, -3, -1, -1, -3, -4, -3, -3, -3, -1, -1, -4, -1, -2, -1, -2, -3, 
    -2,  4, -3,  6,  2, -3, -1, -1, -3, -3, -1, -4, -3,  1, -4, -1,  0, -2,  0, -1, -4, -3, -4, -1, -3,  1, 
    -1,  1, -4,  2,  5, -3, -2,  0, -3, -3,  1, -3, -2,  0, -4, -1,  2,  0,  0, -1, -4, -2, -3, -1, -2,  4, 
    -2, -3, -2, -3, -3,  6, -3, -1,  0,  0, -3,  0,  0, -3, -4, -4, -3, -3, -2, -2, -4, -1,  1, -1,  3, -3, 
     0, -1, -3, -1, -2, -3,  6, -2, -4, -4, -2, -4, -3,  0, -4, -2, -2, -2,  0, -2, -4, -3, -2, -1, -3, -2, 
    -2,  0, -3, -1,  0, -1, -2,  8, -3, -3, -1, -3, -2,  1, -4, -2,  0,  0, -1, -2, -4, -3, -2, -1,  2,  0, 
    -1, -3, -1, -3, -3,  0, -4, -3,  4,  3, -3,  2,  1, -3, -4, -3, -3, -3, -2, -1, -4,  3, -3, -1, -1, -3, 
    -1, -3, -1, -3, -3,  0, -4, -3,  3,  3, -3,  3,  2, -3, -4, -3, -2, -2, -2, -1, -4,  2, -2, -1, -1, -3, 
    -1,  0, -3, -1,  1, -3, -2, -1, -3, -3,  5, -2, -1,  0, -4, -1,  1,  2,  0, -1, -4, -2, -3, -1, -2,  1, 
    -1, -4, -1, -4, -3,  0, -4, -3,  2,  3, -2,  4,  2, -3, -4, -3, -2, -2, -2, -1, -4,  1, -2, -1, -1, -3, 
    -1, -3, -1, -3, -2,  0, -3, -2,  1,  2, -1,  2,  5, -2, -4, -2,  0, -1, -1, -1, -4,  1, -1, -1, -1, -1, 
    -2,  4, -3,  1,  0, -3,  0,  1, -3, -3,  0, -3, -2,  6, -4, -2,  0,  0,  1,  0, -4, -3, -4, -1, -2,  0, 
    -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 
    -1, -2, -3, -1, -1, -4, -2, -2, -3, -3, -1, -3, -2, -2, -4,  7, -1, -2, -1, -1, -4, -2, -4, -1, -3, -1, 
    -1,  0, -3,  0,  2, -3, -2,  0, -3, -2,  1, -2,  0,  0, -4, -1,  5,  1,  0, -1, -4, -2, -2, -1, -1,  4, 
    -1, -1, -3, -2,  0, -3, -2,  0, -3, -2,  2, -2, -1,  0, -4, -2,  1,  5, -1, -1, -4, -3, -3, -1, -2,  0, 
     1,  0, -1,  0,  0, -2,  0, -1, -2, -2,  0, -2, -1,  1, -4, -1,  0, -1,  4,  1, -4, -2, -3, -1, -2,  0, 
     0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, -1,  0, -4, -1, -1, -1,  1,  5, -4,  0, -2, -1, -2, -1, 
    -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1, -4, -4, -4, -4, -4, 
     0, -3, -1, -3, -2, -1, -3, -3,  3,  2, -2,  1,  1, -3, -4, -2, -2, -3, -2,  0, -4,  4, -3, -1, -1, -2, 
    -3, -4, -2, -4, -3,  1, -2, -2, -3, -2, -3, -2, -1, -4, -4, -4, -2, -3, -3, -2, -4, -3, 11, -1,  2, -2, 
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4, -1, -1, -1, -1, -1, -4, -1, -1, -1, -1, -1, 
    -2, -3, -2, -3, -2,  3, -3,  2, -1, -1, -2, -1, -1, -2, -4, -3, -1, -2, -2, -2, -4, -1,  2, -1,  7, -2, 
    -1,  0, -3,  1,  4, -3, -2,  0, -3, -3,  1, -3, -1,  0, -4, -1,  4,  0,  0, -1, -4, -2, -2, -1, -2,  4, 
    };
 
int *create_matrix(int rows, int cols) {
    return (int *)calloc(rows * cols, sizeof(int));
}

void free_matrix(int *matrix) {
    free(matrix);
}

int score(int a, int b) {
    return iBlosum62[a * 26 + b];
}
 
void fill_matrix(int *H, const int *seq1, const int *seq2, int rows, int cols, int *max_i, int *max_j, int *max_score) {
    *max_score = 0;

    for (int i = 1; i < rows; i++) {
        for (int j = 1; j < cols; j++) {
            int match_score = H[(i-1)*cols + (j-1)] + score(seq1[i-1], seq2[j-1]);
            int delete = H[(i-1)*cols + j] + GAP;
            int insert = H[i*cols + (j-1)] + GAP;
            
            // Asegurar que el valor mínimo es 0
            int max_val = 0;
            max_val = (match_score > max_val) ? match_score : max_val;
            max_val = (delete > max_val) ? delete : max_val;
            max_val = (insert > max_val) ? insert : max_val;
            
            H[i*cols + j] = max_val;  // Nunca será negativo
            
            if (H[i*cols + j] > *max_score) {
               *max_score = H[i*cols + j];
               *max_i = i;
               *max_j = j;
            }
        }
    }
}
 
void print_matrix(const int *H, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%3d ", H[i * cols + j]);
        }
        printf("\n");
    }
}
 
 int *optimizeCharSeq(const char *seq, int length)
 {
     int *queryInt = malloc(length * sizeof(int));
 
 #pragma omp parallel for num_threads(n_procs) schedule(static)
     for (int i = 0; i < length; i++)
     {
         queryInt[i] = seq[i] - 'A'; // Convertir el carácter a un índice (0-25)
     }
 
     return queryInt;
 }
 
 int main(int argc, char **argv)
 {
     int opt, rows, cols;
 
     while ((opt = getopt(argc, argv, "m:s:g:q:d:")) != -1)
     {
         switch (opt)
         {
         case 'q':
             strcpy(QUERY, optarg);
             break;
         case 'd':
             strcpy(DATABASE, optarg);
             break;
         default:
             fprintf(stderr, "Usage: %s -q <query_fasta_file> -d <database_fasta_file>\n", argv[0]);
             return 1;
         }
     }
 
     if (QUERY[0] == '\0' || DATABASE[0] == '\0')
     {
         fprintf(stderr, "Usage: %s -q <query_fasta_file> -d <database_fasta_file>\n", argv[0]);
         return 1;
     }
 
     n_procs = omp_get_max_threads();
 
     FASTA_Parser *query_parser = fasta_init(QUERY);
     if (!query_parser)
         return 1;
 
    FASTA_Entry query;
    static int comp_counter = 0;
    while (fasta_next(query_parser, &query)) {
        printf(">%s\n", query.header);
        int max_i = 0, max_j = 0, max_score = 0;
        char *max_db = malloc(1024);
        rows = query.sequence_length + 1;
        int *queryInt = optimizeCharSeq(query.sequence, query.sequence_length);

        int *best_H = NULL;
        int best_rows = 0;
        int best_cols = 0;

        FASTA_Parser *db_parser = fasta_init(DATABASE);
        if (!db_parser)
            return 1;
 
        FASTA_Entry db;
        clock_t start_time = clock();

        while (fasta_next(db_parser, &db)) {
            cols = db.sequence_length + 1;

            int *H = create_matrix(rows, cols);
            int local_max_i = 0, local_max_j = 0, local_max_score = 0;
            int *dbInt = optimizeCharSeq(db.sequence, db.sequence_length);
            fill_matrix(H, queryInt, dbInt, rows, cols, &local_max_i, &local_max_j, &local_max_score);

            if (local_max_score > max_score) {
                max_score = local_max_score;
                strcpy(max_db, db.header);
                
                // Free previous best_H if exists
                if (best_H != NULL)
                    free(best_H);
                
                best_rows = rows;
                best_cols = cols;
                best_H = (int *)malloc(rows * cols * sizeof(int));
                memcpy(best_H, H, rows * cols * sizeof(int));
            }

            comp_counter++;
            printf("\rSequences compared: %d", comp_counter);
            fflush(stdout);

            free_matrix(H);
            free(dbInt);
        }
        clock_t end_time = clock();

        printf("\nMax score: %d at (%d, %d)\n", max_score, max_i, max_j);
        printf("Best match: %s\n\n", max_db);
        printf("Time: %.2f seconds\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
        print_matrix(best_H, best_rows, best_cols);

        comp_counter = 0;
        free(queryInt);
        free(max_db);
        free(best_H);
        fasta_close(db_parser);
    }

    fasta_close(query_parser);
    return 0;
}