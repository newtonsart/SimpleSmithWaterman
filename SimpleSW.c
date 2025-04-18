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
#include "fasta_parser.h"

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

int **create_matrix(int rows, int cols)
{
    int **matrix = (int **)malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++)
    {
        matrix[i] = (int *)calloc(cols, sizeof(int)); // Inicializar con ceros
    }
    return matrix;
}

void free_matrix(int **matrix, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

int score(int a, int b)
{
    return iBlosum62[a * 26 + b];
}

void fill_matrix(int **H, const int *seq1, const int *seq2, int rows, int cols, int *max_i, int *max_j, int *max_score)
{
    *max_score = 0;

    for (int i = 1; i < rows; i++)
    {
        for (int j = 1; j < cols; j++)
        {
            int match_score = H[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1]);
            int delete = H[i - 1][j] + GAP;
            int insert = H[i][j - 1] + GAP;
            H[i][j] = match_score;
            if (delete > H[i][j])
                H[i][j] = delete;
            if (insert > H[i][j])
                H[i][j] = insert;

            if (H[i][j] > *max_score)
            {
                *max_score = H[i][j];
                *max_i = i;
                *max_j = j;
            }
        }
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
    while (fasta_next(query_parser, &query))
    {
        printf(">%s\n", query.header);
        int max_i = 0, max_j = 0, max_score = 0;
        char *max_db = malloc(1024);
        rows = query.sequence_length + 1;
        int *queryInt = optimizeCharSeq(query.sequence, rows);

        FASTA_Parser *db_parser = fasta_init(DATABASE);
        if (!db_parser)
            return 1;

        FASTA_Entry db;

        while (fasta_next(db_parser, &db))
        {
            cols = db.sequence_length + 1;

            int **H = create_matrix(rows, cols);
            int local_max_i = 0, local_max_j = 0, local_max_score = 0;
            int *dbInt = optimizeCharSeq(db.sequence, cols);
            fill_matrix(H, queryInt, dbInt, rows, cols, &local_max_i, &local_max_j, &local_max_score);

            if (local_max_score > max_score)
            {
                max_score = local_max_score;
                max_i = local_max_i;
                max_j = local_max_j;
                strcpy(max_db, db.header);
            }

            // print_results(query.sequence, db.sequence, H, rows, cols);
            comp_counter++;
            printf("\rSequences compared: %d", comp_counter);
            fflush(stdout);

            free_matrix(H, rows);
            free(dbInt);
        }
        printf("\nMax score found for this sequence: %d at (%d, %d)\n", max_score, max_i, max_j);
        printf("Best match found in database:\n %s\n\n", max_db);
        comp_counter = 0;

        free(queryInt);

        fasta_close(db_parser);
    }

    fasta_close(query_parser);
    return 0;
}