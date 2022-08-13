#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calculator.h"
#include "config.h"
#include "input.h"
#include "neb.h"
#include "utils.h"


#ifndef __M_PI__
#define M_PI 3.1415926535897932384626
#endif
int main(int argc, char *argv[])
{
    int i, errno, rank, size;
    long long kmc_step;
    double kmc_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* read input */
    Input *input = (Input *)malloc(sizeof(Input));
    errno = read_input(input, "./INPUT");
    if (errno > 0) {
        if (rank == 0) {
            printf("Check your INPUT file\n");
        }
        free_input(input);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 1;
    }
    if (rank == 0) {
        FILE *fp = fopen("./in.template", "w");
        fputs("LAMMPS input script\n", fp);
        fclose(fp);
    }
    if (size % input->nimages != 0) {
        if (rank == 0) {
            printf("Let NIMAGES be a factor of the number of processor\n");
        }
        free_input(input);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 1;
    }

    Config *ini_config = (Config *)malloc(sizeof(Config));
    read_config(ini_config, input, argv[1]);
    Config *fin_config = (Config *)malloc(sizeof(Config));
    read_config(fin_config, input, argv[2]);
    if (input->init_relax > 0) {
        atom_relax(ini_config, input, MPI_COMM_WORLD);
        atom_relax(fin_config, input, MPI_COMM_WORLD);
    }

    int *mask = (int *)malloc(sizeof(int) * ini_config->tot_num);
    char *rlx_cmd = (char *)malloc(sizeof(char) * ini_config->tot_num * 6);
    char *fix_cmd = (char *)malloc(sizeof(char) * ini_config->tot_num * 6);
    int nummask = get_mask(ini_config, input, rlx_cmd, fix_cmd, mask, atoi(argv[3]), MPI_COMM_WORLD);
    mask = (int *)realloc(mask, sizeof(int) * nummask);

    /* neb */
    neb(ini_config, fin_config, input, rlx_cmd, fix_cmd, mask, nummask);

    if (rank == 0) {
        remove("./output/in.template");
    }

    free_config(ini_config);
    free_config(fin_config);
    free_input(input);
    free(mask);
    free(rlx_cmd);
    free(fix_cmd);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
