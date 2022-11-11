#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calculator.h"
#include "config.h"
#include "input.h"
#include "utils.h"


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

    if (input->oneshot) {
        Config *config = (Config *)malloc(sizeof(Config));
        read_config(config, input, argv[1]);
        oneshot(config, input);
        free_config(config);
    } else if (input->atom_relax) {
        Config *config = (Config *)malloc(sizeof(Config));
        read_config(config, input, argv[1]);
        atom_relax(config, input);
        free_config(config);
    } else if (input->cell_relax) {
        Config *config = (Config *)malloc(sizeof(Config));
        read_config(config, input, argv[1]);
        cell_relax(config, input);
        free_config(config);
    } else if (input->neb) {
        if (size % input->nimages != 0) {
            if (rank == 0) {
                printf("Let NIMAGES be a factor of the number of processor\n");
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
        Config *ini_config = (Config *)malloc(sizeof(Config));
        read_config(ini_config, input, argv[1]);
        atom_relax(ini_config, input);
        Config *fin_config = (Config *)malloc(sizeof(Config));
        read_config(fin_config, input, argv[2]);
        atom_relax(fin_config, input);

        /* neb */
        neb(ini_config, fin_config, input);

        if (rank == 0) {
            remove("./output/in.template");
        }

        free_config(ini_config);
        free_config(fin_config);
    }

    free_input(input);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
