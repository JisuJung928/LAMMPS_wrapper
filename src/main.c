#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calculator.h"
#include "config.h"
#include "input.h"
#include "target.h"
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
        if (rank == 0) {
            write_config(ini_config, "POSCAR_initial", "w");
        }
        Config *fin_config = (Config *)malloc(sizeof(Config));
        read_config(fin_config, input, argv[2]);
        atom_relax(fin_config, input);
        if (rank == 0) {
            write_config(fin_config, "POSCAR_final", "w");
        }
        /* neb */
        neb(ini_config, fin_config, input);
        if (rank == 0) {
            remove("./output/in.template");
        }
        free_config(ini_config);
        free_config(fin_config);
    } else if (input->dynmat) {
        Config *config = (Config *)malloc(sizeof(Config));
        read_config(config, input, argv[1]);
        /* read target */
        int target_num = 0;
        int list_size = 64;
        int *target_list = (int *)malloc(sizeof(int) * list_size);
        errno = read_target(config, input, &target_num, &target_list, &list_size);
        if (errno > 0) {
            printf("ERROR in TARGET FILE!\n");
            free_input(input);
            free_config(config);
            MPI_Finalize();
            return 1;
        }
        dynamical_matrix(config, input, target_num, target_list);
        free_config(config);
        free(target_list);
    } else if (input->nvt_md) {
        Config *config = (Config *)malloc(sizeof(Config));
        read_config(config, input, argv[1]);
        molecular_dynamics(config, input);
    }

    free_input(input);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
