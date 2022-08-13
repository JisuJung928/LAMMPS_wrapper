#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "calculator.h"
#define LAMMPS_LIB_MPI
#include "library.h"
#include "neb.h"
#include "utils.h"


void initial_path(Config *initial, Config *final, Input *input, int *mask, int nummask)
{
    int i, j, k, rank, size;
    double del[3];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* initial linear path */
    int local_size = size / input->nimages;
    int local_rank = rank % local_size;
    int img_index = rank / local_size;
    double ratio = (double)(img_index) / (input->nimages - 1);

    double *tmp_pos = (double *)malloc(sizeof(double) * nummask * 3);
    for (i = 0; i < nummask; ++i) {
        del[0] = final->pos[(mask[i] - 1) * 3 + 0]
               - initial->pos[(mask[i] - 1) * 3 + 0];
        del[1] = final->pos[(mask[i] - 1) * 3 + 1]
               - initial->pos[(mask[i] - 1) * 3 + 1];
        del[2] = final->pos[(mask[i] - 1) * 3 + 2]
               - initial->pos[(mask[i] - 1) * 3 + 2];
        get_minimum_image(del, final->boxlo, final->boxhi,
                          final->xy, final->yz, final->xz);
        tmp_pos[i * 3 + 0] = initial->pos[(mask[i] - 1) * 3 + 0]
                           + ratio * del[0];
        tmp_pos[i * 3 + 1] = initial->pos[(mask[i] - 1) * 3 + 1]
                           + ratio * del[1];
        tmp_pos[i * 3 + 2] = initial->pos[(mask[i] - 1) * 3 + 2]
                           + ratio * del[2];
    }

    /* local MPI communicator */
    int q = nummask / local_size;
    int r = nummask % local_size;
    int begin = local_rank * q + ((local_rank > r) ? r : local_rank);
    int end = begin + q;
    if (r > local_rank) {
        end++;
    }
    MPI_Comm comm_tmp;
    MPI_Comm_split(MPI_COMM_WORLD, img_index, rank, &comm_tmp);
    int count = (end - begin) * 3;
    int *counts = (int *)malloc(sizeof(int) * local_size);
    MPI_Allgather(&count, 1, MPI_INT, counts, 1, MPI_INT, comm_tmp);
    int *displ = (int *)malloc(sizeof(int) * local_size);
    displ[0] = 0;
    if (size > 1) {
        for (i = 1; i < local_size; ++i) {
            displ[i] = displ[i - 1] + counts[i - 1];
        }
    }
    /* maximum trial 1000 */
    for (i = 0; i < 1000; ++i) {
        int near = 0;
        for (j = begin; j < end; ++j) {
            for (k = 0; k < nummask; ++k) {
                if (j == k) {
                    continue;
                }
                del[0] = tmp_pos[k * 3 + 0] - tmp_pos[j * 3 + 0];
                del[1] = tmp_pos[k * 3 + 1] - tmp_pos[j * 3 + 1];
                del[2] = tmp_pos[k * 3 + 2] - tmp_pos[j * 3 + 2];
                get_minimum_image(del, final->boxlo, final->boxhi,
                                  final->xy, final->yz, final->xz);
                double dist = norm(del);
                if (dist < input->bond_length - 0.6) {
                    near = 1;
                    tmp_pos[j * 3 + 0] -= 0.1 * (input->bond_length - dist) * del[0] / dist;
                    tmp_pos[j * 3 + 1] -= 0.1 * (input->bond_length - dist) * del[1] / dist;
                    tmp_pos[j * 3 + 2] -= 0.1 * (input->bond_length - dist) * del[2] / dist;
                }
            }
        }
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       tmp_pos, counts, displ, MPI_DOUBLE, comm_tmp);
        MPI_Allreduce(MPI_IN_PLACE, &near, 1, MPI_INT, MPI_SUM, comm_tmp);
        if (near == 0) {
            break;
        }
    }
    free(counts);
    free(displ);
    MPI_Barrier(comm_tmp);
    MPI_Comm_free(&comm_tmp);

    /* write each replica */
    if (local_rank == 0) {
        char line[128];
        sprintf(line, "replica.%d", img_index);
        FILE *fp = fopen(line, "w");
        sprintf(line, "%d\n", nummask);
        fputs(line, fp);
        for (i = 0; i < nummask; ++i) {
            sprintf(line, "%d %f %f %f\n", mask[i],
                    tmp_pos[i * 3 + 0], tmp_pos[i * 3 + 1], tmp_pos[i * 3 + 2]); 
            fputs(line, fp);
        }
        fclose(fp);
    }
    free(tmp_pos);
}


void neb(Config *initial_config, Config *final_config, Input *input,
         char *rlx_cmd, char *fix_cmd, int *mask, int nummask)
{
    int i, rank, size;
    double del[3];
    char cmd[4096], partition[8];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* initial path */
    initial_path(initial_config, final_config, input, mask, nummask);

    int local_size = size / input->nimages;
    sprintf(partition, "%dx%d", input->nimages, local_size);

    /* create LAMMPS instance */
    void *lmp = NULL;
    //char *lmpargv[] = {"liblammps", "-screen", "none",
    char *lmpargv[] = {"liblammps", "-screen", "none",
                       "-partition", partition, "-in", "in.template"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(initial_config, input, lmpargc, lmpargv, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    
    /* group and delete */
    lammps_command(lmp, rlx_cmd);
    if (nummask < initial_config->tot_num) {
        lammps_command(lmp, fix_cmd);
        lammps_command(lmp, "fix 1 fix setforce 0.0 0.0 0.0");
        lammps_command(lmp, "group del subtract all rlx fix");
        lammps_command(lmp, "delete_atoms group del compress no");
    }
    lammps_command(lmp, "fix 2 rlx neb 5.0 parallel neigh");
    lammps_command(lmp, "timestep 0.01");
    lammps_command(lmp, "min_style fire");

    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* neb */
    lammps_command(lmp, "variable i equal part");
    /*
    lammps_command(lmp, "dump mydump all custom 1 dump.lammps.$i id type x y z");
    lammps_command(lmp, "dump_modify mydump sort id");
    */
    sprintf(cmd, "neb 0.0 %f 1000 0 1 each replica.$i", input->max_force);
    lammps_command(lmp, cmd);

    double *tmp_pos = (double *)malloc(sizeof(double) * nummask * 3);
    lammps_gather_atoms_subset(lmp, "x", 1, 3, nummask, mask, tmp_pos);
    for (i = 0; i < nummask; ++i) {
        final_config->pos[(mask[i] - 1) * 3 + 0] = tmp_pos[i * 3 + 0];
        final_config->pos[(mask[i] - 1) * 3 + 1] = tmp_pos[i * 3 + 1];
        final_config->pos[(mask[i] - 1) * 3 + 2] = tmp_pos[i * 3 + 2];
    }
    free(tmp_pos);

    if (rank % local_size == 0) {
        char filename[128];
        sprintf(filename, "POSCAR_%d", rank / local_size);
        write_config(final_config, filename);
    }

    /* delete LAMMPS instance */
    lammps_close(lmp);
}
