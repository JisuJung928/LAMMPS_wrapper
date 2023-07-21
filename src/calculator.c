#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "calculator.h"
#define LAMMPS_LIB_MPI
#include "library.h"
#include "utils.h"


void *lmp_init(Config *config, Input *input, int lmpargc, char **lmpargv)
{
    /* create LAMMPS instance */
    int i;
    void *lmp;
    char cmd[65536];
    lmp = lammps_open(lmpargc, lmpargv, MPI_COMM_WORLD, NULL);
    if (lmp == NULL) {
        printf("LAMMPS initialization failed");
    }
    /* basic */
    const char *cmds[] = {"units metal",
                          "neigh_modify every 1 delay 0 check yes",
                          "atom_modify map array sort 0 0.0",
                          "box tilt large"};
    lammps_commands_list(lmp, sizeof(cmds) / sizeof(const char *), cmds);
    /* box */
    sprintf(cmd, "region cell prism 0 %f 0 %f 0 %f %f %f %f units box",
            config->cell[0][0], config->cell[1][1], config->cell[2][2],
            config->cell[1][0], config->cell[2][0], config->cell[2][1]);
    lammps_command(lmp, cmd);
    sprintf(cmd, "create_box %d cell", input->nelem);
    lammps_command(lmp, cmd);
    /* atoms */
    lammps_create_atoms(lmp, config->tot_num, config->id,
                        config->type, config->pos, NULL, NULL, 0);
    for (i = 0; i < input->nelem; ++i) {
        sprintf(cmd, "mass %d %f", i + 1,
                get_mass(get_atom_num(input->atom_type[i])));
        lammps_command(lmp, cmd);
    }
    return lmp;
}


void oneshot(Config *config, Input *input)
{
    char cmd[65536];
    void *lmp = NULL;
    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* oneshot */
    lammps_command(lmp, "run 0");
    /* delete LAMMPS instance */
    lammps_close(lmp);
}


void atom_relax(Config *config, Input *input)
{
    int i, rank, size;
    char cmd[65536], tmp_cmd[65536];
    void *lmp = NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* fix */
    int fix = 0;
    for (i = 0; i < config->tot_num; ++i) {
        if (config->fix[i] > 0) {
            fix++;
            break;
        }
    }
    if (fix > 0) {
        sprintf(cmd, "group freeze id");
        for (i = 0; i < config->tot_num; ++i) {
            if (config->fix[i] > 0) {
                sprintf(tmp_cmd, " %d", i + 1);
                strcat(cmd, tmp_cmd);
            }
        }
        lammps_command(lmp, cmd);
        lammps_command(lmp, "fix 1 freeze setforce 0.0 0.0 0.0");
    }
    /* dump */
    lammps_command(lmp, "dump mydump all custom 1 dump.lammps id type x y z fx fy fz");
    lammps_command(lmp, "dump_modify mydump sort id");
    /* minimize */
    sprintf(cmd, "minimize 0 %f 10000 100000", input->max_force);
    lammps_command(lmp, cmd);
    /* update positions */
    lammps_gather_atoms(lmp, "x", 1, 3, config->pos);
    if (rank == 0) {
        write_config(config, "POSCAR_relaxed", "w");
    }
    /* delete LAMMPS instance */
    lammps_close(lmp);
}


void cell_relax(Config *config, Input *input)
{
    int i, rank, size;
    char cmd[65536], tmp_cmd[65536];
    void *lmp = NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* full relax */
    lammps_command(lmp, "fix 1 all box/relax tri 0.0");
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* fix */
    int fix = 0;
    for (i = 0; i < config->tot_num; ++i) {
        if (config->fix[i] > 0) {
            fix++;
            break;
        }
    }
    if (fix > 0) {
        sprintf(cmd, "group freeze id");
        for (i = 0; i < config->tot_num; ++i) {
            if (config->fix[i] > 0) {
                sprintf(tmp_cmd, " %d", i + 1);
                strcat(cmd, tmp_cmd); 
            }
        }
        lammps_command(lmp, cmd);
        lammps_command(lmp, "fix 1 freeze setforce 0.0 0.0 0.0");
    }
    /* dump */
    lammps_command(lmp, "dump mydump all custom 1 dump.lammps id type x y z fx fy fz");
    lammps_command(lmp, "dump_modify mydump sort id");
    /* minimize */
    sprintf(cmd, "minimize 0 %f 10000 100000", input->max_force);
    lammps_command(lmp, cmd);
    double pe = lammps_get_thermo(lmp, "pe") / lammps_get_natoms(lmp);
    /* update positions */
    lammps_gather_atoms(lmp, "x", 1, 3, config->pos);
    if (rank == 0) {
        write_config(config, "POSCAR_relaxed", "w");
    }
    /* delete LAMMPS instance */
    lammps_close(lmp);
}


static void initial_path(Config *initial, Config *final, Input *input)
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

    int tot_num = initial->tot_num;
    double *tmp_pos = (double *)malloc(sizeof(double) * tot_num * 3);
    for (i = 0; i < tot_num; ++i) {
        del[0] = final->pos[i * 3 + 0] - initial->pos[i * 3 + 0];
        del[1] = final->pos[i * 3 + 1] - initial->pos[i * 3 + 1];
        del[2] = final->pos[i * 3 + 2] - initial->pos[i * 3 + 2];
        get_minimum_image(del, final->boxlo, final->boxhi,
                          final->xy, final->yz, final->xz);
        tmp_pos[i * 3 + 0] = initial->pos[i * 3 + 0] + ratio * del[0];
        tmp_pos[i * 3 + 1] = initial->pos[i * 3 + 1] + ratio * del[1];
        tmp_pos[i * 3 + 2] = initial->pos[i * 3 + 2] + ratio * del[2];
    }

    /* local MPI communicator */
    int q = tot_num / local_size;
    int r = tot_num % local_size;
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
            for (k = 0; k < tot_num; ++k) {
                if (j == k) {
                    continue;
                }
                del[0] = tmp_pos[k * 3 + 0] - tmp_pos[j * 3 + 0];
                del[1] = tmp_pos[k * 3 + 1] - tmp_pos[j * 3 + 1];
                del[2] = tmp_pos[k * 3 + 2] - tmp_pos[j * 3 + 2];
                get_minimum_image(del, final->boxlo, final->boxhi,
                                  final->xy, final->yz, final->xz);
                double dist = norm(del);
                if (dist < input->min_dist) {
                    near = 1;
                    tmp_pos[j * 3 + 0] -= 0.1 * (3 - dist) * del[0] / dist;
                    tmp_pos[j * 3 + 1] -= 0.1 * (3 - dist) * del[1] / dist;
                    tmp_pos[j * 3 + 2] -= 0.1 * (3 - dist) * del[2] / dist;
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
        sprintf(line, "%d\n", tot_num);
        fputs(line, fp);
        for (i = 0; i < tot_num; ++i) {
            sprintf(line, "%d %f %f %f\n", i + 1,
                    tmp_pos[i * 3 + 0], tmp_pos[i * 3 + 1], tmp_pos[i * 3 + 2]);
            fputs(line, fp);
        }
        fclose(fp);
    }
    free(tmp_pos);
}


void neb(Config *initial_config, Config *final_config, Input *input)
{
    int i, rank, size;
    double del[3];
    char cmd[4096], tmp_cmd[4096], partition[8];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* initial path */
    initial_path(initial_config, final_config, input);

    int local_size = size / input->nimages;
    sprintf(partition, "%dx%d", input->nimages, local_size);

    /* create LAMMPS instance */
    void *lmp = NULL;
    //char *lmpargv[] = {"liblammps", "-screen", "none",
    char *lmpargv[] = {"liblammps", "-screen", "none",
                       "-partition", partition, "-in", "in.template"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(initial_config, input, lmpargc, lmpargv);

    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);

    /* fix */
    int fix = 0;
    for (i = 0; i < initial_config->tot_num; ++i) {
        if (initial_config->fix[i] > 0) {
            fix++;
            break;
        }
    }
    if (fix > 0) {
        sprintf(cmd, "group freeze id");
        for (i = 0; i < initial_config->tot_num; ++i) {
            if (initial_config->fix[i] > 0) {
                sprintf(tmp_cmd, " %d", i + 1);
                strcat(cmd, tmp_cmd);
            }
        }
        lammps_command(lmp, cmd);
        lammps_command(lmp, "fix 1 freeze setforce 0.0 0.0 0.0");
        lammps_command(lmp, "group neb subtract all freeze");
        lammps_command(lmp, "fix 2 neb neb 5.0 parallel neigh");
    } else {
        lammps_command(lmp, "fix 2 all neb 5.0 parallel neigh");
    }
    lammps_command(lmp, "timestep 0.002");
    lammps_command(lmp, "min_style quickmin");

    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* neb */
    lammps_command(lmp, "variable i equal part");
//    lammps_command(lmp, "dump mydump all custom 1 dump.lammps.$i id type x y z");
//    lammps_command(lmp, "dump_modify mydump sort id");
    sprintf(cmd, "neb 0.0 %f 10000 10000 1 each replica.$i", input->max_force);
//    sprintf(cmd, "neb 0.0 %f 1000 0 1 each replica.$i", input->max_force);
    lammps_command(lmp, cmd);

    /* update positions */
    lammps_gather_atoms(lmp, "x", 1, 3, final_config->pos);
    if (rank % local_size == 0) {
        char filename[128];
        sprintf(filename, "POSCAR_%d", rank / local_size);
        write_config(final_config, filename, "w");
    }

    /* delete LAMMPS instance */
    lammps_close(lmp);
}


void dynamical_matrix(Config *config, Input *input, int target_num, int *target_list)
{
    int i;
    char tmp_cmd[65536], cmd[65536];
    void *lmp = NULL;
    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* target */
    if (target_num > 0) {
        sprintf(cmd, "group target id");
        for (i = 0; i < target_num; ++i) {
            sprintf(tmp_cmd, " %d", target_list[i] + 1);
            strcat(cmd, tmp_cmd);
        }
        lammps_command(lmp, cmd);
    }
    /* dynamical_matrix */
    lammps_command(lmp, "dynamical_matrix target eskm 0.001 file dynmat.dat");
    /* delete LAMMPS instance */
    lammps_close(lmp);
}


void molecular_dynamics(Config *config, Input *input)
{
    int i, rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char tmp_cmd[65536], cmd[65536];
    void *lmp = NULL;
    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* fix */
    int fix = 0;
    for (i = 0; i < config->tot_num; ++i) {
        if (config->fix[i] > 0) {
            fix++;
            break;
        }
    }
    if (fix > 0) {
        sprintf(cmd, "group freeze id");
        for (i = 0; i < config->tot_num; ++i) {
            if (config->fix[i] > 0) {
                sprintf(tmp_cmd, " %d", i + 1);
                strcat(cmd, tmp_cmd);
            }
        }
        lammps_command(lmp, cmd);
        lammps_command(lmp, "fix 1 freeze setforce 0.0 0.0 0.0");
    }
    /* dump */
    lammps_command(lmp, "dump mydump all custom 1 dump.lammps id type x y z");
    lammps_command(lmp, "dump_modify mydump sort id");
    /* md */
    sprintf(cmd, "timestep %f\n", input->timestep / 1000.0);
    lammps_command(lmp, cmd);
    srand(time(NULL));
    int random = rand() % RAND_MAX;
    MPI_Bcast(&random, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (fix > 0) {
        lammps_command(lmp, "group mobile subtract all freeze");
        sprintf(cmd, "velocity mobile create %f %d rot yes dist gaussian",
                input->temperature, random);
        lammps_command(lmp, cmd);
        sprintf(cmd, "fix 2 mobile nvt temp %f %f $(100.0*dt)",
                input->temperature, input->temperature);
        lammps_command(lmp, cmd);
    } else {
        sprintf(cmd, "velocity all create %f %d rot yes dist gaussian",
                input->temperature, random);
        lammps_command(lmp, cmd);
        sprintf(cmd, "fix 2 all nvt temp %f %f $(100.0*dt)",
                input->temperature, input->temperature);
        lammps_command(lmp, cmd);
    }
    sprintf(cmd, "run %d\n", input->max_step);
    lammps_command(lmp, cmd);
    /* update positions */
    lammps_gather_atoms(lmp, "x", 1, 3, config->pos);
    if (rank == 0) {
        write_config(config, "POSCAR_finished", "w");
    }
    /* delete LAMMPS instance */
    lammps_close(lmp);
}
