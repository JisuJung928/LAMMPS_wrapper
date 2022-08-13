#include <stdio.h>
#include <stdlib.h>
#include "calculator.h"
#define LAMMPS_LIB_MPI
#include "library.h"
#include "utils.h"


void *lmp_init(Config *config, Input *input,
               int lmpargc, char **lmpargv, MPI_Comm comm)
{
    /* create LAMMPS instance */
    int i;
    void *lmp;
    char cmd[1024];
    lmp = lammps_open(lmpargc, lmpargv, comm, NULL);
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
        sprintf(cmd, "mass %d %f", i + 1, get_mass(get_atom_num(input->symbol[i])));
        lammps_command(lmp, cmd);
    }
    return lmp;
}


double global_oneshot(Config *config, Input *input, MPI_Comm comm)
{
    char cmd[1024];
    void *lmp = NULL;
    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-log", "none", "-screen", "none"};
    //char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv, comm);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* oneshot */
    lammps_command(lmp, "run 0");
    double pe = lammps_get_thermo(lmp, "pe");
    /* delete LAMMPS instance */
    lammps_close(lmp);

    return pe;
}


double local_oneshot(Config *config, Input *input, int index, MPI_Comm comm)
{
    char cmd[1024];
    void *lmp = NULL;

    char *rlx_cmd = (char *)malloc(sizeof(char) * config->tot_num * 6);
    char *fix_cmd = (char *)malloc(sizeof(char) * config->tot_num * 6);
    /* local mask */
    int *mask = (int *)malloc(sizeof(int) * config->tot_num);
    int nummask = get_mask(config, input, rlx_cmd, fix_cmd, mask, index, comm);
    mask = (int *)realloc(mask, sizeof(int) * nummask);

    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-log", "none", "-screen", "none"};
    //char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv, comm);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* group and delete */
    lammps_command(lmp, rlx_cmd);
    if (nummask < config->tot_num) {
        lammps_command(lmp, fix_cmd);
        lammps_command(lmp, "fix 1 fix setforce 0.0 0.0 0.0");
        lammps_command(lmp, "group del subtract all rlx fix");
        lammps_command(lmp, "delete_atoms group del compress no");
    }
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* oneshot */
    lammps_command(lmp, "run 0");
    double pe = lammps_get_thermo(lmp, "pe");
    /* delete LAMMPS instance */
    lammps_close(lmp);

    free(rlx_cmd);
    free(fix_cmd);
    free(mask);

    return pe;
}


double local_oneshot_xyz(Config *config, Input *input, double *center, MPI_Comm comm)
{
    char cmd[1024];
    void *lmp = NULL;

    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-log", "none", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv, comm);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* group and delete */
    sprintf(cmd, "region double sphere %f %f %f %f",
            center[0], center[1], center[2],
            input->cutoff * 2 + input->bond_length);
    lammps_command(lmp, cmd);
    sprintf(cmd, "region single sphere %f %f %f %f",
            center[0], center[1], center[2],
            input->cutoff + input->bond_length);
    lammps_command(lmp, cmd);
    lammps_command(lmp, "group single region single");
    lammps_command(lmp, "group double region double");
    lammps_command(lmp, "group fix subtract double single");
    lammps_command(lmp, "group del subtract all double");
    lammps_command(lmp, "fix 1 fix setforce 0.0 0.0 0.0");
    lammps_command(lmp, "delete_atoms group del compress no");
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* oneshot */
    lammps_command(lmp, "run 0");
    double pe = lammps_get_thermo(lmp, "pe");
    /* delete LAMMPS instance */
    lammps_close(lmp);

    return pe;
}


double atom_relax(Config *config, Input *input, MPI_Comm comm)
{
    char cmd[1024];
    void *lmp = NULL;
    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-log", "none", "-screen", "none"};
    //char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv, comm);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /*
    lammps_command(lmp, "dump mydump all custom 1 dump.lammps id type x y z");
    lammps_command(lmp, "dump_modify mydump sort id"); 
    */
    /* minimize */
    sprintf(cmd, "minimize 0 %f 10000 100000", input->max_force);
    lammps_command(lmp, cmd);
    double pe = lammps_get_thermo(lmp, "pe");
    /* update positions */
    lammps_gather_atoms(lmp, "x", 2, 3, config->pos);
    /* delete LAMMPS instance */
    lammps_close(lmp);

    return pe;
}


double cell_relax(Config *config, Input *input, MPI_Comm comm)
{
    char cmd[1024];
    void *lmp = NULL;
    /* create LAMMPS instance */
    char *lmpargv[] = {"liblammps", "-log", "none", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(char *);
    lmp = lmp_init(config, input, lmpargc, lmpargv, comm);
    /* potential */
    sprintf(cmd, "pair_style %s", input->pair_style);
    lammps_command(lmp, cmd);
    sprintf(cmd, "pair_coeff %s", input->pair_coeff);
    lammps_command(lmp, cmd);
    /* full relax */
    lammps_command(lmp, "fix 1 all box/relax tri 0.0");
    /* balance */
    lammps_command(lmp, "balance 1.0 shift xyz 10 1.0");
    /* minimize */
    sprintf(cmd, "minimize 0 %f 10000 100000", input->max_force);
    lammps_command(lmp, cmd);
    double pe = lammps_get_thermo(lmp, "pe") / lammps_get_natoms(lmp);
    /* update positions */
    lammps_gather_atoms(lmp, "x", 2, 3, config->pos);
    /* delete LAMMPS instance */
    lammps_close(lmp);

    return pe;
}
