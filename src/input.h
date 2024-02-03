#ifndef __INITIAL_H__
#define __INITIAL_H__

typedef struct _Input
{
    char *pair_style;
    char *pair_coeff;
    char **atom_type;

    int nelem;
    int nimages;
    int oneshot;
    int atom_relax;
    int cell_relax;
    int neb;
    int dynmat;
    int nvt_md;
    int max_step;

    double max_force;
    double min_dist;
    double finite_diff;
    double timestep;
    double temperature;
} Input;

int input_int(int *, char *, char *);
int input_double(double *, char *, char *);
int input_double_arr(double **, char *, int, char *);
int input_longlong(long long *, char *, char *);
int input_char(char **, char *, char *);
int input_char_arr(char ***, char *, int, char *);
int read_input(Input *, char *);
void free_input(Input *);
#endif
