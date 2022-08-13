#ifndef __INITIAL_H__
#define __INITIAL_H__

typedef struct _Input
{
    char *init_config;
    char *pair_style;
    char *pair_coeff;
    char **symbol;

    int nelem;
    int nimages;
    int init_relax;

    double max_force;
    double cutoff;
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
