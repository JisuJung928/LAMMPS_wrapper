#ifndef __CALCULATOR_H__
#define __CALCULATOR_H__
#include "config.h"
#include "input.h"

void *lmp_init(Config *, Input *, int, char **);
void oneshot(Config *, Input *);
void atom_relax(Config *, Input *);
void cell_relax(Config *, Input *);
void neb(Config *, Config *, Input *);
void dynamical_matrix(Config *, Input *, int, int *);
void molecular_dynamics(Config *, Input *);
#endif
