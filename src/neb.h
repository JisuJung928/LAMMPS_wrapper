#ifndef __NEB_H__
#define __NEB_H__
#include "config.h"
#include "input.h"

void initial_path(Config *, Config *, Input *, int *, int);
void neb(Config *, Config *, Input *, char *, char *, int *, int);
#endif
