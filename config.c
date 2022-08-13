#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "config.h"
#include "utils.h"


/* convert into lmp basis */
void convert_basis(Config *config)
{
    int i, j;
    double A_norm, AxB_norm, vol, tmp_value, tmp_pos[3];
    double A[3], B[3], C[3];
    double Ahat[3], AxB[3], AxBhat[3], AhatxB[3], AxBhatxAhat[3];
    double frac[3][3], trans[3][3];

    /* rotate bases into triangular matrix */
    for (i = 0; i < 3; i++) {
        A[i] = config->cell[0][i];
        B[i] = config->cell[1][i];
        C[i] = config->cell[2][i];
    }

    A_norm = norm(A);
    Ahat[0] = A[0] / A_norm;
    Ahat[1] = A[1] / A_norm;
    Ahat[2] = A[2] / A_norm;

    cross(A, B, AxB);
    AxB_norm = norm(AxB);
    AxBhat[0] = AxB[0] / AxB_norm;
    AxBhat[1] = AxB[1] / AxB_norm;
    AxBhat[2] = AxB[2] / AxB_norm;

    cross(Ahat, B, AhatxB);
    cross(AxBhat, Ahat, AxBhatxAhat);

    /* column vector (a b c) */
    config->cell[0][0] = A_norm;
    config->cell[0][1] = dot(B, Ahat);
    config->cell[0][2] = dot(C, Ahat);
    config->cell[1][0] = 0.0;
    config->cell[1][1] = norm(AhatxB);
    config->cell[1][2] = dot(C, AxBhatxAhat);
    config->cell[2][0] = 0.0;
    config->cell[2][1] = 0.0;
    config->cell[2][2] = dot(C, AxBhat);

    /* edge & tilting */
    config->boxlo[0] = 0.0;
    config->boxlo[1] = 0.0;
    config->boxlo[2] = 0.0;
    config->boxhi[0] = config->cell[0][0];
    config->boxhi[1] = config->cell[1][1];
    config->boxhi[2] = config->cell[2][2];
    config->xy = config->cell[0][1];
    config->xz = config->cell[0][2];
    config->yz = config->cell[1][2];

    /* transformation matrix */
    vol = det(config->cell);
    cross(B, C, frac[0]);
    cross(C, A, frac[1]);
    cross(A, B, frac[2]);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            frac[i][j] /= vol;
        }
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            trans[i][j] = config->cell[i][0] * frac[0][j] +
                          config->cell[i][1] * frac[1][j] +
                          config->cell[i][2] * frac[2][j];
        }
    }

    /* convert position */
    for (i = 0; i < config->tot_num; i++) {
        for (j = 0; j < 3; j++) {
            tmp_pos[j] = trans[j][0] * config->pos[i * 3 + 0] +
                         trans[j][1] * config->pos[i * 3 + 1] +
                         trans[j][2] * config->pos[i * 3 + 2];
        }
        config->pos[i * 3 + 0] = tmp_pos[0];
        config->pos[i * 3 + 1] = tmp_pos[1];
        config->pos[i * 3 + 2] = tmp_pos[2];
    }

    /* transpose lattice vector */
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 3; j++) {
            tmp_value = config->cell[i][j];
            config->cell[i][j] = config->cell[j][i];
            config->cell[j][i] = tmp_value;
        }
    }
}


/* remove one atom from config */
void extract_atom(Config *config, int idx)
{
    int i, j, ntype, tot_num;
    int acc_num = 0;
    for (i = 0; i < config->ntype; ++i) {
        acc_num += config->each_num[i];
        if (acc_num > idx) {
            config->each_num[i]--;
            if (config->each_num[i] == 0) {
                for (j = i; j < config->ntype - 1; ++j) {
                    config->atom_num[i] = config->atom_num[i + 1];
                    config->each_num[i] = config->each_num[i + 1];
                }
                config->ntype--;
                ntype = config->ntype;
                config->atom_num = (int *)realloc(config->atom_num, sizeof(int) * ntype);
                config->each_num = (int *)realloc(config->each_num, sizeof(int) * ntype);
            }
            break;
        }
    }

    for (i = idx; i < config->tot_num - 1; ++i) {
        config->type[i] = config->type[i + 1];
        config->pos[i * 3 + 0] = config->pos[(i + 1) * 3 + 0];
        config->pos[i * 3 + 1] = config->pos[(i + 1) * 3 + 1];
        config->pos[i * 3 + 2] = config->pos[(i + 1) * 3 + 2];
    }
    config->tot_num--;
    tot_num = config->tot_num;
    config->id = (int *)realloc(config->id, sizeof(int) * tot_num);
    config->type = (int *)realloc(config->type, sizeof(int) * tot_num);
    config->pos = (double *)realloc(config->pos, sizeof(double) * tot_num * 3);
}


#define MAXLINE 128
int read_config(Config *config, Input *input, char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        return 1;
    }
    int i, j, k;
    double tmp_pos[3];
    char line[MAXLINE], tmp_line[MAXLINE], *ptr;

    /* system name */
    ptr = fgets(line, MAXLINE, fp);

    /* scale */
    ptr = fgets(line, MAXLINE, fp);
    double scale = atof(line);

    /* lattice vector */
    for (i = 0; i < 3; ++i) {
        ptr = fgets(line, MAXLINE, fp);
        config->cell[i][0] = atof(strtok(line, " \n")) * scale;
        config->cell[i][1] = atof(strtok(NULL, " \n")) * scale;
        config->cell[i][2] = atof(strtok(NULL, " \n")) * scale;
    }

    /* the number of type */
    ptr = fgets(line, MAXLINE, fp);
    strncpy(tmp_line, line, MAXLINE);
    config->ntype = 0;
    ptr = strtok(line, " \r\n");
    while (ptr != NULL) {
        if (strlen(ptr) > 0) {
            config->ntype++;
        }
        ptr = strtok(NULL, " \r\n");
    }

    /* atomic number type */
    ptr = strtok(tmp_line, " \n");
    config->atom_num = (int *)calloc(config->ntype, sizeof(int));
    for (i = 0; i < config->ntype; ++i) {
        config->atom_num[i] = get_atom_num(ptr);
        ptr = strtok(NULL, " \n");
    }

    /* each number of type */
    config->each_num = (int *)calloc(config->ntype, sizeof(int));
    int *start_idx = (int *)calloc(config->ntype, sizeof(int));
    config->tot_num = 0;
    ptr = fgets(line, MAXLINE, fp);
    config->each_num[0] = atoi(strtok(line, " \n"));
    config->tot_num += config->each_num[0];
    for (i = 1; i < config->ntype; ++i) {
        config->each_num[i] = atoi(strtok(NULL, " \n"));
        config->tot_num += config->each_num[i];
        start_idx[i] = config->each_num[i - 1] + start_idx[i - 1];
    }

    /* type index of atom */
    int count = 0;
    config->id = (int *)malloc(sizeof(int) * config->tot_num);
    config->type = (int *)malloc(sizeof(int) * config->tot_num);
    for (i = 0; i < config->ntype; ++i) {
        for (j = 0; j < input->nelem; ++j) {
            if (config->atom_num[i] == get_atom_num(input->symbol[j])) {
                for (k = 0; k < config->each_num[i]; ++k) {
                    config->type[count] = j + 1;
                    config->id[count] = count + 1;
                    count++;
                }
            }
        }
    }
    free(start_idx);

    /* positions and constraint */
    ptr = fgets(line, MAXLINE, fp);
    config->pos = (double *)malloc(sizeof(double) * config->tot_num * 3);
    if (strncasecmp(line, "S", 1) == 0) {
        ptr = fgets(line, MAXLINE, fp);
    }
    if (strncasecmp(line, "D", 1) == 0) {
        for (i = 0; i < config->tot_num; ++i) {
            ptr = fgets(line, MAXLINE, fp);
            tmp_pos[0] = atof(strtok(line, " \n"));
            tmp_pos[1] = atof(strtok(NULL, " \n"));
            tmp_pos[2] = atof(strtok(NULL, " \n"));
            for (j = 0; j < 3; ++j) {
                config->pos[i * 3 + j] = tmp_pos[0] * config->cell[0][j]
                                       + tmp_pos[1] * config->cell[1][j]
                                       + tmp_pos[2] * config->cell[2][j];
            }
        }
    } else {
        for (i = 0; i < config->tot_num; ++i) {
            ptr = fgets(line, MAXLINE, fp);
            config->pos[i * 3 + 0] = atof(strtok(line, " \n"));
            config->pos[i * 3 + 1] = atof(strtok(NULL, " \n"));
            config->pos[i * 3 + 2] = atof(strtok(NULL, " \n"));
        }
    }
    convert_basis(config);
    fclose(fp);
    return 0;
}


void write_config(Config *config, char *filename)
{
    int i;
    char line[MAXLINE];
    FILE *fp = fopen(filename, "a");

    /* title */
    fputs("Configuration in MINK\n", fp);

    /* scale */
    fputs("1.0\n", fp);
    
    /* lattice vector */
    for (i = 0; i < 3; ++i) {
        sprintf(line, " %19.16f", config->cell[i][0]);
        fputs(line, fp);
        sprintf(line, " %19.16f", config->cell[i][1]);
        fputs(line, fp);
        sprintf(line, " %19.16f\n", config->cell[i][2]);
        fputs(line, fp);
    }

    /* symbols */
    for (i = 0; i < config->ntype; ++i) {
        sprintf(line, " %s", get_symbol(config->atom_num[i]));
        fputs(line, fp);
    }
    fputs("\n", fp);

    /* the number of each type */
    for (i = 0; i < config->ntype; ++i) {
        sprintf(line, " %d", config->each_num[i]);
        fputs(line, fp);
    }
    fputs("\n", fp);

    double cross[3][3]; 
    cross[0][0] = config->cell[1][1] * config->cell[2][2]
                - config->cell[1][2] * config->cell[2][1];
    cross[0][1] = config->cell[1][2] * config->cell[2][0]
                - config->cell[1][0] * config->cell[2][2];
    cross[0][2] = config->cell[1][0] * config->cell[2][1]
                - config->cell[1][1] * config->cell[2][0];
    cross[1][0] = config->cell[2][1] * config->cell[0][2]
                - config->cell[2][2] * config->cell[0][1];
    cross[1][1] = config->cell[2][2] * config->cell[0][0]
                - config->cell[2][0] * config->cell[0][2];
    cross[1][2] = config->cell[2][0] * config->cell[0][1]
                - config->cell[2][1] * config->cell[0][0];
    cross[2][0] = config->cell[0][1] * config->cell[1][2]
                - config->cell[0][2] * config->cell[1][1];
    cross[2][1] = config->cell[0][2] * config->cell[1][0]
                - config->cell[0][0] * config->cell[1][2];
    cross[2][2] = config->cell[0][0] * config->cell[1][1]
                - config->cell[0][1] * config->cell[1][0];

    double vol = cross[0][0] * config->cell[0][0]
               + cross[0][1] * config->cell[0][1]
               + cross[0][2] * config->cell[0][2];

    double inv[3][3];
    inv[0][0] = cross[0][0] / vol;
    inv[0][1] = cross[1][0] / vol;
    inv[0][2] = cross[2][0] / vol;
    inv[1][0] = cross[0][1] / vol;
    inv[1][1] = cross[1][1] / vol;
    inv[1][2] = cross[2][1] / vol;
    inv[2][0] = cross[0][2] / vol;
    inv[2][1] = cross[1][2] / vol;
    inv[2][2] = cross[2][2] / vol;

    /* positions */
    fputs("Direct\n", fp);
    for (i = 0; i < config->tot_num; ++i) {
        double pos[3];
        pos[0] = config->pos[i * 3 + 0] * inv[0][0]
               + config->pos[i * 3 + 1] * inv[1][0]
               + config->pos[i * 3 + 2] * inv[2][0];
        pos[1] = config->pos[i * 3 + 0] * inv[0][1]
               + config->pos[i * 3 + 1] * inv[1][1]
               + config->pos[i * 3 + 2] * inv[2][1];
        pos[2] = config->pos[i * 3 + 0] * inv[0][2]
               + config->pos[i * 3 + 1] * inv[1][2]
               + config->pos[i * 3 + 2] * inv[2][2];
        sprintf(line, "  %19.16f", pos[0]);
        fputs(line, fp);
        sprintf(line, "  %19.16f", pos[1]);
        fputs(line, fp);
        sprintf(line, "  %19.16f\n", pos[2]);
        fputs(line, fp);
    }
    fclose(fp);
}


void copy_config(Config *config2, Config *config1)
{
    int i;

    config2->ntype = config1->ntype;
    config2->tot_num = config1->tot_num;
    config2->cell[0][0] = config1->cell[0][0];
    config2->cell[0][1] = config1->cell[0][1];
    config2->cell[0][2] = config1->cell[0][2];
    config2->cell[1][0] = config1->cell[1][0];
    config2->cell[1][1] = config1->cell[1][1];
    config2->cell[1][2] = config1->cell[1][2];
    config2->cell[2][0] = config1->cell[2][0];
    config2->cell[2][1] = config1->cell[2][1];
    config2->cell[2][2] = config1->cell[2][2];

    config2->boxlo[0] = config1->boxlo[0];
    config2->boxlo[1] = config1->boxlo[1];
    config2->boxlo[2] = config1->boxlo[2];
    config2->boxhi[0] = config1->boxhi[0];
    config2->boxhi[1] = config1->boxhi[1];
    config2->boxhi[2] = config1->boxhi[2];
    config2->xy = config1->xy;
    config2->xz = config1->xz;
    config2->yz = config1->yz;

    config2->atom_num = (int *)malloc(sizeof(int) * config1->ntype);
    config2->each_num = (int *)malloc(sizeof(int) * config1->ntype);
    config2->id = (int *)malloc(sizeof(int) * config1->tot_num);
    config2->type = (int *)malloc(sizeof(int) * config1->tot_num);
    config2->pos = (double *)malloc(sizeof(double) * config1->tot_num * 3);

    for (i = 0; i < config1->ntype; ++i) {
        config2->atom_num[i] = config1->atom_num[i];
        config2->each_num[i] = config1->each_num[i];
    }
    for (i = 0; i < config1->tot_num; ++i) {
        config2->id[i] = config1->id[i];
        config2->type[i] = config1->type[i];
        config2->pos[i * 3 + 0] = config1->pos[i * 3 + 0];
        config2->pos[i * 3 + 1] = config1->pos[i * 3 + 1];
        config2->pos[i * 3 + 2] = config1->pos[i * 3 + 2];
    }
}


void free_config(Config *config)
{
    free(config->atom_num);
    free(config->each_num);
    free(config->id);
    free(config->type);
    free(config->pos);
    free(config);
}
