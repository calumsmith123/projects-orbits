//--------------------------------------------------//--------------------------------------------------//

/* Calum Smith - 07/11/2017                                                                             *
 *
 * Reads input data for 'N' bodies in cartesian coordinates. Centres position around system centre of   *
 * mass and calculates acceleration. Iterates finding new velocity, position and acceleration for user  *
 * defined number of time steps. Initial and final total energies (kinetic and potential) printed to    *
 * verify simulation accuracy.                                                                          *
 *
 * Object positions are output to orbit.dat and plotted by SEM_script.txt to create a 3D orbital plot   *
 * of all objects.                                                                                      *
 *
 * System of equations can be solved by Verlet integration or by using existing GSL libraries           *
 *
 * To compile use: gcc main.c -Wall -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm
 */

//--------------------------------------------------//--------------------------------------------------//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#include "./definitions/body.h"

#include "./configuration/validate.c"
#include "./configuration/iofile.c"
#include "./calculation/calculation.c"

//--------------------------------------------------//--------------------------------------------------//

/* prints system energies */

static void print_energy(struct energy i_energy, struct energy f_energy) {
    printf("Initial energy\ttot %e\tke %e\tpe %e\n", i_energy.tot, i_energy.ke, i_energy.pe);
    printf("Final energy\ttot %e\tke %e\tpe %e\n", f_energy.tot, f_energy.ke, f_energy.pe);
    printf("Energy loss tot\tpercentage %e\n", 100.0 * fabs((i_energy.tot - f_energy.tot) / i_energy.tot));
}

/* initialises motion of objects in the system */

static void orbits(const gsl_odeiv2_step_type *alg, struct Body *bodies, double run_time) {
    FILE *wr_out = fopen(OUT_FILE, "w");
    if (wr_out == NULL) {
        fprintf(stderr, "ERROR: Can't open %s\n", OUT_FILE);
        exit(ERROR);
    }
    int obj_num = rd_obj(bodies), succeed;
    int centre_object = centre(bodies, obj_num);
    double mass_sum = tot_mass(bodies, obj_num);
    // Normalises input data around centre of mass
    centre_of_mass(centre_object, bodies, obj_num, mass_sum);

    struct energy i_energy = system_energy(bodies, obj_num);
    if (alg == NULL) {
        succeed = verlet(centre_object, bodies, wr_out, obj_num, mass_sum, run_time);
    } else {
        succeed = gsl(centre_object, alg, bodies, wr_out, obj_num, mass_sum, run_time);
    }
    fclose(wr_out);
    struct energy f_energy = system_energy(bodies, obj_num);
    print_energy(i_energy, f_energy);

    if (!succeed) {
        gnuplot(bodies, obj_num);
    } else {
        fprintf(stderr, "ERROR: program failed\n");
        exit(ERROR);
    }
}

/* Solves for orbits using the specified method from the command line */

int main(int argc, char *argv[]) {
    if (argc == 3 && number(argv[2]) &&
        (strcmp(argv[1], "GSL_rkf45") == 0 || strcmp(argv[1], "GSL_rk2") == 0 || strcmp(argv[1], "Verlet") == 0)) {
        // Sets memory for pointer
        struct Body *bodies = malloc(MAX_OBJ * sizeof(body));
        if (bodies == NULL) {
            fprintf(stderr, "ERROR: Memory not allocated\n");
            exit(ERROR);
        }
        double run_time = strtol(argv[2], NULL, 10);
        if (strcmp(argv[1], "Verlet") == 0) {
            orbits(NULL, bodies, run_time);
        } else if (strcmp(argv[1], "GSL_rkf45") == 0) {
            orbits(gsl_odeiv2_step_rkf45, bodies, run_time);
        } else {
            orbits(gsl_odeiv2_step_rk2, bodies, run_time);
        }

        free(bodies);
    } else {
        printf("Two argument expected\nProgram type: 'GSL_rk2/rkf45' or 'Verlet'\nRun time 'number of days'\n");
    }
    return NO_ERROR;
}

//--------------------------------------------------//--------------------------------------------------//
