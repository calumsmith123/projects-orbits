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

#define NO_ERROR 0
#define ERROR 1

#define IN_FILE "./inputData/N_body_data.txt"
#define OUT_FILE "./outputData/orbit.dat"
#define GNUPLOT_EXE "gnuplot -p"
#define GNUPLOT_SH "./src/SEM_script.txt"

// maximum length of individual file lines to read in
#define MAX_FILE_LINE 250
// name mass x y x vx vy vz
#define ARG_PER_LINE 8
#define MAX_OBJ 10
// maximum characters in object names
#define MAX_NAME 32
// In seconds
#define TIME_STEP 86400
#define G 6.67408e-11
// precision controller solves equations to using GSL algorithms
#define GSL_ACC 1e-6

typedef enum coords {
    X, Y, Z, N_COORDS
} Coords;

struct object {
    char name[MAX_NAME];
    double mass;
    double r[N_COORDS];
    double v[N_COORDS];
    double a[N_COORDS];
} Object;

struct ptr_array {
    int centre_object;
    int k;
    int obj_num;
    double mass_sum;
    struct object ptr[MAX_OBJ];
};

struct energy {
    double ke;
    double pe;
    double tot;
};

//--------------------------------------------------//--------------------------------------------------//

/* checks character string is a positive integer */

static int number(const char *num) {
    int i = 0;
    for (; num[i] != 0; ++i) {
        if (!isdigit(num[i])) {
            i = 0;
            break;
        }
    }
    return i;
}

/* CITE{http://vle.exeter.ac.uk/pluginfile.php/558109/mod_resource/content/2/ReadOrbitsFileExample.pdf}
 * Reads in data from file for lines which satisfy the correct format and saves to pointer */

static int rd_obj(struct object *obj_ptr) {
    FILE *rd_in = fopen(IN_FILE, "r");
    if (rd_in == NULL) {
        fprintf(stderr, "ERROR: Can't open %s\n", IN_FILE);
        exit(ERROR);
    }
    char line[MAX_FILE_LINE], name_buf[MAX_FILE_LINE];
    int obj_num = 0, arg_num, k = 0;
    while (obj_num < MAX_OBJ && fgets(line, MAX_FILE_LINE, rd_in)) {
        ++k; // Tracks file line number
        if (line[0] != '#') {
            arg_num = sscanf(line, "%s %lg %lg %lg %lg %lg %lg %lg", name_buf, &(obj_ptr + obj_num)->mass,
                             &(obj_ptr + obj_num)->r[X], &(obj_ptr + obj_num)->r[Y], &(obj_ptr + obj_num)->r[Z],
                             &(obj_ptr + obj_num)->v[X], &(obj_ptr + obj_num)->v[Y], &(obj_ptr + obj_num)->v[Z]);
            // If correct format on line, save name to pointer
            if (arg_num == ARG_PER_LINE) {
                strncpy((obj_ptr + obj_num)->name, name_buf, MAX_NAME);
                ++obj_num;
            } else {
                fprintf(stderr, "ERROR: line %d unreadable\n%s\n", k, line);
            }
        }
    }
    fclose(rd_in);
    return obj_num;
}

/* select centre of system */

static int centre(struct object *obj_ptr, int obj_num) {
    printf("Centre system on object, %d defaults to centre of mass\n", obj_num);
    char input;
    while (1) {
        for (int i = 0; i < obj_num; ++i) {
            printf("%d\t%s\n", i, (obj_ptr + i)->name);
        }
        scanf("%s", &input);
        if (number(&input)) {
            int object = (int) strtol(&input, NULL, 10);
            if (object >= 0 && object <= obj_num) {
                return object;
            }
        }
        printf("input must be <= %d and >= 0\n", obj_num);
    }
}

/* Finds total system mass */

static double tot_mass(struct object *obj_ptr, int obj_num) {
    double mass_sum = 0.0;
    for (int i = 0; i < obj_num; ++i) {
        mass_sum += (obj_ptr + i)->mass;
    }
    if (mass_sum == 0.0) {
        fprintf(stderr, "ERROR: no mass found in system\n");
        exit(ERROR);
    } else {
        return mass_sum;
    }
}

/* Finds the centre of mass of the system then recalculates relative position */

static void centre_of_mass(int centre_object, struct object *obj_ptr, int obj_num, double mass_sum) {
    double centre[N_COORDS], rel_mass[N_COORDS];
    for (int j = 0; j < N_COORDS; ++j) {
        if (!centre_object) {
            rel_mass[j] = 0.0;
            for (int i = 0; i < obj_num; ++i) {
                rel_mass[j] += (obj_ptr + i)->mass * (obj_ptr + i)->r[j];
            }
            centre[j] = rel_mass[j] / mass_sum;
        } else {
            centre[j] = (obj_ptr + centre_object)->r[j];
        }
        for (int i = 0; i < obj_num; ++i) {
            // translate position vectors to normalised coordinates
            (obj_ptr + i)->r[j] -= centre[j];

        }
    }

}

/* Calculates energy of the system */

static struct energy system_energy(struct object *obj_ptr, int obj_num) {
    double m, vsq, rsq;
    struct energy energy;
    energy.ke = 0.0, energy.pe = 0.0;
    for (int i = 0; i < obj_num; ++i) {
        vsq = 0.0;
        for (int j = 0; j < N_COORDS; ++j) {
            vsq += pow((obj_ptr + i)->v[j], 2);
        }
        // Sum kinetic energy for each object
        energy.ke += 0.5 * (obj_ptr + i)->mass * vsq;
        m = 0.0;
        for (int k = 0; k < obj_num; ++k) {
            if (k > i) {
                rsq = 0.0;
                for (int j = 0; j < N_COORDS; ++j) {
                    rsq += pow((obj_ptr + i)->r[j] - (obj_ptr + k)->r[j], 2);
                }
                m += (obj_ptr + k)->mass / sqrt(rsq);
            }
        }
        // Sum gravitational potential energy for each object
        energy.pe -= G * (obj_ptr + i)->mass * m;
    }
    energy.tot = energy.ke + energy.pe;
    return energy;
}

/* Calculates acceleration of objects in the system */

static void acceleration(struct object *obj_ptr, int obj_num) {
    double rsq, rn[N_COORDS];
    for (int i = 0; i < obj_num; ++i) {
        for (int j = 0; j < N_COORDS; ++j) {
            // Reset object accelerations ready for recalculation
            (obj_ptr + i)->a[j] = 0.0;
        }
        for (int k = 0; k < obj_num; ++k) {
            if (k != i) {
                rsq = 0.0;
                for (int j = 0; j < N_COORDS; ++j) {
                    // convert x y z to r
                    rn[j] = (obj_ptr + i)->r[j] - (obj_ptr + k)->r[j];
                    rsq += pow(rn[j], 2);
                }
                if (rsq != 0.0) {
                    for (int j = 0; j < N_COORDS; ++j) {
                        (obj_ptr + i)->a[j] -= G * (obj_ptr + k)->mass * rn[j] / fabs(pow(rsq, 1.5));
                    }
                } else {
                    fprintf(stderr, "ERROR: %s has inf acceleration\n", (obj_ptr + i)->name);
                    exit(ERROR);
                }
            }
        }
    }
}

/* prints object positions to file */

static void print_file(FILE *wr_out, struct object *obj_ptr, int obj_num, double time) {
    fprintf(wr_out, "%e\t", time);
    for (int i = 0; i < obj_num; ++i) {
        for (int j = 0; j < N_COORDS; ++j) {
            fprintf(wr_out, "%e\t", (obj_ptr + i)->r[j]);
        }
    }
    fprintf(wr_out, "\n");
}

/* Calculates velocity of objects in the system */

static void velocity(struct object *obj_ptr, int obj_num) {
    for (int i = 0; i < obj_num; ++i) {
        for (int j = 0; j < N_COORDS; ++j) {
            (obj_ptr + i)->v[j] += (obj_ptr + i)->a[j] * TIME_STEP / 2;
        }
    }
}

/* Calculates position of objects in the system */

static void position(struct object *obj_ptr, int obj_num) {
    for (int i = 0; i < obj_num; ++i) {
        for (int j = 0; j < N_COORDS; ++j) {
            (obj_ptr + i)->r[j] += (obj_ptr + i)->v[j] * TIME_STEP;
        }
    }
}

/* updates motion of objects in the system for each time step using verlet method */

static int
verlet(int centre_object, struct object *obj_ptr, FILE *wr_out, int obj_num, double mass_sum, double run_time) {
    acceleration(obj_ptr, obj_num);

    double time = 0.0;
    for (int t = 1; t <= run_time; ++t) {
        print_file(wr_out, obj_ptr, obj_num, time);
        velocity(obj_ptr, obj_num);
        position(obj_ptr, obj_num);
        centre_of_mass(centre_object, obj_ptr, obj_num, mass_sum);
        acceleration(obj_ptr, obj_num);
        velocity(obj_ptr, obj_num);
        time = t * TIME_STEP;
    }
    print_file(wr_out, obj_ptr, obj_num, time);
    return NO_ERROR;
}

/* solves system of first order differential equations for object motion */

static int ode_solve(double t, const double y[], double dydt[], void *params) {
    (void) (t);
    struct ptr_array array = *(struct ptr_array *) params;
    centre_of_mass(array.centre_object, array.ptr, array.obj_num, array.mass_sum);
    acceleration(array.ptr, array.obj_num);

    for (int j = 0; j < N_COORDS; ++j) {
        dydt[j] = y[j + N_COORDS];
        dydt[j + N_COORDS] = array.ptr[array.k].a[j];
    }

    return GSL_SUCCESS;
}

/* updates motion of objects in the system for each time step using gnu scientific libraries */

static int gsl(int centre_object, const gsl_odeiv2_step_type *alg, struct object *obj_ptr, FILE *wr_out, int obj_num,
               double mass_sum, double run_time) {
    struct ptr_array array;
    array.centre_object = centre_object;
    array.obj_num = obj_num;
    array.mass_sum = mass_sum;
    for (int i = 0; i < obj_num; ++i) {
        // save data to structure to be passed to ode_solve
        array.ptr[i] = *(obj_ptr + i);
    }
    gsl_odeiv2_system sys = {ode_solve, NULL, 2 * N_COORDS, &array};
    gsl_odeiv2_driver *drive = gsl_odeiv2_driver_alloc_y_new(&sys, alg, GSL_ACC, GSL_ACC, 0.0);

    double time = 0.0, y[2 * N_COORDS];
    for (int t = 1; t <= run_time; ++t) {
        double t_next = t * TIME_STEP;
        print_file(wr_out, obj_ptr, obj_num, time);

        for (int i = 0; i < obj_num; ++i) {
            for (int j = 0; j < N_COORDS; ++j) {
                // copy position and velocity into a single array
                y[j] = (obj_ptr + i)->r[j];
                y[j + N_COORDS] = (obj_ptr + i)->v[j];
            }
            time = t_next - TIME_STEP;
            array.k = i;

            // call ode_solve within the time step
            int status = gsl_odeiv2_driver_apply(drive, &time, t_next, y);
            if (status != GSL_SUCCESS) {
                fprintf(stderr, "ERROR: return value %d\n", status);
                exit(ERROR);
            }
            for (int j = 0; j < N_COORDS; ++j) {
                // copy updated array back into position and velocity pointers
                (obj_ptr + i)->r[j] = y[j];
                (obj_ptr + i)->v[j] = y[j + N_COORDS];
            }
            array.ptr[i] = *(obj_ptr + i);
        }
        time = t_next;
    }
    gsl_odeiv2_driver_free(drive);
    print_file(wr_out, obj_ptr, obj_num, time);
    return NO_ERROR;
}

/* prints system energies */

static void print_energy(struct energy i_energy, struct energy f_energy) {
    printf("Initial energy\ttot %e\tke %e\tpe %e\n", i_energy.tot, i_energy.ke, i_energy.pe);
    printf("Final energy\ttot %e\tke %e\tpe %e\n", f_energy.tot, f_energy.ke, f_energy.pe);
    printf("Energy loss tot\tpercentage %e\n", 100.0 * fabs((i_energy.tot - f_energy.tot) / i_energy.tot));
}

/* writes gnuplot script to plot data file */

static void gnuplot(struct object *obj_ptr, int obj_num) {
    FILE *wr_script = fopen(GNUPLOT_SH, "w");
    if (wr_script == NULL) {
        fprintf(stderr, "ERROR: Can't open %s\n", GNUPLOT_SH);
        exit(ERROR);
    }
    int dim = 0, line;
    for (int j = 0; j < N_COORDS; ++j) {
        for (int i = 0; i < obj_num; ++i) {
            // detect dimensions of data
            if ((obj_ptr + i)->r[j] != 0) {
                ++dim;
                break;
            }
        }
    }
    if (dim == 2 || dim == 3) {
        fprintf(wr_script, "# %s\nset title 'Orbits'\n", GNUPLOT_SH);
        //  -----   CODE USED TO VIEW MOON ORBIT WHEN CENTERED ON EARTH -----
        // fprintf(wr_script, "set xrange[-5e8:5e8]\nset yrange[-5e8:5e8]\n");
        if (dim == 2) {
            // define range in empty dimension as default is zero
            fprintf(wr_script, "set zrange[-1.0:1.0]\n");
        }
        fprintf(wr_script, "splot ");
    } else {
        fprintf(stderr, "ERROR: data must be 2D or 3D\n");
    }

    for (int i = 0; i < obj_num; ++i) {
        // lines in data file corresponding to each object
        line = N_COORDS * i + 2;
        fprintf(wr_script, "'%s' using %d", OUT_FILE, line);
        for (int j = 1; j < N_COORDS; ++j) {
            fprintf(wr_script, ":%d", line + j);
        }
        fprintf(wr_script, " title '%s' with linespoints pt 7, ", (obj_ptr + i)->name);
    }
    fprintf(wr_script, "\nexit");

    fclose(wr_script);
    char command[PATH_MAX];
    snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SH);
    system(command);
}

/* initialises motion of objects in the system */

static void orbits(const gsl_odeiv2_step_type *alg, struct object *obj_ptr, double run_time) {
    FILE *wr_out = fopen(OUT_FILE, "w");
    if (wr_out == NULL) {
        fprintf(stderr, "ERROR: Can't open %s\n", OUT_FILE);
        exit(ERROR);
    }
    int obj_num = rd_obj(obj_ptr), succeed;
    int centre_object = centre(obj_ptr, obj_num);
    double mass_sum = tot_mass(obj_ptr, obj_num);
    // Normalises input data around centre of mass
    centre_of_mass(centre_object, obj_ptr, obj_num, mass_sum);

    struct energy i_energy = system_energy(obj_ptr, obj_num);
    if (alg == NULL) {
        succeed = verlet(centre_object, obj_ptr, wr_out, obj_num, mass_sum, run_time);
    } else {
        succeed = gsl(centre_object, alg, obj_ptr, wr_out, obj_num, mass_sum, run_time);
    }
    fclose(wr_out);
    struct energy f_energy = system_energy(obj_ptr, obj_num);
    print_energy(i_energy, f_energy);

    if (!succeed) {
        gnuplot(obj_ptr, obj_num);
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
        struct object *obj_ptr = malloc(MAX_OBJ * sizeof(Object));
        if (obj_ptr == NULL) {
            fprintf(stderr, "ERROR: Memory not allocated\n");
            exit(ERROR);
        }
        double run_time = strtol(argv[2], NULL, 10);
        if (strcmp(argv[1], "Verlet") == 0) {
            orbits(NULL, obj_ptr, run_time);
        } else if (strcmp(argv[1], "GSL_rkf45") == 0) {
            orbits(gsl_odeiv2_step_rkf45, obj_ptr, run_time);
        } else {
            orbits(gsl_odeiv2_step_rk2, obj_ptr, run_time);
        }

        free(obj_ptr);
    } else {
        printf("Two argument expected\nProgram type: 'GSL_rk2/rkf45' or 'Verlet'\nRun time 'number of days'\n");
    }
    return NO_ERROR;
}

//--------------------------------------------------//--------------------------------------------------//
