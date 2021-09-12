#define NO_ERROR 0
#define ERROR 1

// In seconds
#define TIME_STEP 86400
#define G 6.67408e-11

// precision controller solves equations to using GSL algorithms
#define GSL_ACC 1e-6

struct ptr_array {
    int centre_object;
    int k;
    int obj_num;
    double mass_sum;
    struct Body ptr[MAX_OBJ];
};

struct energy {
    double ke;
    double pe;
    double tot;
};

/* select centre of system */

static int centre(struct Body *bodies, int obj_num) {
    printf("Centre system on object, %d defaults to centre of mass\n", obj_num);
    char input;
    while (1) {
        for (int i = 0; i < obj_num; ++i) {
            printf("%d\t%s\n", i, (bodies + i)->name);
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

static double tot_mass(struct Body *bodies, int obj_num) {
    double mass_sum = 0.0;
    for (int i = 0; i < obj_num; ++i) {
        mass_sum += (bodies + i)->mass;
    }
    if (mass_sum == 0.0) {
        fprintf(stderr, "ERROR: no mass found in system\n");
        exit(1);
    } else {
        return mass_sum;
    }
}

/* Finds the centre of mass of the system then recalculates relative position */

static void centre_of_mass(int centre_object, struct Body *bodies, int obj_num, double mass_sum) {
    double centre[N], rel_mass[N];
    for (int j = 0; j < N; ++j) {
        if (!centre_object) {
            rel_mass[j] = 0.0;
            for (int i = 0; i < obj_num; ++i) {
                rel_mass[j] += (bodies + i)->mass * (bodies + i)->radius[j];
            }
            centre[j] = rel_mass[j] / mass_sum;
        } else {
            centre[j] = (bodies + centre_object)->radius[j];
        }
        for (int i = 0; i < obj_num; ++i) {
            // translate position vectors to normalised coordinates
            (bodies + i)->radius[j] -= centre[j];

        }
    }

}

/* Calculates energy of the system */

static struct energy system_energy(struct Body *bodies, int obj_num) {
    double m, vsq, rsq;
    struct energy energy;
    energy.ke = 0.0, energy.pe = 0.0;
    for (int i = 0; i < obj_num; ++i) {
        vsq = 0.0;
        for (int j = 0; j < N; ++j) {
            vsq += pow((bodies + i)->velocity[j], 2);
        }
        // Sum kinetic energy for each object
        energy.ke += 0.5 * (bodies + i)->mass * vsq;
        m = 0.0;
        for (int k = 0; k < obj_num; ++k) {
            if (k > i) {
                rsq = 0.0;
                for (int j = 0; j < N; ++j) {
                    rsq += pow((bodies + i)->radius[j] - (bodies + k)->radius[j], 2);
                }
                m += (bodies + k)->mass / sqrt(rsq);
            }
        }
        // Sum gravitational potential energy for each object
        energy.pe -= G * (bodies + i)->mass * m;
    }
    energy.tot = energy.ke + energy.pe;
    return energy;
}

/* Calculates acceleration of objects in the system */

static void acceleration(struct Body *bodies, int obj_num) {
    double rsq, rn[N];
    for (int i = 0; i < obj_num; ++i) {
        for (int j = 0; j < N; ++j) {
            // Reset object accelerations ready for recalculation
            (bodies + i)->acceleration[j] = 0.0;
        }
        for (int k = 0; k < obj_num; ++k) {
            if (k != i) {
                rsq = 0.0;
                for (int j = 0; j < N; ++j) {
                    // convert x y z to r
                    rn[j] = (bodies + i)->radius[j] - (bodies + k)->radius[j];
                    rsq += pow(rn[j], 2);
                }
                if (rsq != 0.0) {
                    for (int j = 0; j < N; ++j) {
                        (bodies + i)->acceleration[j] -= G * (bodies + k)->mass * rn[j] / fabs(pow(rsq, 1.5));
                    }
                } else {
                    fprintf(stderr, "ERROR: %s has inf acceleration\n", (bodies + i)->name);
                    exit(ERROR);
                }
            }
        }
    }
}

/* Calculates velocity of objects in the system */

static void velocity(struct Body *bodies, int obj_num) {
    for (int i = 0; i < obj_num; ++i) {
        for (int j = 0; j < N; ++j) {
            (bodies + i)->velocity[j] += (bodies + i)->acceleration[j] * TIME_STEP / 2;
        }
    }
}

/* Calculates position of objects in the system */

static void position(struct Body *bodies, int obj_num) {
    for (int i = 0; i < obj_num; ++i) {
        for (int j = 0; j < N; ++j) {
            (bodies + i)->radius[j] += (bodies + i)->velocity[j] * TIME_STEP;
        }
    }
}

/* updates motion of objects in the system for each time step using verlet method */

static int
verlet(int centre_object, struct Body *bodies, FILE *wr_out, int obj_num, double mass_sum, double run_time) {
    acceleration(bodies, obj_num);

    double time = 0.0;
    for (int t = 1; t <= run_time; ++t) {
        print_file(wr_out, bodies, obj_num, time);
        velocity(bodies, obj_num);
        position(bodies, obj_num);
        centre_of_mass(centre_object, bodies, obj_num, mass_sum);
        acceleration(bodies, obj_num);
        velocity(bodies, obj_num);
        time = t * TIME_STEP;
    }
    print_file(wr_out, bodies, obj_num, time);
    return NO_ERROR;
}

/* solves system of first order differential equations for object motion */

static int ode_solve(double t, const double y[], double dydt[], void *params) {
    (void) (t);
    struct ptr_array array = *(struct ptr_array *) params;
    centre_of_mass(array.centre_object, array.ptr, array.obj_num, array.mass_sum);
    acceleration(array.ptr, array.obj_num);

    for (int j = 0; j < N; ++j) {
        dydt[j] = y[j + N];
        dydt[j + N] = array.ptr[array.k].acceleration[j];
    }

    return GSL_SUCCESS;
}

/* updates motion of objects in the system for each time step using gnu scientific libraries */

static int gsl(int centre_object, const gsl_odeiv2_step_type *alg, struct Body *bodies, FILE *wr_out, int obj_num,
               double mass_sum, double run_time) {
    struct ptr_array array;
    array.centre_object = centre_object;
    array.obj_num = obj_num;
    array.mass_sum = mass_sum;
    for (int i = 0; i < obj_num; ++i) {
        // save data to structure to be passed to ode_solve
        array.ptr[i] = *(bodies + i);
    }
    gsl_odeiv2_system sys = {ode_solve, NULL, 2 * N, &array};
    gsl_odeiv2_driver *drive = gsl_odeiv2_driver_alloc_y_new(&sys, alg, GSL_ACC, GSL_ACC, 0.0);

    double time = 0.0, y[2 * N];
    for (int t = 1; t <= run_time; ++t) {
        double t_next = t * TIME_STEP;
        print_file(wr_out, bodies, obj_num, time);

        for (int i = 0; i < obj_num; ++i) {
            for (int j = 0; j < N; ++j) {
                // copy position and velocity into a single array
                y[j] = (bodies + i)->radius[j];
                y[j + N] = (bodies + i)->velocity[j];
            }
            time = t_next - TIME_STEP;
            array.k = i;

            // call ode_solve within the time step
            int status = gsl_odeiv2_driver_apply(drive, &time, t_next, y);
            if (status != GSL_SUCCESS) {
                fprintf(stderr, "ERROR: return value %d\n", status);
                exit(ERROR);
            }
            for (int j = 0; j < N; ++j) {
                // copy updated array back into position and velocity pointers
                (bodies + i)->radius[j] = y[j];
                (bodies + i)->velocity[j] = y[j + N];
            }
            array.ptr[i] = *(bodies + i);
        }
        time = t_next;
    }
    gsl_odeiv2_driver_free(drive);
    print_file(wr_out, bodies, obj_num, time);
    return NO_ERROR;
}
