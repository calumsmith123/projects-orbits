#define IN_FILE "./inputData/N_body_data.txt"

// maximum length of individual file lines to read in
#define MAX_FILE_LINE 250

// name mass x y x vx vy vz
#define ARG_PER_LINE 8
#define MAX_OBJ 10

#define OUT_FILE "./outputData/orbit.dat"
#define GNUPLOT_EXE "gnuplot -p"
#define GNUPLOT_SH "./src/display/SEM_script.txt"

/* CITE{http://vle.exeter.ac.uk/pluginfile.php/558109/mod_resource/content/2/ReadOrbitsFileExample.pdf}
 * Reads in data from file for lines which satisfy the correct format and saves to pointer */

static int rd_obj(struct Body *bodies) {
    FILE *rd_in = fopen(IN_FILE, "r");
    if (rd_in == NULL) {
        fprintf(stderr, "ERROR: Can't open %s\n", IN_FILE);
        exit(1);
    }
    char line[MAX_FILE_LINE], name_buf[MAX_FILE_LINE];
    int obj_num = 0, arg_num, k = 0;
    while (obj_num < MAX_OBJ && fgets(line, MAX_FILE_LINE, rd_in)) {
        ++k; // Tracks file line number
        if (line[0] != '#') {
            arg_num = sscanf(line, "%s %lg %lg %lg %lg %lg %lg %lg", name_buf, &(bodies + obj_num)->mass,
                             &(bodies + obj_num)->radius[X], &(bodies + obj_num)->radius[Y], &(bodies + obj_num)->radius[Z],
                             &(bodies + obj_num)->velocity[X], &(bodies + obj_num)->velocity[Y], &(bodies + obj_num)->velocity[Z]);
            // If correct format on line, save name to pointer
            if (arg_num == ARG_PER_LINE) {
                strncpy((bodies + obj_num)->name, name_buf, MAX_NAME);
                ++obj_num;
            } else {
                fprintf(stderr, "ERROR: line %d unreadable\n%s\n", k, line);
            }
        }
    }
    fclose(rd_in);
    return obj_num;
}

/* prints object positions to file */

static void print_file(FILE *wr_out, struct Body *bodies, int obj_num, double time) {
    fprintf(wr_out, "%e\t", time);
    for (int i = 0; i < obj_num; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(wr_out, "%e\t", (bodies + i)->radius[j]);
        }
    }
    fprintf(wr_out, "\n");
}

/* writes gnuplot script to plot data file */

static void gnuplot(struct Body *bodies, int obj_num) {
    FILE *wr_script = fopen(GNUPLOT_SH, "w");
    if (wr_script == NULL) {
        fprintf(stderr, "ERROR: Can't open %s\n", GNUPLOT_SH);
        exit(1);
    }
    int dim = 0, line;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < obj_num; ++i) {
            // detect dimensions of data
            if ((bodies + i)->radius[j] != 0) {
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
        line = N * i + 2;
        fprintf(wr_script, "'%s' using %d", OUT_FILE, line);
        for (int j = 1; j < N; ++j) {
            fprintf(wr_script, ":%d", line + j);
        }
        fprintf(wr_script, " title '%s' with linespoints pt 7, ", (bodies + i)->name);
    }
    fprintf(wr_script, "\nexit");

    fclose(wr_script);
    char command[PATH_MAX];
    snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SH);
    system(command);
}
