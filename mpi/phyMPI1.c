#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <tgmath.h>
#include "mpi.h"
#include "ppp_pnm.h"

// Needed for MPI
int np, self;

/*
 * Die Gravitationskonstante in m^3/(kg*s^2).
 */
static const long double G = 6.674e-11;

/*
 * Datentyp zur Beschreibung eines Koerpers.
 * (Die Repraesentation der Koerper kann natuerlich bei Bedarf so
 * geaendert werden, wie es fuer die angestrebte Loesung
 * der Aufgabe erforderlich ist.)
 */
typedef struct {
    long double mass;    /* Masse in kg */
    long double x, y;    /* x- und y-Position in m */
    long double vx, vy;  /* x- und y-Geschwindigkeit in m/s */  
} body;

/* Return seconds passed since midnight on 1970-01-01 */
static double seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + ((double)tv.tv_usec)/1000000.0;
}

/*
 * Kommentare (mit "# ...") in ".dat" Dateien ueberspringen.
 */
static void skipComments(FILE *f) {
    int n;
    int dummy; /* um "unused result" Warnungen zu unterdruecken */
    do {
	n=0;
        dummy = fscanf(f, " #%n", &n);
        if (n > 0) {
            dummy += fscanf(f, "%*[^\n]"); 
            dummy += fscanf(f, "\n");
        }
    } while (n>0);
}

/*
 * Eine ".dat" Datei mit Beschreibungen der Koerper einlesen.
 * (Format siehe Uebungsblatt).
 *    f: Dateihandle, aus dem gelesen wird
 *    n: Output-Parameter fuer die Anzahl der gelesenen Koerper
 * Die Koerper werden in einem Array von body-Strukturen
 * zurueckgeliefert. Im Fehlerfall wird NULL zurueckgegeben.
 */
body* readBodies(FILE *f, int *n) {
    int i, conv;
    body *bodies;

    skipComments(f);
    if (fscanf(f, " %d", n) != 1)
        return NULL;
    bodies = (body *) malloc(sizeof(body) * *n);
    if (bodies == NULL)
	return NULL;

    for (i=0; i<*n; i++) {
	skipComments(f);
	conv = fscanf(f, " %Lf %Lf %Lf %Lf %Lf",
		      &(bodies[i].mass),
		      &(bodies[i].x), &(bodies[i].y),
		      &(bodies[i].vx), &(bodies[i].vy));
	if (conv != 5) {
	    free(bodies);
	    return NULL;
	}
    }
    return bodies;
}

/*
 * Schreibe 'n' Koerper aus dem Array 'bodies' in die
 * durch das Dateihandle 'f' bezeichnete Datei im ".dat" Format.
 */
void writeBodies(FILE *f, const body *bodies, int n) {
    int i;
    fprintf(f, "%d\n", n);
    for (i=0; i<n; i++) {
	fprintf(f, "% 10.4Lg % 10.4Lg % 10.4Lg % 10.4Lg % 10.4Lg\n",
		bodies[i].mass, bodies[i].x, bodies[i].y,
		bodies[i].vx, bodies[i].vy);
    }
}

/*
 * Berechne den Gesamtimpuls des Systems.
 *   bodies:  Array der Koerper
 *   nBodies: Anzahl der Koerper
 *   (px,py): Output-Parameter fuer den Gesamtimpuls
 */
void totalImpulse(const body *bodies, int nBodies,
                  long double *px, long double *py)
{
    long double px_=0, py_=0;
    int i;

    for (i=0; i<nBodies; i++) {
	px_ += bodies[i].mass * bodies[i].vx;
	py_ += bodies[i].mass * bodies[i].vy;
    }
    *px = px_;
    *py = py_;
}


/*
 * Parameter fuer saveImage().
 *   width, height: Breite und Hoehe des Raumausschnitts in Metern
 *         der im Bild abgespeichert wird. (0,0) liegt im Zentrum.
 *   imgWidth, imgHeight: Breite und Hoehe des Bilds in Pixel, in dem
 *         der Raumausschnitt abgespeichert wird.
 *   imgFilePrefix: Praefix des Dateinamens fuer die Bilder. An
 *         das Praefix wird 00000.pbm, 00001.pbm, 00002.pbm, etc.
 *         angehaengt.
 */
struct ImgParams {
    long double width, height;
    int imgWidth, imgHeight;
    char *imgFilePrefix;
};

/*
 * Einfache Routine zur Ausgabe der Koerper als Bild.
 * Legt ein PBM (portable bitmap) Bild mit einem weissen
 * Pixel fuer jeden Koerper an.
 *   imgNum:  Nummer des Bildes (geht in den Dateinamen ein)
 *   bodies:  Array der Koerper
 *   nBodies: Anzahl der Koerper
 *   params:  Parameter fuer das Bild
 */
void saveImage(int imgNum, const body *bodies, int nBodies,
               const struct ImgParams *params)
{
    int i, x, y;
    const int pixels = params->imgWidth * params->imgHeight;
    char name[strlen(params->imgFilePrefix)+10];
    uint8_t *img = (uint8_t *) malloc(sizeof(uint8_t) * pixels);

    if (img == NULL) {
        fprintf(stderr, "Oops: could not allocate memory for image\n");
	return;
    }

    sprintf(name, "%s%05d.pbm", params->imgFilePrefix, imgNum);
    for (i=0; i<pixels; i++)
	img[i] = 0;

    for (i=0; i<nBodies; i++) {
	x = params->imgWidth/2  + bodies[i].x*params->imgWidth/params->width;
	y = params->imgHeight/2 - bodies[i].y*params->imgHeight/params->height;

	if (x >= 0 && x < params->imgWidth && y >= 0 && y < params->imgHeight)
	    img[y*params->imgWidth + x] = 1;
    }

    if (ppp_pnm_write(name, PNM_KIND_PBM, params->imgHeight, params->imgWidth,
                      1, img) != 0) {
        fprintf(stderr, "Error writing image\n");
    }
    free(img);
}
//TODO: statt bodies zwei pointer mit zwei count angeben
void computeNextStep(body *bodies, int count, long double deltaT, long double G) {
    long double accelerationX, accelerationY;
    long double hX, hY;
    long double denom;
    long double oldVX, oldVY;
    for(int i = 0; i < count; i++) {
        accelerationX = accelerationY = 0.00e-5;
        for(int j = 0; j < count; j++) {
            if(i != j) {
                hX = bodies[j].x - bodies[i].x;
                hY = bodies[j].y - bodies[i].y;
                denom = powl(hypotl(hX, hY), 3.0);
                /*
                printf("accX: %Lg \t accY: %Lg\n", accelerationX, accelerationY);
                printf("nennerX: %Lg \t nennerY: %Lg\n", powl(fabsl(hX), 3.0), powl(fabsl(hY), 3.0));
                printf("zaehlerX: %Lg \t zaehlerY: %Lg\n", (bodies[j].mass * hX), (bodies[j].mass * hY));
                printf("gesX: %Lg\tgesY: %Lg\n", (long double)G * (((long double)bodies[j].mass * (long double)hX) /  (long double)powl((long double)fabsl((long double)hX), (long double)3.0)), (long double)G * (((long double)bodies[j].mass * (long double)hY) /  (long double)powl((long double)fabsl((long double)hY), (long double)3.0)));
                */

                accelerationX += G * ((bodies[j].mass * hX) /  denom);
                accelerationY += G * ((bodies[j].mass * hY) /  denom);

            }

            //printf("%d body: %d step\taccX: %Lg\taccY: %Lg\n",i+1, j+1, accelerationX, accelerationY);
        }

        hX = accelerationX * deltaT;
        hY = accelerationY * deltaT;
        oldVX = bodies[i].vx;
        oldVY = bodies[i].vy;
        bodies[i].vx = bodies[i].vx + hX;
        bodies[i].vy = bodies[i].vy + hY;
        bodies[i].x = bodies[i].x + oldVX * deltaT + 0.5 * hX * deltaT;
        bodies[i].y = bodies[i].y + oldVY * deltaT + 0.5 * hY * deltaT;
    }
}



/*
 * Testprogramm fuer readBodies.
 */
int main(int argc, char *argv[]) {
    FILE *f;
    body *bodies;
    // n is number of bodies
    int i, n;
    // px, py is total impulse
    long double px, py;
    // the image to sve the progress
    struct ImgParams params;

    int steps = 366;
    // TODO: rausfinden, was deltaT fuer Werte annimmt und in getopt() auch aendern
    int deltaT = 86400;
    // if program saved images for each step
    bool picture = false;
    // the file to write the bodies in after the final calculation
    FILE *o = NULL;
    char *o_path;

    long double impulse_start_px, impulse_start_py;

    double time_start, time_loaded, time_computed, time_finished;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &self);

    int option;
    bool steps_opt, deltaT_opt, picture_opt, input_opt, output_opt;
    steps_opt = deltaT_opt = picture_opt = input_opt = output_opt = false;

    while((option = getopt(argc, argv, "S:d:pf:o:")) != -1) {
        switch(option) {
            case 'S': steps_opt = true; steps = atoi(optarg); break;
            case 'd': deltaT_opt = true; deltaT = atoi(optarg); break;
            case 'p': picture_opt = true; picture = true;break;
            case 'f':
                input_opt = true;
                f = fopen(optarg, "r");
                if (f == NULL) {
                    fprintf(stderr, "Could not open file '%s'.\n", optarg);
                    MPI_Finalize();
                    return 1;
                }
                break;
            case 'o':
                output_opt = true;
                o = fopen(optarg, "a");
                if(o == NULL) {
                    fprintf(stderr, "Could not open file '%s'.\tStandard used: \toutputXXXXX.dat\n", optarg);
                    o = fopen("output", "w");
                    if(o == NULL) {
                        printf("Cannot write in directory. Console log only. Proceed anyway? (y or n)");
                        char yesno = getchar();
                        if(yesno != 'y') {
                            printf("\nEnd program.\n");
                            MPI_Finalize();
                            return 1;
                        }
                    }
                    o_path = "output";
                } else {
                    o_path = optarg;
                    printf("Output saved in: \t\t%s\n", optarg);
                }
                break;
        }
    }

    if(!input_opt) {if(self == 0) {fprintf(stderr, "Need argument input-file: a .dat file-\n");} MPI_Finalize(); return 1;}
    if(self == 0) {
        if(!output_opt) {
            printf("No output-file set.  \tStandard used: \toutputXXXXX.dat\n");
            o = fopen("output", "w");
            if(o == NULL) {
                printf("Cannot write in directory. Proceed anyway?[y/n]");
                if(getchar() != 'y') {
                    printf("\nEnd program.\n");
                    MPI_Finalize();
                    return 1;
                }
                printf("\n");
            }
            o_path = "output";
        }
        (!steps_opt) ? printf("No number of steps set.\tStandard used: \t%d\n", steps) : printf("Number of steps:      \t\t\t%d\n", steps);
        (!deltaT_opt) ? printf("No length of steps set.\tStandard used: \t%d\n", deltaT) : printf("Length of steps:    \t\t\t%d\n", deltaT);
        (!picture_opt) ? printf("Picture Mode:         \t\t\tfalse\n") : printf("Picture Mode:       \t\t\ttrue\n");
    }

    time_start = seconds();

    //TODO: readParf fuer jeden n
    bodies = readBodies(f, &n);
    if (bodies == NULL) {
	    fprintf(stderr, "Error reading .dat file\n");
	    fclose(f);
        MPI_Finalize();
	    return 1;
    }

    // Gelesenen Koerper ausgeben.
    if(self == 0) {
        if(n > 20) {
            printf("Display all %d inputs? (y or n)", n);
            if(getchar() == 'y') {
                for (i=0; i<n; i++) {
                    printf("Body %d: mass = %Lg, x = %Lg, y = %Lg, vx = %Lg, vy = %Lg\n",
                           i, bodies[i].mass, bodies[i].x, bodies[i].y,
                           bodies[i].vx, bodies[i].vy);
                }
            }
        } else {
            for (i=0; i<n; i++) {
                printf("Body %d: mass = %Lg, x = %Lg, y = %Lg, vx = %Lg, vy = %Lg\n",
                       i, bodies[i].mass, bodies[i].x, bodies[i].y,
                       bodies[i].vx, bodies[i].vy);
            }
        }
    }

    //TODO: Impuls reduzieren
    // Calculate and print the impulse.
    totalImpulse(bodies, n, &px, &py);
    impulse_start_px = px;
    impulse_start_py = py;
    printf("Total impulse: px=%20.20Lg, py=%20.20Lg\n", px, py);

    time_loaded = seconds();
    if(self == 0) {printf("Computing...\n");}

    //TODO: picture-modus auf self == 0 und bodies mit parts ersetzen, methode auf parts abändern
    if(!picture) {
        for(int k = 0; k < steps; k++) {
            computeNextStep(bodies, n, deltaT, G);
        }
    } else {
        // Gelesene Koerper als PBM-Bild abspeichern.
        // Parameter passen zu "twogalaxies.dat".
        params.imgFilePrefix = "output";
        params.imgWidth = params.imgHeight = 200;
        params.width = params.height = 2.0e21;
        for(int k = 0; k < steps; k++) {
            saveImage(k, bodies, n, &params);
            computeNextStep(bodies, n, deltaT, G);
        }
        // TODO: nur wenn self == 0
        saveImage(steps, bodies, n, &params);
    }

    if(self == 0) {printf("Computing finished.\n");}

    time_computed = seconds();

    //TODO: Impuls reduzieren
    // Calculate and print the impulse.
    totalImpulse(bodies, n, &px, &py);
    printf("Total impulse: px=%20.20Lg, py=%20.20Lg\n", px, py);

    //TODO: nach dem Reduzieren mit self == 0 datei schreiben
    // Write output in output-file.
    if(self == 0) {
        int output_size = 70 * n * sizeof(char);
        char *output = malloc(output_size);
        snprintf(output, output_size, "# %s nach %d Schritt(en) der Laenge %d\n%d\n", o_path, steps, deltaT, n);
        for(int i=0; i<n; i++){
            int extend_size = 70 * sizeof(char);
            char *extend = malloc(extend_size);
            snprintf(extend, extend_size,
                     "% -8.5Le    % -8.5Le    % -8.5Le    % -8.5Le\n",
                     bodies[i].x, bodies[i].y, bodies[i].vx, bodies[i].vy);
            strcat(output, extend);
        }
        fputs(output, o);
    }

    time_finished = seconds();

    //TODO: Zeiten reduzieren und ausgeben lassen
    //TODO: interDev gathern oder reduzieren
    if(self == 0) {
        double interRate = time_computed-time_loaded;
        printf("Times:\n");
        printf("  Load:    %.6f s\n", time_loaded-time_start);
        printf("  Compute: %.6f s\n", interRate);
        printf("  Collect: %.6f s\n", time_finished-time_computed);
        printf("  TOTAL:   %.6f s\n", time_finished-time_start);
        interRate = n * (n - 1) * steps / interRate;
        printf("Interaction-Rate: %f\n",  interRate);
        long double interDev = impulse_start_px + impulse_start_py;
        interDev = ((interDev) - (px + py)) / fabsl(interDev);
        printf("Total interaction-deviation: %10.10f%%\n", (double)interDev);
    }

    MPI_Finalize();
    return 0;
}
