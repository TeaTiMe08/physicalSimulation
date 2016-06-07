#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "ppp_pnm.h"


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

void computeNextStep(body *bodies, int count, long double deltaT, const long double G) {
    long double accelerationX, accelerationY;
    accelerationX = accelerationY = 0.0;
    long double hX, hY;
    for(int i = 0; i < count; i++) {
        for(int j = 0; j < count; j++) {
            if(i != j) {
                hX = bodies[j].x - bodies[i].x;
                hY = bodies[j].y - bodies[i].y;
                accelerationX += G * (bodies[j].mass * hX) /  pow(fabs(hX), 3);
                accelerationY += G * (bodies[j].mass * hY) /  pow(fabs(hY), 3);
            }
        }

        bodies[i].vx = bodies[i].vx + accelerationX * deltaT;
        bodies[i].vy = bodies[i].vy + accelerationY * deltaT;
        bodies[i].x = bodies[i].x + bodies[i].x * deltaT + 0.5 * accelerationX * deltaT * deltaT;
        bodies[i].y = bodies[i].y + bodies[i].y * deltaT + 0.5 * accelerationY * deltaT * deltaT;
    }
}



/*
 * Testprogramm fuer readBodies.
 */
int main(int argc, char *argv[]) {
    FILE *f;
    char *f_path;
    body *bodies;
    // n is number of bodies
    int i, n;
    // px, py is total impulse
    long double px, py;
    // the image to sve the progress
    struct ImgParams params;

    int steps = 10;
    // TODO: rausfinden, was deltaT fuer Werte annimmt und in getopt() auch aendern
    int deltaT = 1000;
    // if program saved images for each step
    bool picture = false;
    // the file to write the bodies in after the final calculation
    FILE *o = NULL;
    char *o_path;

    double time_start, time_loaded, time_computed, time_finished;

    /*
    if (argc != 3) {
	fprintf(stderr, "Need exactly two arguments: "
                "a .dat file and an image file\n");
	return 1;
    }
    */

    int option;
    bool steps_opt, deltaT_opt, picture_opt, file_opt, outfile_opt;
    steps_opt = deltaT_opt = picture_opt = file_opt = outfile_opt = false;
    while((option = getopt(argc, argv, "S:d:pf:o:")) != 1) {
        switch(option) {
            case 'S': steps_opt = true; steps = atoi(optarg); break;
            case 'd': deltaT_opt = true; deltaT = atoi(optarg); break;
            case 'p': picture_opt = true; picture = true; break;
            case 'f':
                file_opt = true;

                f = fopen((const char*)optarg, "r");
                if(f == NULL) {
                    fprintf(stderr, "Could not open input-file '%s'.\n", optarg);
                    return 1;
                }
                break;
            case 'o':
                outfile_opt = true;
                o = fopen(optarg, "w");
                if(o == NULL) {
                    fprintf(stderr, "Could not open output-file '%s'.\n", optarg);
                    return 1;
                }
                break;
        }
    }

    /*
     * case 'f': file_opt = true; f = fopen(optarg, "r");
                f_path = malloc(sizeof(optarg)); strcpy(f_path, optarg); break;
                // open file with mode "w" for writing
            case 'o': outfile_opt = true; o = fopen(optarg, "w");
                o_path = malloc(sizeof(optarg)); strcpy(o_path, optarg); break;
                // open file with mode "r" for readin
     */

    if(!file_opt) {printf("No input-file set.  \t\t\n"); return 1;} else {printf("Input located at: \t\t%s\n", f_path);}
    (!steps_opt) ? printf("No number of steps set.\tStandard used: \t%d\n", steps) : printf("Number of steps:      \t\t%d\n", steps);
    (!deltaT_opt) ? printf("No length of steps set.\tStandard used: \t%d\n", deltaT) : printf("Length of steps:    \t\t%d\n", deltaT);
    (!picture_opt) ? printf("Picture Mode:         \t\tfalse\n") : printf("Picture Mode:       \t\ttrue\n");
    (!outfile_opt) ? printf("No output-file set.  \tStandard used: \toutputXXXXX.dat\n") : printf("Output saved in: \t\t%s\n", o_path);

    time_start = seconds();

    bodies = readBodies(f, &n);
    if (bodies == NULL) {
	fprintf(stderr, "Error reading .dat file\n");
	fclose(f);
	return 1;
    }

    /*
     * Gelesene Koerper ausgeben.
     */
    for (i=0; i<n; i++) {
	printf("Body %d: mass = %Lg, x = %Lg, y = %Lg, vx = %Lg, vy = %Lg\n",
	       i, bodies[i].mass, bodies[i].x, bodies[i].y,
           bodies[i].vx, bodies[i].vy);
    }

    time_loaded = seconds();

    for(int k = 0; k < steps; k++) {
        computeNextStep(bodies, n, deltaT, G);
        //TODO: hier noch ausgabe und bild laden
    }


    time_computed = seconds();

    /*
     * Impuls berechnen und ausgeben.
     */
    totalImpulse(bodies, n, &px, &py);
    printf("Total impulse: px=%Lg, py=%Lg\n", px, py);

    /*
    // Write output in output-file
    int output_size = 100 * n * sizeof(char);
    char *output = malloc(output_size);
    //snprintf(output, output_size, "# %s nach %d Schritt(en) der Laenge %d", o_path, steps, deltaT);
    for(i=0; i<n; i++){
        int extend_size = 110 * sizeof(char);
        char *extend = malloc(extend_size);
        //snprintf(extend + n * sizeof(extend), extend_size, "%Lg\t%Lg\t%Lg\t%Lg\t%Lg\n", bodies[i].mass, bodies[i].x, bodies[i].y, bodies[i].vx, bodies[i].vy);
        //strcat(output, extend);
    }
    const char *out = output;
    fputs(out, f);
     */

    if(picture) {
        /*
        * Gelesene Koerper als PBM-Bild abspeichern.
        * Parameter passen zu "twogalaxies.dat".
        */
        params.imgFilePrefix = "output";
        params.imgWidth = params.imgHeight = 200;
        params.width = params.height = 2.0e21;
        //TODO: Hier statt 0 noch den Index der S-Schleife einfuegen
        saveImage(0, bodies, n, &params);
    }

    time_finished = seconds();

    printf("Times:\n");
    printf("  Load:    %.6f s\n", time_loaded-time_start);
    printf("  Compute: %.6f s\n", time_computed-time_loaded);
    printf("  Collect: %.6f s\n", time_finished-time_computed);
    printf("  TOTAL:   %.6f s\n", time_finished-time_start);

    return 0;
}
