#include <stdio.h>
#include <proj_api.h>

int main(int argc, char **argv) {
	if(argc < 2) {
		printf("Usage: proj_ellips <proj string>\n");
		printf("Prints out the ellipsoid parameters for a projection as represented by proj4.\n");
		return 1;
	}
	char * p = argv[1];
	double a, e2;
	projPJ pj = pj_init_plus(p);
	pj_get_spheroid_defn(pj, &a, &e2);
	printf("%f %f\n", a, e2);
	return 0;
}