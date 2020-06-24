/*************************************************************************
  scecperf.c

  test the performace, esp. the bootstrapping time of SCEC 3D velocity model
  (a FORTRAN program). 

*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/time.h>

static char *modeldir;
static int iterations;

static void runscec();
static void usage();

static const double dist1lat = 110922;
static const double dist1lon = 92382;

/*------------------------------------------------------------------------
  int main(int argc, char **argv)
  
  -------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int iter, max;
    struct timeval t_start, t_end;
    int queries;

    if (argc != 3) 
        usage();

    modeldir = argv[1];
    sscanf(argv[2], "%d", &iterations);


    chdir(modeldir);

    max = 10;
    for (iter = 0; iter < iterations; iter++) {
        int i, j, k;
        FILE *infp;
        double lon, lat, dep;
        double ms;

        infp = fopen("btestin","w+");
        if (infp == NULL){
            perror("fopen btestin");
            exit(-1);
        }

        queries = (max + 1) * (max + 1) * (max + 1);

        fprintf(infp, "%d\n",  queries);
        for (j = 0; j <= max ; j++ ) {
            lat = 33.6 + (20 * j) / dist1lat;
            
            for (i = 0; i <= max ; i++) {
                lon = -118.68 + (20 * i) / dist1lon;
                
                for (k = 0; k <= max; k++) {
                    dep = 20 * k;
                    fprintf(infp, "%f %f %f\n", lat, lon, dep);
                }
            }
        }

        fclose(infp);

        gettimeofday(&t_start, NULL);
        runscec();
        gettimeofday(&t_end, NULL);


        ms = (t_end.tv_sec - t_start.tv_sec) * 1000.0 +
            (t_end.tv_usec - t_start.tv_usec) / 1000.0;


        fprintf(stderr, 
                "%d. Points: %d Time(msec):%f Throughput(Points/sec:%f\n", 
                iter, queries, ms, queries/ms * 1000);

        max += 10;
    }

    return 0;
}


/*-----------------------------------------------------------------------
  void usage()

  -------------------------------------------------------------------------*/
void usage()
{
    fprintf(stderr, "Usage: scecperf modeldir numoftests output\n");
    exit(-1);
}

/*-----------------------------------------------------------------------
  void runscec()

  - wrapper function to call scec fortan program
  - return if OK, exit -1 on error
  -------------------------------------------------------------------------*/
void runscec()
{
    pid_t pid;

    pid = fork();
    if (pid == -1) {
        perror("fork");
        exit(-1);
    } else if (pid == 0) {
        unlink("btestout");
        execl("cvmscec", "cvmscec", (char *)0);
        perror("excel"); /* shall never get to here */
        exit(-1);
    } else {
        int status;
        waitpid(pid, &status, 0);
        if (WIFEXITED(status)) {
            return;
        }
        else {
            fprintf(stderr, "SCEC model exited abnormally\n");
            exit(-1);
        }
    }
}
