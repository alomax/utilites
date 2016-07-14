
//#define PI 3.14159265359
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif
#define PI M_PI


/* functions */

float* copyData(float* sample, int num_samples);
float* applySqr(float* sample, int num_samples);
float* applyAbs(float* sample, int num_samples);
float* applyTriangleSmoothing(double deltaTime, float* sample, int num_samples, double windowHalfWidthSec,
        timedomain_memory** pmem, int useMemory);
float* applyBoxcarSmoothing(double deltaTime, float* sample, int num_samples, double windowHalfWidthSec,
        timedomain_memory** pmem, int useMemory);
float* applyInstantPeriodWindowed(double deltaTime, float* sample, int num_samples, double windowWidth,
        timedomain_memory** pmem, timedomain_memory** pmem_dval, int useMemory);
float* applyIntegral(double deltaTime, float* sample, int num_samples, timedomain_memory** pmem,  int useMemory);
