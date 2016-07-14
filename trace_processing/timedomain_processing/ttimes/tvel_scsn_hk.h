
//scsn_hk - depth, Vp, Vs, rho from scsn_hk.tvel in iaspei-tau ttimes
// K Hutton et. al, BSSA (2010) 
// from 0 to 800km depth only

#define NUM_TVEL_DEPTH 26

static double depth_Vp_Vs_rho[NUM_TVEL_DEPTH][4] = {
    { 0.000, 5.5000, 3.1792, 2.7200},
    { 5.000, 5.5000, 3.1792, 2.7200},
    { 5.000, 6.3000, 3.6416, 2.9200},
    { 16.000, 6.3000, 3.6416, 2.9200},
    { 16.000, 6.7000, 3.8728, 3.3198},
    { 32.000, 6.7000, 3.8728, 3.3198},
    { 32.000, 7.8000, 4.5087, 3.3198},
    { 77.50, 8.05, 4.49, 3.35},
    { 120.00, 8.05, 4.50, 3.37},
    { 165.00, 8.18, 4.51, 3.40},
    { 210.00, 8.30, 4.52, 3.43},
    { 210.00, 8.30, 4.52, 3.43},
    { 260.00, 8.48, 4.61, 3.46},
    { 310.00, 8.67, 4.70, 3.49},
    { 360.00, 8.85, 4.78, 3.52},
    { 410.00, 9.03, 4.87, 3.55},
    { 410.00, 9.36, 5.08, 3.76},
    { 460.00, 9.53, 5.19, 3.82},
    { 510.00, 9.70, 5.29, 3.88},
    { 560.00, 9.86, 5.40, 3.94},
    { 610.00, 10.03, 5.50, 4.00},
    { 660.00, 10.20, 5.61, 4.06},
    { 660.00, 10.79, 5.96, 4.37},
    { 710.00, 10.92, 6.09, 4.40},
    { 760.00, 11.06, 6.21, 4.43},
    { 809.50, 11.14, 6.24, 4.46},
};
