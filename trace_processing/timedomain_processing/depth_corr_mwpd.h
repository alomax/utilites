
// PREM depth corrections for Mwpd
/*
 f = Vp ,rho corr for Mo
 f = [Vp**5/2 * rho**1/2 at depth] / [Vp**5/2 * rho**1/2 at 0-71km depth]
 corr =log10(f)/1.5 -> corr for Mwpd
 */
// from 0 to 800km depth only

#define NUM_DEPTH_CORR_MWPD_PREM_DEPTH 11

static double depth_corr_mwpd_prem[NUM_DEPTH_CORR_MWPD_PREM_DEPTH][2] = {
    {0, 0.0},
    {71, 0.00},
    {80, 0.00},
    {171, -0.01},
    {220, 0.05},
    {271, 0.06},
    {371, 0.07},
    {400, 0.12},
    {471, 0.15},
    {571, 0.18},
    {600, 0.22}
};
