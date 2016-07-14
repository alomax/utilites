


/* functions */

float* filter_bp_bu_co_1_5_n_4(double deltaTime, float* sample, int num_samples, timedomain_memory** pmem, int useMemory);
float* filter_hp_bu_co_0_005_n_2(double deltaTime, float* sample, int num_samples, timedomain_memory** pmem, int useMemory);
float* filter_lp_bu_co_5_n_6(double deltaTime, float* sample, int num_samples, timedomain_memory** pmem, int useMemory);
float* filter_bp_bu_co_0_333_5_n_4(double deltaTime, float* sample, int num_samples, timedomain_memory** pmem, int useMemory);
float* filter_WWSSN_SP_VEL(double deltaTime, float* sample, int num_samples, timedomain_memory** pmem, int useMemory);
float* filter_WWSSN_SP_DISP(double deltaTime, float* sample, int num_samples, timedomain_memory** pmem, int useMemory);


