

/** time domain memory structure/class */

// $DOC =============================
// $DOC time domain memory for realtime processing of packets of data

typedef struct
{

    double deltaTime;

	double* input;
	double* output;

	int numInput;
	int numOutput;

}
timedomain_memory;




/* functions */

timedomain_memory* init_timedomain_memory(const int nInput, const double inputInitValue, const int nOuput, const double outputInitValue, double deltaTime);
void free_timedomain_memory(timedomain_memory** ptd_memory);
void update(double* array, int arrayLength, float* sample, int sampleLength);
void updateOutput(timedomain_memory* td_memory, float* sample, int sampleLength);
void updateInput(timedomain_memory* td_memory, float* sample, int sampleLength);
void update_d(double* array, int arrayLength, double* sample, int sampleLength);
void updateOutput_d(timedomain_memory* td_memory, double* sample, int sampleLength);
void updateInput_d(timedomain_memory* td_memory, double* sample, int sampleLength);


