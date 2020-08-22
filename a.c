#include <math.h>
#include <stdio.h>
#include <stdint.h>

static const float π = M_PI;

int main()
{
	float F = 40000;  // sampling frequency
	float f = 440;    // base frequency
	int n = 3*F;      // total number of samples
	float x[n];       // waveform

	for (int i = 0; i < n; i++)
	{
		float t = i/F;    // time in seconds
		float ω = 2*π*f;  // angular frequency
		x[i] = sin(ω*t);  // fill waveform
		x[i] += sin(1.9*ω*t)/3;
	}

	uint8_t X[n];     // discrete waveform
	for (int i = 0; i < n; i++)
	{
		X[i] = 80*x[i] + 127;
		if (i < 20)
			fprintf(stderr, "i,X,x = %d %d %g\n", i, X[i], x[i]);
	}

	fwrite(X, 1, n, stdout);

	return 0;
}
