// wavecanvas: a disturbingly simple synthesizer
//
// Due to lack of expertise with common synthesizer language, I use a
// painting-based language.
//
// 	sound file   <->  canvas
// 	notes        <->  colors
// 	instruments  <->  brushes
//
// Thus, you can paint (play) the same color (note) using different
// brushes (instruments) and the texture (timbre) is different.
// There is no notion of "time", the whole sound file is seen at once and the
// note can be played at any time, not necessarily in order.


#include <assert.h>
#include <math.h>     // sin, exp, fabs, fmax
#include <stdlib.h>   // malloc, free
#include <stdint.h>   // uint8_t uint16_t
#include <stdio.h>    // fwrite, stdout


static const float π = M_PI;


struct wave_canvas {  // data structure to hold the whole canvas (sound file)
	float F;      // sampling frequency in hertz
	int n;        // total number of samples
	float *x;     // array of n raw samples
};

struct wave_brush {      // data structure to hold the timbre of an instrument
	int n;           // number of partials including the fundamental
	float f[0x100];  // relative position of the partials (f[0]=1)
	float a[0x100];  // relative amplitude of the partials
	float λ;         // exponential decay rate (in hertz)
	                 // note: for more realistic timbres, each partial
	                 // should have a different decay rate, etc
};

struct wave_orchestra { // a set of instruments (brushes)
	int n;          // number of instruments
	struct wave_brush t[0x100];
};

#define MAX_NOTES 40000
struct wave_score {          // a set of notes
	int n;               // total number of notes
	float f[MAX_NOTES];  // pitch (fundamental frequency)
	float t[MAX_NOTES];  // start of attack (time)
	float l[MAX_NOTES];  // length
	int   i[MAX_NOTES];  // instrument
};


// this is the sole function that paints a brush stroke on the canvas
static void wave_play(                  // "play" note f between T[0] and T[1]
		struct wave_canvas *w,  // sound canvas
		struct wave_brush *b,   // instrument
		float f,                // fundamental frequency (note)
		float T[2]              // time interval
		)
{
	assert(1 == b->f[0]);
	int n[2] = { T[0] * w->F, T[1] * w->F }; // discrete time interval
	fprintf(stderr, "note %g from %g to %g (%d harmonics)\n",
			f, T[0], T[1], b->n);
	for (int i = 0; i < n[1] - n[0]; i++)
	{
		int j = i + n[0]; // index in the w->n table
		if (j < 0 || j >= w->n) continue; // avoid outside of canvas

		float t = i / w->F;  // time since the start of the note

		float x = 0;        // accumulator for the j-th sample
		for (int k = 0; k < b->n; k++)
			if (f * b->f[k] * 2 < w->F) // avoid aliasing
				x += b->a[k] * sin(f * b->f[k] * 2*π*t);
		w->x[j] += exp(- b->λ * f * t) * x;
	}
}

// bureaucracy
static void wave_canvas_init(
		struct wave_canvas *w,  // sound canvas to initialize
		float s,                // number of seconds
		float F                 // frequency in hertz
		)
{
	w->n = F * s;
	w->x = malloc(w->n * sizeof*w->x);
	w->F = F;
	for (int i = 0; i < w->n; i++)
		w->x[i] = 0;
}

// pure sine wave
static void wave_brush_init_pure(struct wave_brush *b)
{
	b->n = 1;
	b->f[0] = 1;  // singular spectrum
	b->a[0] = 1;
	b->λ = 0;
}

// full harmonic spectrum
static void wave_brush_init_full(struct wave_brush *b)
{
	b->n = 0x100;
	for (int i = 0; i < b->n; i++)
	{
		b->f[i] = i + 1; // harmonic spectrum
		b->a[i] = 1;     // no decay
	}
	b->λ = 0;
}

// triangular, soft-sounding wave
static void wave_brush_init_triangular(struct wave_brush *b)
{
	b->n = 0x100;
	for (int i = 0; i < b->n; i++)
	{
		b->f[i] = i + 1;      // harmonic spectrum
		b->a[i] = 1/b->f[i];  // inverse frequency decay
	}
	b->λ = 0;
}

static void wave_quantized_stdout(struct wave_canvas *w)
{
	float M = 0;
	for (int i = 0; i < w->n; i++)
		M = fmax(M, fabs(w->x[i]));
	int16_t *x = malloc(w->n * sizeof*x);
	for (int i = 0; i < w->n; i++)
		x[i] = 10000*(w->x[i]/M);
	fwrite(x, sizeof*x, w->n, stdout);
	free(x);
}

int main()
{
	struct wave_canvas w[1];
	wave_canvas_init(w, 7, 40000);

	struct wave_brush b[1];
	wave_brush_init_triangular(b);
	b->λ = 0.004;

	for (int i = 0; i <= 12; i++)
		wave_play(w, b, 440*pow(2,i/12.0), (float[]){i/2.0, (i+1)/2.0});
	wave_play(w, b, 1*440, (float[]){0, 1});
	wave_play(w, b, 2*440, (float[]){1, 2});
	wave_play(w, b, 3*440, (float[]){2, 3});
	wave_play(w, b, 4*440, (float[]){3, 4});

	wave_quantized_stdout(w);
	return 0;
}
