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
//
// to find better "painting" names
// 	orchestra    <-> palette/brush set?
// 	score        <-> drawing instructions?
//
//   Play a score using an orchestra :: draw a painting using these brushes


#include <assert.h>
#include <math.h>     // sin, exp, fabs, fmax
#include <stdlib.h>   // malloc, free
#include <stdint.h>   // uint8_t uint16_t
#include <stdio.h>    // fwrite, stdout (plus fprintf, stderr for debugging)
#include <ctype.h>    // isalpha, tolower


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
struct wave_score {          // a set of notes (in no particular order)
	int n;               // total number of notes
	float f[MAX_NOTES];  // pitch (fundamental frequency)
	float t[MAX_NOTES];  // start of attack (time)
	float l[MAX_NOTES];  // length
	int   i[MAX_NOTES];  // instrument (index on the orchestra)
};

//struct wave_temperament { // gives the fundamental pitch of each note
//
//};

// BASIC SCALE
// C is middle C
// A is about 440
// C, D, E, F, G, A, B, C D E F G A B c d e f g a b c' d' e' f' g' a' b'
// (can add more modifiers)
// accidentals (to be notated explicitly at each note):
// _B : B flat
// ^F : F sharp

static int parse_note_name(
		int *o,        // output octave (0 is that of middle C)
		int *n,        // output note from 0=C to 6=B
		int *d,        // output accidental deviation
		float *t,        // output note length
		char *a        // input note string
		)
{
	fprintf(stderr, "parsing note string \"%s\"...", a);
	*o = *n = *d = 0;
	*t = 1;
	int i = 0;

	// prefix: accidentals
	while (a[i] && !isalpha(a[i]))
	{
		if (a[i] == '^') *d += 1;
		if (a[i] == '_') *d -= 1;
		i += 1;
	}

	// the actual note
	switch(tolower(a[i])) {
	case 'c': *n = 0; break;
	case 'd': *n = 1; break;
	case 'e': *n = 2; break;
	case 'f': *n = 3; break;
	case 'g': *n = 4; break;
	case 'a': *n = 5; break;
	case 'b': *n = 6; break;
	case 'z': *n = -1; break;
	default: fprintf(stderr, "ERROR: unrecognized note %c\n", a[i]);
	}
	if (islower(a[i])) *o += 1;
	i += 1;

	// first suffix: octave modifier
	while (a[i] && (a[i]=='\'' || a[i]==','))
	{
		if (a[i] == '\'') *o += 1;
		if (a[i] == ',' ) *o -= 1;
		i += 1;
	}

	// second suffix: note duration
	float p = 0; // number before the slash
	float q = 0; // number after the slash
	while (a[i] && isdigit(a[i]))
	{
		p = 10*p + a[i]-'0';
		i += 1;
	}
	if (a[i] == '/')
	{
		i += 1;
		while (a[i] && isdigit(a[i]))
		{
			q = 10*q + a[i]-'0';
			i += 1;
		}
	}
	if (!p) p = 1;
	if (!q) q = 1;
	*t = p/q;

	// debug
	fprintf(stderr, "\tn=%d o=%d d=%d t=%g\n", *n, *o, *d, *t);
	return i;
}

static int test_parse_note_name(void)
{
	int o, n, d;
	float t;
	parse_note_name(&o, &n, &d, &t, "a");
	parse_note_name(&o, &n, &d, &t, "A");
	parse_note_name(&o, &n, &d, &t, "^A");
	parse_note_name(&o, &n, &d, &t, "A,");
	parse_note_name(&o, &n, &d, &t, "a'");
	parse_note_name(&o, &n, &d, &t, "a'3");
	parse_note_name(&o, &n, &d, &t, "a'3/2");
	parse_note_name(&o, &n, &d, &t, "a'355/113");
	parse_note_name(&o, &n, &d, &t, "a'/4");
	parse_note_name(&o, &n, &d, &t, "z4");
}



static int chromatic_note(  // compute a chromatic note from major+accidentals
	int n,              // input note in major scale
	int d               // input accidental deviation
	)
{
	if (n < 0) return -1;
	assert(n >= 0 && n < 7);
	int m[7] = {0, 2, 4, 5, 7, 9, 11};
	return m[n] + d;
}

//static float equal_temperament_pitch_from_note_name(char *a)
//{
//	int o; // octave (0=octave that contains middle C)
//	int n; // note within the octave (from 0=c to 6=b)
//	int d; // accidental increment
//	float t; // note length
//	parse_note_name(&o, &n, &d, &t, a);
//	int x = chromatic_note(n, d) + 12*o + 9;
//	return d >= 0 ? 440 * pow(2, x/12) : 0;
//}

//static void get_pitch_and_duration(float *f, float *t, char *a)
//{
//	int o; // octave (0=octave that contains middle C)
//	int n; // note within the octave (from 0=c to 6=b)
//	int d; // accidental increment
//	//float t; // note length
//	parse_note_name(&o, &n, &d, t, a);
//	float x = chromatic_note(n, d) + 12*o - 9;
//	*f = d >= 0 ? 440 * pow(2, x/12) : 0;
//	fprintf(stderr, "gpd \"%s\" x=%g f=%g t=%g\n", a, x, *f, *t);
//}

// return a pointer to the remaning part of the string
static char *parse_pitch_and_duration(float *f, float *t, char *a)
{
	int o; // octave (0=octave that contains middle C)
	int n; // note within the octave (from 0=c to 6=b)
	int d; // accidental increment
	//float t; // note length
	int i = parse_note_name(&o, &n, &d, t, a);
	float x = chromatic_note(n, d) + 12*o - 9;
	*f = n >= 0 ? 440 * pow(2, x/12) : 0;
	//fprintf(stderr, "gpd \"%s\" x=%g f=%g t=%g\n", a, x, *f, *t);
	return a + i;
}

static int test_get_pitch_and_duration(void)
{
	float f, t;
	parse_pitch_and_duration(&f, &t, "A,,");
	parse_pitch_and_duration(&f, &t, "A,");
	parse_pitch_and_duration(&f, &t, "A");
	parse_pitch_and_duration(&f, &t, "a");
	parse_pitch_and_duration(&f, &t, "a'");
	parse_pitch_and_duration(&f, &t, "^A");
	parse_pitch_and_duration(&f, &t, "a'3");
	parse_pitch_and_duration(&f, &t, "a'3/2");
	parse_pitch_and_duration(&f, &t, "a'355/113");
	parse_pitch_and_duration(&f, &t, "a'/4");
	parse_pitch_and_duration(&f, &t, "z4");
	parse_pitch_and_duration(&f, &t, "C,,");
	parse_pitch_and_duration(&f, &t, "C,");
	parse_pitch_and_duration(&f, &t, "C");
	parse_pitch_and_duration(&f, &t, "c");
	parse_pitch_and_duration(&f, &t, "c'");
}

static void test_parser(void)
{
	char s[] = "zCDE FDEC G2c2B2C2";
	char *t = s;
	while (*t)
	{
		float f, l;
		t = parse_pitch_and_duration(&f, &l, t);
		fprintf(stderr, "f=%g l=%g\n", f, l);
	}
}

static void add_abc_chunk_into_score(
		struct wave_score *s, // output score
		char *a,              // input ABC string
		float b;              // input beats per minute of a unit note
		float t;              // time offset at the beginning
		)
{
	s->n = 0; // TODO: remove this line
	while (*a)
	{
		float ω; // note pitch
		float λ; // note length (in "abc" units)
		a = parse_pitch_and_duration(&ω, &λ, a);
		s->f[s->n] = ω;    // pitch
		s->t[s->n] = t;    // attack time
		s->l[s->n] = λ/b;  // duration
		s->i[s->n] = 1;    // instrument
		t += t + λ/b;
		s->n += 1;
	}
}


// this is the sole function that paints a brush stroke on the canvas
static void wave_play(                  // "play" note f between T[0] and T[1]
		struct wave_canvas *w,  // sound canvas
		struct wave_brush *b,   // instrument
		float f,                // fundamental frequency (note)
		float Ta, float Tb      // time interval
		)
{
	assert(1 == b->f[0]);
	int n[2] = { Ta * w->F, Tb * w->F }; // discrete time interval
	fprintf(stderr, "note %g from %g to %g (%d harmonics)\n",
			f, Ta, Tb, b->n);
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

// smoother, softer-sounding wave
static void wave_brush_init_smoother(struct wave_brush *b)
{
	b->n = 0x100;
	//b->n = 4;
	for (int i = 0; i < b->n; i++)
	{
		b->f[i] = i + 1;      // harmonic spectrum
		b->a[i] = 1/pow(b->f[i],2);  // inverse square frequency decay
	}
	b->λ = 0;
}
static void wave_brush_init_smoother3(struct wave_brush *b)
{
	b->n = 0x100;
	for (int i = 0; i < b->n; i++)
	{
		b->f[i] = i + 1;      // harmonic spectrum
		b->a[i] = 1/pow(b->f[i],2.1);  // inverse square frequency decay
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

static void test_pipeline(void)
{
	char a[] = "zCDE FDEC G2c2B2C2";  // the ABC string
	struct wave_score s[1];           // the wave score
	add_abc_chunk_into_score(s, a, 60, 0);
}


//static void wave_bwc

int main_no()
{
	//test_parse_note_name();
	test_get_pitch_and_duration();
	return 0;
}
int main_no2()
{
	test_parser();
	return 0;
}
int main_no3()
{
	test_pipeline();
	return 0;
}
int main_yes()
{
	struct wave_canvas w[1];
	wave_canvas_init(w, 7, 44000);

	struct wave_brush b[1];
	wave_brush_init_smoother(b);
	//wave_brush_init_smoother3(b);
	//wave_brush_init_triangular(b);
	//wave_brush_init_full(b);
	//wave_brush_init_pure(b);
	b->λ = 0.001;

	for (int i = 0; i <= 12; i++)
		wave_play(w, b, 440*pow(2,i/12.0), i/2.0, (i+1)/2.0);
	wave_play(w, b, 1*440, 0, 1);
	wave_play(w, b, 2*440, 1, 2);
	wave_play(w, b, 3*440, 2, 3);
	wave_play(w, b, 4*440, 3, 4);

	wave_quantized_stdout(w);
	return 0;
}
int main(){return main_no3();}
