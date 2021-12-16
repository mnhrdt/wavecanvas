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

static void debug_score(struct wave_score *s)
{
	fprintf(stderr, "score with %d notes\n", s->n);
	for (int i = 0; i < s->n; i++)
		fprintf(stderr, "\tf=%g t=%g l=%g i=%d\n",
				s->f[i], s->t[i], s->l[i], s->i[i]);
}

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
	//fprintf(stderr, "parsing note string \"%s\"...", a);
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
	//fprintf(stderr, "\tn=%d o=%d d=%d t=%g\n", *n, *o, *d, *t);
	return i;
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

// return a pointer to the remaning part of the string
static char *parse_pitch_and_duration(float *f, float *t, char *a)
{
	int o; // octave (0=octave that contains middle C)
	int n; // note within the octave (from 0=c to 6=b)
	int d; // accidental increment
	int i = parse_note_name(&o, &n, &d, t, a);
	float x = chromatic_note(n, d) + 12*o - 9;
	*f = n >= 0 ? 440 * pow(2, x/12) : 0; // equal temperament
	//fprintf(stderr, "gpd \"%s\" x=%g f=%g t=%g\n", a, x, *f, *t);
	return a + i;
}

static float add_abc_chunk_into_score(
		struct wave_score *s, // output score
		char *a,              // input ABC string
		float b,              // input beats per minute of a unit note
		float t               // time offset at the beginning
		) // return value = end of timespan
{
	b /= 60;
	int c = 0; // chord state
	while (*a)
	{
		float ω; // note pitch
		float λ; // note length (in "abc" units)
		if (isspace(*a)) { a += 1; continue; }   // eat spaces
		if (*a == '[') { c = 1; a += 1; }        // begin chord
		a = parse_pitch_and_duration(&ω, &λ, a);
		s->f[s->n] = ω;    // pitch
		s->t[s->n] = t;    // attack time
		s->l[s->n] = λ/b;  // duration
		s->i[s->n] = 0;    // instrument
		if (!c) t += λ/b;  // advance counter if outside chord
		if (*a == ']') { c = 0; a += 1; }       // end chord
		s->n += 1;
	}
	return t;
}

static float decay_exp(float λ, float t)
{
	return exp(-λ * t);
}

static float decay_planck(float λ, float t)
{
	λ *= 0.01;
	t *= 0.001;
	return (t*t*t/(exp(t/λ)-1));
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
	//fprintf(stderr, "note %g from %g to %g (%d harmonics)\n",
	//		f, Ta, Tb, b->n);
	float φ = 20; // least note
	float A = (220-φ)/(f-φ); // volume equalizer
	for (int i = 0; i < n[1] - n[0]; i++)
	{
		int j = i + n[0]; // index in the w->n table
		if (j < 0 || j >= w->n)
		{fprintf(stderr,"ξ");continue;} // avoid outside of canvas

		float t = i / w->F;  // time since the start of the note

		float x = 0;        // accumulator for the j-th sample
		for (int k = 0; k < b->n; k++)
		{
			// random phase
			float ξ = k*1.273+0.3*b->n - f * (1+Ta) / (1+k+Tb*Tb);

			if (f * b->f[k] * 2 < w->F) // avoid aliasing
				x += A * b->a[k] * sin(f * b->f[k] * 2*π*t + ξ);
		}
		//w->x[j] += exp(- b->λ * f * t) * x;
		//w->x[j] += x * decay_exp(b->λ * f, t);
		w->x[j] += x * decay_exp(b->λ*100, t);
		//w->x[j] += x * decay_planck(b->λ, t);
	}
}

static void wave_play_score_using_orchestra(
		struct wave_canvas *c,    // canvas to fill
		struct wave_score *s,     // list of notes to play
		struct wave_orchestra *o  // instruments to play the notes with
		)
{
	for (int i = 0; i < s->n; i++)
	{
		struct wave_brush *b = o->t + s->i[i];
		float f = s->f[i];
		float t = s->t[i];
		float T = s->l[i] + t;
		wave_play(c, b, f, t, T);
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
	b->λ = 0.001;
}

static void wave_brush_init_smoother3(struct wave_brush *b)
{
	float T[] = {1, .5, .35, 0.25, 0.15, 0.1, 0.05};
	b->n = sizeof T / sizeof *T;
	for (int i = 0; i < b->n; i++)
	{
		b->f[i] = (i + 1);
		//if (i>0) b->f[i] -= 0.06;
		b->a[i] = T[i];//1/pow(b->f[i],2.1);  // inverse square frequency decay
	}
	b->λ = 0.005;
}

//#include "random.c"
static void wave_quantized_stdout(struct wave_canvas *w)
{
	float M = 0;
	for (int i = 0; i < w->n; i++)
		M = fmax(M, fabs(w->x[i]));
	int16_t *x = malloc(w->n * sizeof*x);
	for (int i = 0; i < w->n; i++)
		x[i] = 8000*((w->x[i]/M));// + 0.00002*random_cauchy());
	fwrite(x, sizeof*x, w->n, stdout);
	free(x);
}

// inventio 1, unit = sixteenth note, no trills, separate voices
// pagination by Nels Drue's "urtext" from IMSLP 70230
static char *bwv_772_stimme1 =
	"zCDE FDEC G2c2B2c2       dGAB cABG d2g2f2g2"
	"eagf egfa gfed cedf      edcB AcBd cBAG ^FAGB"
	"A2D2 c3d BAG^F EG^FA     GBAc Bdce dB/2c/2dg B2AG"
	"G2zz zzzz zGAB cABG      ^F2zz zzzz zABc dBcA"
	"B2zz zzzz zdcB AcBd      c2zz zzzz zedc Bd^ce"

	"d2^c2d2e2 f2A2B2c2       d2^F2^G2A2 B2c2 d5"
	" E^F^G A^F^GE edce dcBd  ca^gb aefd ^Gfed c2BA"
	"Aagf egfa g9              efg afge f9"
	" gfe dfeg f9              def gefd e9"
	" cde fdec defg afge      fgab c'abg c'2g2 e2dc"
	"c_BAG FEG_B ABCE DcFB    [E16G16c16]"
	;
static char *bwv_772_stimme2 =
	"zzzz zzzz zC,D,E, F,D,E,C,            G,2G,,2 zzzz zG,A,B, CA,B,G,"
	"C2B,2C2D2 E2G,2A,2B,2                 C2E,2^F,2G,2 A,2B,2 C5"
	" D,E,^F, G,E,^F,D, G,2B,,2C,2D,2      E,2^F,2G,2E,2 B,,3C, D,2D,,2"
	"zG,,A,,B,, C,A,,B,,G,, D,2G,2^F,2G,2  A,D,E,^F, G,E,^F,D, A,2D2C2D2"
	"G,GFE DFEG F2E2F2D2                   EAGF EGFA G2F2G2E2"

	"F_BAG FAG_B AGFE DFEG                 FEDC B,DCE DCB,A, ^G,B,A,C"
	"B,2E,2 D3E CB,A,G, ^F,A,^G,B,         A,CB,D CEDF E2A,2E2E,2"
	"A,2A,,2 zzzz zEDC B,D^CE              D9 A,B,C DB,CA,"
	"B,9 DCB, A,CB,D                       C9 G,A,_B, CA,B,G,"
	"A,2_B,2A,2G,2 F,2D2C2_B,2             A,2F2E2D2 ED,E,F, G,E,F,D,"
	"E,2C,2D,2E,2 F,D,E,F, G,2G,,2         [C,,16C,16]"
	;

static void test_score(void)
{
	//char a[] = "zCDE FDEC G2c2B2c2";  // the ABC string
	char *a = bwv_772_stimme1;
	char *b = bwv_772_stimme2;
	struct wave_score s[1];           // the wave score
	s->n = 0;
	float ta = add_abc_chunk_into_score(s, a, 80*4, 0);
	float tb = add_abc_chunk_into_score(s, b, 80*4, 0);
	fprintf(stderr, "ta=%g tb=%g\n", ta, tb);
	//debug_score(s);

	struct wave_canvas w[1];
	wave_canvas_init(w, 66, 44000);

	// setup an "orchestra" with one instrument
	struct wave_orchestra o[1];
	o->n = 1;
	wave_brush_init_smoother3(o->t + 0);

	wave_play_score_using_orchestra(w, s, o);

	wave_quantized_stdout(w);
}

static void test_waveplay(void)
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
}

int main_yes()
{
	test_score();
	return 0;
}

int main(){return main_yes();}
