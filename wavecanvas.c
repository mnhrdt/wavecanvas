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
// TODO: find better "painting" names
// 	orchestra    <-> palette/brush set?
// 	score        <-> drawing instructions? recipe?
//
//   Play a score using an orchestra :: draw a painting using these brushes


#include <assert.h>
#include <math.h>     // sin, exp, fabs, fmax
#include <stdlib.h>   // malloc, free
#include <stdint.h>   // uint8_t uint16_t
#include <stdio.h>    // fwrite, stdout (plus fprintf, stderr for debugging)
#include <ctype.h>    // isalpha, tolower


static const float π = M_PI;
static const float SAMPLING_RATE_IN_HERTZ = 44100;


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


#define MAX_NOTES 40000 // should be enough for anybody
struct wave_score {          // a set of notes (in no particular order)
	int n;               // total number of notes
	float f[MAX_NOTES];  // pitch (fundamental frequency)
	float t[MAX_NOTES];  // start of attack (time)
	float l[MAX_NOTES];  // length
};

static void debug_score(struct wave_score *s)
{
	fprintf(stderr, "score with %d notes\n", s->n);
	for (int i = 0; i < s->n; i++)
		fprintf(stderr,"\tf=%g t=%g l=%g\n", s->f[i], s->t[i], s->l[i]);
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
// z  : rest

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
	// note: this loop also eats leading spaces and extraneous characters
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
	default: fprintf(stderr, "ERROR: unrecognized note '%c' (%d)\n",
				 a[i], a[i]);
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
	int c = 0;   // chord state
	float C = 0; // held voices state
	while (*a)
	{
		float ω; // note pitch
		float λ; // note length (in "abc" units)
		if (isspace(*a)) { a += 1; continue; }   // eat spaces
		if (*a == '[') { c = 1; a += 1; }        // begin chord
		if (*a == '{') { C = t; a += 1; }        // begin holding voices
		if (*a == ';') { t = C; a += 1; }        // next voice
		a = parse_pitch_and_duration(&ω, &λ, a);
		s->f[s->n] = ω;    // pitch
		s->t[s->n] = t;    // attack time
		s->l[s->n] = λ/b;  // duration
		if (!c) t += λ/b;  // advance counter if outside chord
		if (*a == ']') { c = 0; a += 1; }        // end chord
		if (*a == '}') { C = 0; a += 1; }        // end hold
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

static float decay_hermite(float λ, float t)
{
	return t*exp(-λ * t);
}

static float decay_sigma(float λ, float t)
{
	return tanh(200*t)*exp(-λ * t);
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
			float ξ = 0;///3;//k*1.273+0.3*b->n - f * (1+Ta) / (1+k+Tb*Tb);

			if (f * b->f[k] * 2 < w->F) // avoid aliasing
				x += A * b->a[k] * sin(f * b->f[k] * 2*π*t + ξ);
		}
		//w->x[j] += exp(- b->λ * f * t) * x;
		//w->x[j] += x * decay_exp(b->λ * f, t);
		//w->x[j] += x * decay_exp(b->λ*100, t);
		//w->x[j] += x * decay_planck(b->λ, t);
		w->x[j] += x * decay_sigma(b->λ*f, t);
	}
}

static void wave_play_score_using_single_instrument(
		struct wave_canvas *c,    // canvas to fill
		struct wave_score *s,     // list of notes to play
		struct wave_brush *b      // instruments to play the notes with
		)
{
	for (int i = 0; i < s->n; i++)
	{
		float f = s->f[i];
		float t = s->t[i];
		float T = s->l[i] + t;
		f = 1.00*(f - 261.63) + 261.63;
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
	b->n = 10;
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
	b->n = 10;
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
		//if (i>0) b->f[i] -= 0.04;
		b->a[i] = T[i];//1/pow(b->f[i],2.1);  // inverse square frequency decay
	}
	b->λ = 0.005;
}

static void wave_triangular_filter(struct wave_canvas *w, int d)
{
	float *y = malloc(w->n * sizeof*y);
	for (int i = 0; i < w->n - d; i++)
	{
		y[i] = w->x[i];
		for (int j = 0; j < d; j++)
			y[i] += (d - j)*w->x[i+j]/d;
	}
	free(w->x);
	w->x = y;
}

//#include "random.c"
static void wave_quantized_stdout(struct wave_canvas *w)
{
	wave_triangular_filter(w, 40);
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

// prelude in C, unit = sixteenth note
static char *bwv_846 =
	// CEGce       CDAdf         B,DGdf       CEGce
	// CEAea       CD^FAd        B,DGdg       B,CEGc
	// A,CEGc      D,A,D^Fc      G,B,DGB      G,_B,EG^c
	// F,A,DGd     F,_A,DFB      E,G,CGc      E,F,A,CF
	// D,F,A,CF    G,,D,G,B,F    C,E,G,CE     C,G,_B,CE
	// F,,F,A,CE   ^F,,C,A,C_E   _A,,F,B,CD   G,,F,G,B,D
	// G,,E,G,CE   G,,D,G,CF     G,,D,G,B,F   G,,_E,A,C^F
	// G,,E,G,CG   G,,D,G,CF     G,,D,G,B,F   C,,C,G,_B,E

	// cat wavecanvas.c|sed -e '1,/bwv_846/ d'|head -n 8|cut -c5-|tr -s ' '|tr ' ' '\n'|sed 's/\([_^]*[[:alpha:]][,]*\)/\t\1/g'|awk '{print "{" $1 "8;z" $2 "7;zz" $3$4$5$3$4$5 "}"}'|awk '{print "\t\""$0, $0"\""}'
	"{C8;zE7;zzGceGce} {C8;zE7;zzGceGce}"
	"{C8;zD7;zzAdfAdf} {C8;zD7;zzAdfAdf}"
	"{B,8;zD7;zzGdfGdf} {B,8;zD7;zzGdfGdf}"
	"{C8;zE7;zzGceGce} {C8;zE7;zzGceGce}"
	"{C8;zE7;zzAeaAea} {C8;zE7;zzAeaAea}"
	"{C8;zD7;zz^FAd^FAd} {C8;zD7;zz^FAd^FAd}"
	"{B,8;zD7;zzGdgGdg} {B,8;zD7;zzGdgGdg}"
	"{B,8;zC7;zzEGcEGc} {B,8;zC7;zzEGcEGc}"
	"{A,8;zC7;zzEGcEGc} {A,8;zC7;zzEGcEGc}"
	"{D,8;zA,7;zzD^FcD^Fc} {D,8;zA,7;zzD^FcD^Fc}"
	"{G,8;zB,7;zzDGBDGB} {G,8;zB,7;zzDGBDGB}"
	"{G,8;z_B,7;zzEG^cEG^c} {G,8;z_B,7;zzEG^cEG^c}"
	"{F,8;zA,7;zzDGdDGd} {F,8;zA,7;zzDGdDGd}"
	"{F,8;z_A,7;zzDFBDFB} {F,8;z_A,7;zzDFBDFB}"
	"{E,8;zG,7;zzCGcCGc} {E,8;zG,7;zzCGcCGc}"
	"{E,8;zF,7;zzA,CFA,CF} {E,8;zF,7;zzA,CFA,CF}"
	"{D,8;zF,7;zzA,CFA,CF} {D,8;zF,7;zzA,CFA,CF}"
	"{G,,8;zD,7;zzG,B,FG,B,F} {G,,8;zD,7;zzG,B,FG,B,F}"
	"{C,8;zE,7;zzG,CEG,CE} {C,8;zE,7;zzG,CEG,CE}"
	"{C,8;zG,7;zz_B,CE_B,CE} {C,8;zG,7;zz_B,CE_B,CE}"
	"{F,,8;zF,7;zzA,CEA,CE} {F,,8;zF,7;zzA,CEA,CE}"
	"{^F,,8;zC,7;zzA,C_EA,C_E} {^F,,8;zC,7;zzA,C_EA,C_E}"
	"{_A,,8;zF,7;zzB,CDB,CD} {_A,,8;zF,7;zzB,CDB,CD}"
	"{G,,8;zF,7;zzG,B,DG,B,D} {G,,8;zF,7;zzG,B,DG,B,D}"
	"{G,,8;zE,7;zzG,CEG,CE} {G,,8;zE,7;zzG,CEG,CE}"
	"{G,,8;zD,7;zzG,CFG,CF} {G,,8;zD,7;zzG,CFG,CF}"
	"{G,,8;zD,7;zzG,B,FG,B,F} {G,,8;zD,7;zzG,B,FG,B,F}"
	"{G,,8;z_E,7;zzA,C^FA,C^F} {G,,8;z_E,7;zzA,C^FA,C^F}"
	"{G,,8;zE,7;zzG,CGG,CG} {G,,8;zE,7;zzG,CGG,CG}"
	"{G,,8;zD,7;zzG,CFG,CF} {G,,8;zD,7;zzG,CFG,CF}"
	"{G,,8;zD,7;zzG,B,FG,B,F} {G,,8;zD,7;zzG,B,FG,B,F}"
	"{C,,8;zC,7;zzG,_B,EG,_B,E} {C,,8;zC,7;zzG,_B,EG,_B,E}"

	"{C,,8C,,8;zC,15;zzF,A,CFCA,CA,F,A,F,D,F,D,}"
	"{C,,8C,,8;zB,,15;zzGBdfdBdBGBDFED}"
	"[C,,16C,16E16G16c16]"
;

static void test_score(void)
{
	char *x = bwv_772_stimme1;
	char *y = bwv_772_stimme2;
	struct wave_score s[1];           // the wave score
	s->n = 0;
	float tx = add_abc_chunk_into_score(s, x, 80*4, 0);
	float ty = add_abc_chunk_into_score(s, y, 80*4, 0);
	fprintf(stderr, "tx=%g ty=%g\n", tx, ty);
	//debug_score(s);

	struct wave_canvas w[1];
	wave_canvas_init(w, 66, SAMPLING_RATE_IN_HERTZ);

	struct wave_brush b[1];
	wave_brush_init_smoother3(b);
	b->λ = 0.005; // decay rate

	wave_play_score_using_single_instrument(w, s, b);

	wave_quantized_stdout(w);
}

static void test_chords(void)
{
	char *x = bwv_846;
	struct wave_score s[1];
	s->n = 0;
	float t = add_abc_chunk_into_score(s, x, 90*4, 0);

	struct wave_canvas w[1];
	wave_canvas_init(w, t+4, SAMPLING_RATE_IN_HERTZ);
	//debug_score(s);

	struct wave_brush b[1];
	wave_brush_init_smoother3(b);
	b->λ = 0.005; // decay rate

	wave_play_score_using_single_instrument(w, s, b);

	wave_quantized_stdout(w);
}

static void test_waveplay(void)
{
	struct wave_canvas w[1];
	wave_canvas_init(w, 7, SAMPLING_RATE_IN_HERTZ);

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
	test_chords();
	return 0;
}

int main(){return main_yes();}
