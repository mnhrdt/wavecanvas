\version "2.22"

\header {
	title = "couronnée d'étoiles"
	composer = "M. Dannaud"
}

%\transpose a g
\relative
{
	\clef treble
	\key a \minor
	\time 4/4
	a'8 b c a c4. d8 |
	b( a) b g e4 e4  |
	a8 b c a c4 c8 c |
	b a b c b4. r8   |
	c8 d e c c4 c8 c |
	f e d c d4 r8 d  |
	d e c b a4. c8   |
	b a g e a2       \bar "||"
	e'8 e16 e16~ 16 f16 e8 d4 d4   |
	c8 c16 c16~ 16 d16 c8 b2       \bar "||"
	a8 a16 b16~ 16 a16 b8 c4 c8 a8 |
	d d16 d16~ 8 c8 d2             \bar "||"
}


% vim:set ai tw=77 filetype=lilypond:
