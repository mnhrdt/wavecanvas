\version "2.22"

\header {
	title = "couronnée d'étoiles"
	composer = "M. Dannaud"
	tagline = ""
}

\transpose a g
\relative
{
	\clef treble
	\key a \minor
	\time 4/4
	\override Lyrics.LyricSpace.minimum-distance = #3.0
	a'8 b c a c4. d8 |
	b( a) b g e4 e4  |
	%\break
	a8 b c a c4 c8 c |
	b a b c b4. r8   |
	%\break
	c8 d e c c4 c8 c |
	f e d c d4 r8 d  |
	%\break
	d e c b a4. c8   |
	b a g e a2       \bar "||"
	\break
	e'8 e16 e16~ 16 f16 e8 d4 d4   |
	c8 c16 c16~ 16 d16 c8 b2       \bar "||"
	%\break
	a8 a16 b16~ 16 a16 b8 c4 c8 a8 |
	d d16 d16~ 8 c8 d2             \bar "||"
	%\break
	e8 e16 e16~ 16 f16 e8 d4 d4    |
	c8 c16 c16~ 16 d16 c8 b2       \bar "||"
	%\break
	a8 a16 b16~ 16 a16 b8 c4. c8   |
	b a g e a2                     \bar ":|."
}
\addlyrics
{
	Nous te sa -- lu -- ons
	ô toi, No -- tre Da -- me
	Ma -- rie Vier -- ge Sain -- te que
	dra -- pe le so -- leil
	cou -- ron -- née d'é -- toi -- les, la
	lune est sous tes pas, en
	toi nous est don -- née l'au --
	ro -- re du sa -- lut.
	Ma -- rie, È -- ve nou -- vel -- le
	et joie de ton Sei -- gneur,
	tu as don né nais -- san -- ce à
	Jé -- sus le Sau -- veur
	Par toi nous sont ou -- ver -- tes
	les por -- tes du jar -- din
	Gui -- de nous en che -- min
	É -- toi -- le du ma -- tin
}

% vim:set ai tw=77 filetype=lilypond:
