cc wavecanvas.c -lm &&
	./a.out > x.raw &&
	sox -v 1.5 -r 44100 -e signed-integer -b 16 -c 1 x.raw x.wav &&
	sox x.wav x2.wav reverb 45 95 lowpass -2 2300 pitch 0 tempo 1.0 &&
	aplay x2.wav

