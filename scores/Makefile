PDF   = $(shell ls *.ly | sed 's/\.ly/\.pdf/')
all   : $(PDF)
%.pdf : %.ly ; lilypond -s --pdf $^
clean :      ; $(RM) $(PDF)
