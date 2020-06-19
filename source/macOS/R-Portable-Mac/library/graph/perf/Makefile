BASE = graphperf
R ?= R-devel

all: $(BASE).pdf
	@echo "building perf report using $(R)"

$(BASE).pdf: $(BASE).tex
	$(R) CMD texi2dvi --pdf -c -q $<

$(BASE).tex: $(BASE).Rnw
	$(R) CMD Sweave $<

clean:
	@rm -fv $(BASE).pdf $(BASE).tex
