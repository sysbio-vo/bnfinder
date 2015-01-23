ALL:
		rubber -d tut_art
		pdflatex tut
		mv tut_art.bbl tut.bbl
		pdflatex tut

clean:
		rubber --clean tut.tex
		rubber --clean tut_art.tex

clean_all:	
		rm tut.pdf
		rm tut_art.pdf
		rubber --clean tut.tex
		rubber --clean tut_art.tex