include ../../Makefile.config
all:
	make Math_
	make NBody_
	make KDTree_
	make Analysis_
	make Cosmology_
	make InitCond_

clean: Mathc_ NBodyc_ KDTreec_ Analysisc_ Cosmologyc_ InitCondc_ libandincc_

Math_: ../../Makefile.config
	cd Math; make
NBody_: ../../Makefile.config
	cd NBody; make
KDTree_: ../../Makefile.config
	cd KDTree; make
Cosmology_: ../../Makefile.config
	cd Cosmology; make
InitCond_: ../../Makefile.config
	cd InitCond; make
Analysis_: ../../Makefile.config
	cd Analysis; make

Mathc_:
	cd Math; make clean
NBodyc_:
	cd NBody; make clean
KDTreec_:
	cd KDTree; make clean
Cosmologyc_:
	cd Cosmology; make clean
InitCondc_:
	cd InitCond; make clean
Analysisc_:
	cd Analysis; make clean

libandincc_:
	rm -f $(NBODYINCLUDEDIR)/*.h
	rm -f $(NBODYLIBDIR)/*.a
