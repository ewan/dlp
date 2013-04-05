all:
	cd java && $(MAKE) all;
	sh make_R_package.sh;

clean: 
	cd java && $(MAKE) clean;
	rm -rf bin;

install: all
	cd bin && R CMD INSTALL dlp
