# Pipeline to build Penalized time-varying COX model

#setssh:
#	git remote set-url origin git@github.com:CYGUBICKO/pcoxtime-pkg.git

######################################################################

## Create and install the package
build-package:
	R CMD build .

check-win-release:
	echo "devtools::check_win_release('.')" | R --slave

check-win-dev:
	echo "devtools::check_win_devel('.')" | R --slave

install-package:
	R CMD INSTALL pcoxtime_1.0.4.*

check-package:
	echo "devtools::check('.')" | R --slave

update-export:
	echo "Rcpp::compileAttributes('.')" | R --slave

update-doc:
	echo "devtools::document('.')" | R --slave

install:
	make update-export && make update-doc && make check-package && make build-package && make install-package

update-oldrepo:
	cp -r DESCRIPTION man NAMESPACE R src tests ../../promote_science/learning_codes/pcoxtime/pcoxtime/

######################################################################
clean: 
	find . \( -name "\.#*" -o -name "*~" -o -name ".Rhistory" \) -exec rm {} \;

######################################################################

