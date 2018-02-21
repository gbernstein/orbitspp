KPL/MK

A "furnsh kernel" for use with SPICE software, commanding input of
the kernel (=input data) files of use for TNO orbit fitting.
The PATH_VALUES needs to be changed to site-specific location of the repo.

\begindata

PATH_VALUES = ('/Users/garyb/CODE/orbitspp/data')
PATH_SYMBOLS= ('DATA')
KERNELS_TO_LOAD = ( '$DATA/naif0012.tls',
                    '$DATA/de430_tno.bsp',
		    '$DATA/pck0010.tpc',
		    '$DATA/earth_latest_high_prec.bpc' )

\begintext
