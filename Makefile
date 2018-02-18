#	$Id: Makefile,v 2.2 2002/01/22 03:18:41 garyb Exp $	
# Make file for orbit fitting programs.  Change flags below as required.

CC = gcc
CFLAGS = -O
LIBS = -lm
BINDIR = ~/bin.linux/

PROGS = fit_radec predict abg_to_aei postdict scatter AllMPC scenario2 planner
ALL:  $(PROGS)

clean:
	rm *.o *.dvi *.aux *.log core $(PROGS)
install:
	mv -f $(PROGS) $(BINDIR)
tar:
	tar -cvf orbfit.tar *.c *.h Makefile TODO README observatories.dat \
	 qb1.a* pluto14yr.a* update.tex
OBJS  = fit_radec.o predict.o abg_to_aei.o abg_to_xyz.o

NR   = nrutil.o ludcmp.o lubksb.o gaussj.o mrqmin_orbit.o mrqcof_orbit.o \
	covsrt.o gasdev.o ran1.o
ORBSUBS = orbfit1.o transforms.o dms.o ephem_earth.o


ephem_earth.o: ephem_read.h
orbfit1.o ephem_earth.o: ephem_types.h
$(OBJS) $(ORBSUBS) orbfit2.o mrqmin_orbit.o mrqcof_orbit.o :  orbfit.h Makefile

# Note that aeiderivs is very slow to compile optimized so I'll give
# it its own compilation here.
aeiderivs.o: orbfit.h
	$(CC) -O0 -c aeiderivs.c -o aeiderivs.o

fit_radec: fit_radec.o orbfit2.o $(ORBSUBS) $(NR)
	$(CC) fit_radec.o orbfit2.o $(ORBSUBS) $(NR) $(LIBS) -o fit_radec
predict: predict.o $(ORBSUBS) $(NR)
	$(CC) predict.o $(ORBSUBS) $(NR) $(LIBS) -o predict
postdict: postdict.o $(ORBSUBS) $(NR)
	$(CC) postdict.o $(ORBSUBS) $(NR) $(LIBS) -o postdict
abg_to_aei: abg_to_aei.o $(ORBSUBS)  aeiderivs.o $(NR)
	$(CC)  abg_to_aei.o aeiderivs.o $(ORBSUBS) $(NR) $(LIBS) -o abg_to_aei
abg_to_xyz: abg_to_xyz.o $(ORBSUBS) $(NR)
	$(CC)  abg_to_xyz.o $(ORBSUBS) $(NR) $(LIBS) -o abg_to_xyz
scenario: scenario.o orbfit2.o aeiderivs.o $(ORBSUBS) $(NR)
	$(CC) scenario.o orbfit2.o aeiderivs.o \
	$(ORBSUBS) $(NR) $(LIBS) -o scenario
scenario2: scenario2.o orbfit2.o aeiderivs.o $(ORBSUBS) $(NR)
	$(CC) scenario2.o orbfit2.o aeiderivs.o \
	$(ORBSUBS) $(NR) $(LIBS) -o scenario2
scatter: scatter.o orbfit2.o aeiderivs.o $(ORBSUBS) $(NR)
	$(CC) scatter.o orbfit2.o aeiderivs.o \
	$(ORBSUBS) $(NR) $(LIBS) -o scatter
AllMPC: AllMPC.o orbfit2.o aeiderivs.o $(ORBSUBS) $(NR)
	$(CC) AllMPC.o orbfit2.o aeiderivs.o \
	$(ORBSUBS) $(NR) $(LIBS) -o AllMPC
planner: planner.o $(ORBSUBS) $(NR)
	$(CC) planner.o $(ORBSUBS) $(NR) $(LIBS) -o planner

# Programs for debugging tests ( am not maintaining all these)
check_posn: check_posn.o orbfit1.o transforms.o dms.o ephem_earth.o $(NR)
	$(CC) check_posn.o orbfit1.o transforms.o dms.o ephem_earth.o $(NR) \
	$(LIBS) -o check_posn

check_orbit: check_orbit.o orbfit1.o transforms.o dms.o ephem_earth.o $(NR)
	$(CC) check_orbit.o orbfit1.o transforms.o dms.o ephem_earth.o $(NR) \
	$(LIBS) -o check_orbit
fit_amoeba: fit_amoeba.o orbfit1.o transforms.o dms.o ephem_earth.o \
	$(NR) amoeba.o amotry.o	
	$(CC) fit_amoeba.o orbfit1.o \
	transforms.o dms.o ephem_earth.o $(NR) amoeba.o amotry.o\
	$(LIBS) -o fit_amoeba
testang: testang.o orbfit1.o transforms.o dms.o ephem_earth.o $(NR)
	$(CC) testang.o orbfit1.o transforms.o dms.o ephem_earth.o $(NR) \
	$(LIBS) -o testang


# HST-planning programs
hststuff1: hststuff1.o $(ORBSUBS) $(NR)
	$(CC) hststuff1.o $(ORBSUBS) $(NR) $(LIBS) -o hststuff1
hststuff2: hststuff2.o $(ORBSUBS) $(NR)
	$(CC) hststuff2.o $(ORBSUBS) $(NR) $(LIBS) -o hststuff2
hststuff3: hststuff3.o $(ORBSUBS) $(NR)
	$(CC) hststuff3.o $(ORBSUBS) $(NR) $(LIBS) -o hststuff3
hststuff4: hststuff4.o $(ORBSUBS) aeiderivs.o orbfit2.o $(NR)
	$(CC) hststuff4.o $(ORBSUBS) aeiderivs.o orbfit2.o $(NR) $(LIBS) -o hststuff4
radec_to_invar: radec_to_invar.o $(ORBSUBS) $(NR)
	$(CC) radec_to_invar.o $(ORBSUBS) $(NR) $(LIBS) -o radec_to_invar
invdump: invdump.o $(ORBSUBS) $(NR)
	$(CC) invdump.o $(ORBSUBS) $(NR) $(LIBS) -o invdump
testeph: testeph.o $(ORBSUBS) $(NR)
	$(CC) testeph.o $(ORBSUBS) $(NR) $(LIBS) -o testeph
plutinos: plutinos.o $(ORBSUBS) $(NR)
	$(CC) plutinos.o $(ORBSUBS) $(NR) $(LIBS) -o plutinos
