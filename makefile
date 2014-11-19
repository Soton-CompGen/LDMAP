#  @(#)makefile	1.15 

SOURCES= \
	simboot.o file_inputg2.o ran1.o intfile7.o getfreqsh2.o getfreqsh.o holes.o intfile6.o directmap2.o updat1.o writeterfinseg.o keep_it2.o ldmapseg.o reorderall.o intfile4.o reorder.o writeterfin.o writeter.o printinputmap.o getfreqs.o intfileseg.o segments.o quicklike.o hapallele2.o diplo1.o multallele1.o haplo.o hapallele.o printmap.o locus_list.o fast.o directmap.o one_interval.o pool.o ludcmp.o keep_it.o update.o final.o get_predl.o nrutil.o ldmap3.o fill_str2.o mergeloc.o meankchi.o multallele.o file_inputg.o runewt3.o ldmap2.o bound.o intfile3.o metricho.o diplo.o file_input2.o func.o dfunc.o dfpmin.o lnsrch.o cal_SE2.o getpa.o jobin.o runewt2.o ldmapper.o

PROGRAM= ldmapper1 

# gcc can take "-Wall -Werror" for checking
CFLAGS= -g -O
CPPFLAGS= $(INCLUDEDIRS)
LDFLAGS=
LINTFLAGS=
LIBRARIES= -lm

OBJECTS= $(SOURCES:.c=.o)

$(PROGRAM): $(OBJECTS)
	$(LINK.c) -o $@ $(OBJECTS) $(ULIBS) $(LIBRARIES)

.PHONY: clean

clean:
	$(RM) *.o ldmapper1
