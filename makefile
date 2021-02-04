#  @(#)makefile	1.15 

SOURCES= \
	chi.o ran1.o newentry2.o  newentry.o file_inputnew.o fill_gai.o intfile6.o directmap2.o updat1.o writeterfinseg.o keep_it2.o ldmapseg.o reorderall.o intfile4.o reorder.o getfreqs.o intfileseg.o segments.o quicklike.o locus_list.o fast.o one_interval.o ludcmp.o update.o get_predl.o nrutil.o meankchi.o multallele.o runewt3.o bound.o metricho.o diplo.o func.o dfunc.o dfpmin.o lnsrch.o cal_SE2.o getpa.o jobin.o runewt2.o ldmapper.o

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
