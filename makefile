#  @(#)makefile	1.15 

SRC_DIR = src
SOURCES = \
	$(SRC_DIR)/chi.o \
	$(SRC_DIR)/ran1.o \
	$(SRC_DIR)/newentry2.o \
	$(SRC_DIR)/newentry.o \
	$(SRC_DIR)/file_inputnew.o \
	$(SRC_DIR)/fill_gai.o \
	$(SRC_DIR)/intfile6.o \
	$(SRC_DIR)/directmap2.o \
	$(SRC_DIR)/updat1.o \
	$(SRC_DIR)/writeterfinseg.o \
	$(SRC_DIR)/keep_it2.o \
	$(SRC_DIR)/ldmapseg.o \
	$(SRC_DIR)/reorderall.o \
	$(SRC_DIR)/intfile4.o \
	$(SRC_DIR)/reorder.o \
	$(SRC_DIR)/getfreqs.o \
	$(SRC_DIR)/intfileseg.o \
	$(SRC_DIR)/segments.o \
	$(SRC_DIR)/quicklike.o \
	$(SRC_DIR)/locus_list.o \
	$(SRC_DIR)/fast.o \
	$(SRC_DIR)/one_interval.o \
	$(SRC_DIR)/ludcmp.o \
	$(SRC_DIR)/update.o \
	$(SRC_DIR)/get_predl.o \
	$(SRC_DIR)/nrutil.o \
	$(SRC_DIR)/meankchi.o \
	$(SRC_DIR)/multallele.o \
	$(SRC_DIR)/runewt3.o \
	$(SRC_DIR)/bound.o \
	$(SRC_DIR)/metricho.o \
	$(SRC_DIR)/diplo.o \
	$(SRC_DIR)/func.o \
	$(SRC_DIR)/dfunc.o \
	$(SRC_DIR)/dfpmin.o \
	$(SRC_DIR)/lnsrch.o \
	$(SRC_DIR)/cal_SE2.o \
	$(SRC_DIR)/getpa.o \
	$(SRC_DIR)/jobin.o \
	$(SRC_DIR)/runewt2.o \
	$(SRC_DIR)/ldmapper.o

PROGRAM = ldmap

# gcc can take "-Wall -Werror" for checking
CFLAGS = -g -O
CPPFLAGS = $(INCLUDEDIRS)
LDFLAGS =
LINTFLAGS =
LIBRARIES = -lm

OBJECTS = $(SOURCES:.c=.o)

$(PROGRAM): $(OBJECTS)
	$(LINK.c) -o $@ $(OBJECTS) $(ULIBS) $(LIBRARIES)

.PHONY: clean

clean:
	$(RM) src/*.o ldmap
