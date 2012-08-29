OBJS = \
	   cut.o \
	   electron.o \
	   event.o \
	   jet.o \
	   lepton.o \
	   muon.o \
	   options.o \
	   photon.o \
	   preselectionCMG.o \
	   toolbox.o \
	   toolbox.o \
	   toolsCMG.o \
	   triggerinfo.o \
	   variableGetter.o \
	   Dictionary.o
OUT = libHZZ2l2nu.a
HEADERS = \
		  cut.h \
		  electron.h \
		  event.h \
		  jet.h \
		  lepton.h \
		  muon.h \
		  options.h \
		  photon.h \
		  preselectionCMG.h \
		  toolbox.h \
		  toolsCMG.h \
		  triggerinfo.h \
		  variableGetter.h
CC = g++
CFLAGS = -O3 -ansi -std=c++0x -Wall `root-config --cflags`
LFLAGS = `root-config --libs`
NAME = libHZZ2l2nu.a

$(NAME) : $(OBJS)
	ar rcs $(NAME) $(OBJS)

cut.o : cut.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c cut.cpp

event.o : event.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c event.cpp

options.o : options.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c options.cpp

preselectionCMG.o : preselectionCMG.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c preselectionCMG.cpp

toolbox.o : toolbox.cpp $(HEADERS) 
	$(CC) $(CFLAGS) -c toolbox.cpp

triggerinfo.o : triggerinfo.cpp $(HEADERS) 
	$(CC) $(CFLAGS) -c triggerinfo.cpp

variableGetter.o : variableGetter.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c variableGetter.cpp

lepton.o : lepton.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c lepton.cpp

electron.o : electron.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c electron.cpp

muon.o : muon.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c muon.cpp

jet.o : jet.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c jet.cpp

photon.o : photon.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c photon.cpp

toolsCMG.o : toolsCMG.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c toolsCMG.cpp

Dictionary.cpp : triggerinfo.h LinkDef.h
	rootcint -v -f Dictionary.cpp -c triggerinfo.h LinkDef.h
Dictionary.o : Dictionary.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c Dictionary.cpp -o Dictionary.o

clean :
	rm -f *.o $(NAME)
	rm -f Dictionary.*
