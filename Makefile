OBJS = \
	   cut.o \
	   electron.o \
	   event.o \
	   eventPrinter.o \
	   jet.o \
	   lepton.o \
	   muon.o \
	   options.o \
	   photon.o \
	   photonPrescale.o \
	   preselectionCMG.o \
	   toolbox.o \
	   toolbox.o \
	   toolsCMG.o \
	   triggerinfo.o \
	   variableGetter.o \
	   RooZPtPdf.o \
	   Dictionary.o
OUT = libHZZ2l2nu.a
HEADERS = \
		  cut.h \
		  electron.h \
		  event.h \
		  eventPrinter.h \
		  jet.h \
		  lepton.h \
		  muon.h \
		  options.h \
		  photon.h \
		  photonPrescale.h \
		  preselectionCMG.h \
		  toolbox.h \
		  toolsCMG.h \
		  triggerinfo.h \
		  variableGetter.h \
		  RooZPtPdf.h
CC = g++
CFLAGS = -ansi -std=c++0x -Wall -O3 `root-config --cflags` -I$(ROOFITSYS)/include
LFLAGS = `root-config --libs` -O3 -ansi -std=c++0x -Wall
NAME = libHZZ2l2nu.a

$(NAME) : $(OBJS)
	ar rcs $(NAME) $(OBJS)

cut.o : cut.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c cut.cpp

event.o : event.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c event.cpp

eventPrinter.o : eventPrinter.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c eventPrinter.cpp

options.o : options.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c options.cpp

preselectionCMG.o : preselectionCMG.cpp $(HEADERS) eventPrinter.h
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

photonPrescale.o : photonPrescale.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c photonPrescale.cpp

toolsCMG.o : toolsCMG.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c toolsCMG.cpp

RooZPtPdf.o : RooZPtPdf.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c RooZPtPdf.cpp

Dictionary.cpp : triggerinfo.h RooZPtPdf.h LinkDef.h
	rootcint -v -f Dictionary.cpp -c -I$(ROOFITSYS)/include triggerinfo.h RooZPtPdf.h LinkDef.h
Dictionary.o : Dictionary.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c Dictionary.cpp -o Dictionary.o

clean :
	rm -f *.o $(NAME)
	rm -f Dictionary.*
