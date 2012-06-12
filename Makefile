OBJS = cut.o event.o options.o preselection.o preselectionCMG.o toolbox.o triggerinfo.o toolbox.o variableGetter.o lepton.o electron.o muon.o jet.o toolsCMG.o Dictionary.o
OUT = libHZZ2l2nu.a
HEADERS = cut.h event.h options.h preselection.h preselectionCMG.h toolbox.h triggerinfo.h variableGetter.h lepton.h electron.h muon.h jet.h toolsCMG.h
CC = g++
CFLAGS = -g -ansi -std=c++0x -Wall `root-config --cflags` -pg
LFLAGS = `root-config --libs` -pg
NAME = libHZZ2l2nu.a

$(NAME) : $(OBJS)
	ar rcs $(NAME) $(OBJS)

cut.o : cut.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c cut.cpp

event.o : event.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c event.cpp

options.o : options.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c options.cpp

preselection.o : preselection.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c preselection.cpp

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

toolsCMG.o : toolsCMG.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c toolsCMG.cpp

Dictionary.cpp : triggerinfo.h LinkDef.h
	rootcint -v -f Dictionary.cpp -c triggerinfo.h LinkDef.h
Dictionary.o : Dictionary.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c Dictionary.cpp -o Dictionary.o

clean :
	rm -f *.o $(NAME)
	rm -f Dictionary.*
