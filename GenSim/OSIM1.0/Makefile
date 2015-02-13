#############################################################################
# Makefile for building: OSIM1.0

INCPATH       = -I/Users/erxingfangshui/Dropbox/CODE/boost \
		-I/Users/erxingfangshui/Dropbox/CODE/eigen3
CXXFLAGS      = -O3

OSIM1.0:main.o cohort.o animal_class.o genome_info.o
	g++ -o OSIM1.0 main.o cohort.o animal_class.o genome_info.o

main.o: main.cpp cohort.h
	g++ -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp 

cohort.o: cohort.cpp cohort.h animal_class.h
	g++ -c $(CXXFLAGS) $(INCPATH) -o cohort.o cohort.cpp 

animal_class.o: animal_class.cpp animal_class.h genome_info.h
	g++ -c $(CXXFLAGS) $(INCPATH) -o animal_class.o animal_class.cpp 

genome_info.o:genome_info.cpp genome_info.h
	g++ -c $(CXXFLAGS) $(INCPATH) -o genome_info.o genome_info.cpp 

clean:
	\rm *.o OSIM1.0


