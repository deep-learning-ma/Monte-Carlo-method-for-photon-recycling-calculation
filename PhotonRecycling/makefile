CXX = g++
ROOT_FLAGS = $(shell root-config --cflags --libs) -lMinuit
CXXFLAGS = -lpthread -lm -O3
objects = Material.o Simulation.o Run.o Run.exe

buildobj : $(objects) 

Material.o : Material.cpp Material.h
	$(CXX) $(ROOT_FLAGS) $(CXXFLAGS) -c $< -o $@

Simulation.o : Simulation.cpp Simulation.cpp Material.cpp Material.h
	$(CXX) $(ROOT_FLAGS) $(CXXFLAGS) -c $< -o $@

Run.o : Run.cpp Simulation.cpp Simulation.h Material.cpp Material.h
	$(CXX) $(ROOT_FLAGS) $(CXXFLAGS) -c $< -o $@

Run.exe : Run.o Simulation.o Material.o
	$(CXX) $^ $(ROOT_FLAGS) $(CXXFLAGS) -o $@

.PHONY : clean
clean :
	rm *.o *.exe
