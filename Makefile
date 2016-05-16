
CXX? = g++
CFLAGS = -Wall -Wconversion -fPIC 
# -O3
SHVER = 1

all: svm-train svm-predict 

svm-predict: svm-predict.c solve_linear_system.o common.o Solver.o Solver_NU.o Solver_plus.o Cache.o Kernel.o
	$(CXX) $(CFLAGS) svm-predict.c solve_linear_system.o common.o Solver.o Solver_NU.o Solver_plus.o Cache.o Kernel.o -o svm-predict -lm
svm-train: svm-train.c solve_linear_system.o common.o Solver.o Solver_NU.o Solver_plus.o  Cache.o Kernel.o
	$(CXX) $(CFLAGS) svm-train.c solve_linear_system.o common.o Solver.o Solver_NU.o Solver_plus.o Cache.o Kernel.o -o svm-train -lm
solve_linear_system.o: solve_linear_system.c solve_linear_system.h 
	$(CXX) $(CFLAGS) -c solve_linear_system.c
common.o: common.c common.h Solver.h Solver_NU.h Solver_plus.h Kernel.h
	$(CXX) $(CFLAGS) -c common.c
Solver.o: Solver.cpp Solver.h Kernel.h common.h solve_linear_system.h
	$(CXX) $(CFLAGS) -c Solver.cpp
Solver_NU.o: Solver_NU.cpp Solver_NU.h Solver.h Kernel.h common.h
	$(CXX) $(CFLAGS) -c Solver_NU.cpp
Solver_plus.o: Solver_plus.cpp Solver_plus.h Kernel.h common.h solve_linear_system.h
	$(CXX) $(CFLAGS) -c Solver_plus.cpp
Kernel.o: Kernel.cpp Kernel.h Cache.h common.h
	$(CXX) $(CFLAGS) -c Kernel.cpp
Cache.o: Cache.cpp Cache.h common.h
	$(CXX) $(CFLAGS) -c Cache.cpp
clean:
	rm -f *~ common.o solve_linear_system.o Cache.o Kernel.o Solver.o Solver_NU.o Solver_plus.o svm-train svm-predict svm-scale
