all: run

EXEC_NAME=solver
CMP=g++-12
CMP_KEYS=-Wall -Wextra -Wfloat-equal -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual -pedantic
CMP_OPT=-O3

run: build $(wildcard data.dat)
	./$(EXEC_NAME) 20

build: $(wildcard *.cpp *.h Makefile)
	$(CMP) -std=c++20 $(CMP_KEYS) $(CMP_OPT) *.cpp -o $(EXEC_NAME)

clear:
	rm -rf fft result.dat abs.dat
