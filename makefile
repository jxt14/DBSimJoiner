main: main.o SimJoiner.o
	g++ main.o SimJoiner.o -o main -g

SimJoiner.o:SimJoiner.cpp SimJoiner.h
	g++ SimJoiner.cpp -c
main.o:main.cpp SimJoiner.h
	g++ main.cpp -c

clean:
	rm main.o SimJoiner.o main