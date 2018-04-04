#pour compiler le programme en entier : make all
#pour compiler un sous-programme : make sous-programme.o

run : main.cc read_data.cpp maillage.cpp diffusion.cpp plic.cpp
	g++ -std=c++11 -I EigenLibrary/Eigen -o run main.cc read_data.cpp maillage.cpp diffusion.cpp plic.cpp recul.cpp

all : main.o diffusion.o maillage.o plic.o read_data.o recul.o
	g++ -std=c++11 -I EigenLibrary/Eigen -o run build/main.o build/diffusion.o build/maillage.o build/plic.o build/read_data.o build/recul.o

main.o : main.cc read_data.h maillage.h diffusion.h plic.h recul.h
	g++ -std=c++11 -I EigenLibrary/Eigen -c -o build/main.o main.cc

diffusion.o : diffusion.cpp diffusion.h
	g++ -std=c++11 -I EigenLibrary/Eigen -c -o build/diffusion.o diffusion.cpp

maillage.o : maillage.cpp maillage.h
	g++ -std=c++11 -I EigenLibrary/Eigen -c -o build/maillage.o maillage.cpp

plic.o : plic.cpp plic.h
	g++ -std=c++11 -I EigenLibrary/Eigen -c -o build/plic.o plic.cpp

read_data.o : read_data.cpp read_data.h
	g++ -std=c++11 -I EigenLibrary/Eigen -c -o build/read_data.o read_data.cpp

recul.o : recul.cpp recul.h
	g++ -std=c++11 -I EigenLibrary/Eigen -c -o build/recul.o recul.cpp

recul3D.o : recul3D.cpp recul3D.h
	g++ -std=c++11 -I EigenLibrary/Eigen -c -o build/recul3D.o recul3D.cpp
