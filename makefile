run : main.cc read_data.cpp maillage.cpp diffusion.cpp plic.cpp
	g++ -std=c++11 -I EigenLibrary/Eigen -o run main.cc read_data.cpp maillage.cpp diffusion.cpp plic.cpp recul.cpp
