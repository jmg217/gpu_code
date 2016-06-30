all: executable

executable: inner_control_mesh.o meshgen.o meshestimator.o pathestimator.o plot_high_bias.o Payoff.o 
	g++ inner_control_mesh.o meshgen.o meshestimator.o pathestimator.o plot_high_bias.o Payoff.o -o executable

meshgen.o: meshgen.cpp
	g++ -c meshgen.cpp -std=c++0x

meshestimator.o: meshestimator.cpp
	g++ -c meshestimator.cpp -std=c++0x

pathestimator.o: pathestimator.cpp
	g++ -c pathestimator.cpp -std=c++0x

inner_control_mesh.o: inner_control_mesh.cpp
	g++ -c inner_control_mesh.cpp -std=c++0x

Payoff.o: Payoff.cpp
	g++ -c Payoff.cpp -std=c++0x

plot_high_bias.o: plot_high_bias.cpp
	g++ -c plot_high_bias.cpp -std=c++0x

clean:
	rm -rf *o mesh

