all: executable

executable: convert_to_array.o cuda_density.o cuda_meshestimator.o cuda_meshgen.o cudapathestimator.o index.o inner_control_mesh.o cuda_meshweights.o Payoff.o plot_high_bias.o 
	nvcc -G -lineinfo convert_to_array.o cuda_density.o cuda_meshestimator.o cuda_meshgen.o cudapathestimator.o index.o inner_control_mesh.o cuda_meshweights.o Payoff.o plot_high_bias.o -o executable

cuda_meshgen.o: cuda_meshgen.cpp
	g++ -c cuda_meshgen.cpp -std=c++0x 

convert_to_array.o: convert_to_array.cpp
	g++ -c convert_to_array.cpp

cuda_density.o: cuda_density.cpp
	g++ -c cuda_density.cpp

cuda_meshestimator.o: cuda_meshestimator.cpp
	g++ -c cuda_meshestimator.cpp 

cudapathestimator.o: cudapathestimator.cu
	nvcc -G -lineinfo -c cudapathestimator.cu -arch=sm_30

index.o: index.cpp
	g++ -c index.cpp

inner_control_mesh.o: inner_control_mesh.cpp
	g++ -c inner_control_mesh.cpp 

Payoff.o: Payoff.cpp
	g++ -c Payoff.cpp -std=c++0x

plot_high_bias.o: plot_high_bias.cpp
	g++ -c plot_high_bias.cpp -std=c++0x

cuda_meshweights.0: cuda_meshweights.cpp
	g++ -c cuda_meshweights.cpp

clean:
	rm -rf *o mesh

