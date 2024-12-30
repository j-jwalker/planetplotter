all: final_project

final_project: planet_visualization.o central_pressure_est.o
	g++ -o final_project planet_visualization.o central_pressure_est.o -lpgplot -lcpgplot -lX11 -lm

planet_visualization.o: planet_visualization.cpp final_project_header.h
	g++ -c planet_visualization.cpp -Wall -O2 -std=c++11

central_pressure_est.o: central_pressure_est.cpp final_project_header.h
	g++ -c central_pressure_est.cpp -Wall -O2 -std=c++11

# Clean build files
clean:
	rm -f final_project planet_visualization.o central_pressure_est.o
