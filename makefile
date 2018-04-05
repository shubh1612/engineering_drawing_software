generate2D_v2.0: generate2D_v2.0.o 
	g++ -o generate2D_v2.0.exe generate2D_v2.0.o

generate2D_v2.0.o: generate2D_v2.0.cpp
	g++ -c generate2D_v2.0.cpp -std=c++11
clean:
	rm -f generate2D_v2.0 generate2D_v2.0.o