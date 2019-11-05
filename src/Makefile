FLAGS=-Wall -g -O3 -std=c++11 -fopenmp

%.o: %.cpp
	g++ $(FLAGS) -c $^

gbkc: gbkc.o check-profiles.o count.o distance.o
	g++ $(FLAGS) -o $@ $^ -lz

clean:
	rm gbkc *.o
