FLAGS=-Wall -g -O3 -std=c++11

%.o: %.cpp
	g++ $(FLAGS) -c $^

gbkc: gbkc.o
	g++ $(FLAGS) -o $@ $^ -lz

clean:
	rm gbkc *.o
