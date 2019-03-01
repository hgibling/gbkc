FLAGS=-Wall -g -O3

%.o: %.cpp
	g++ $(FLAGS) -c $^

gbkc: gbkc.o
	g++ $(FLAGS) -o $@ $^ -lz

clean:
	rm gbkc *.o
