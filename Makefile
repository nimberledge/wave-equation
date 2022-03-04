parwave: src/parwave.cpp src/Subdomain.h src/Subdomain.cpp src/initialconditions.h
	mpicxx -o bin/parwave src/parwave.cpp src/Subdomain.h src/Subdomain.cpp src/initialconditions.h src/config.h

clean:
	rm bin/*