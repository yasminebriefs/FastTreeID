FILES := main.cpp ../src/identification.cpp ../src/algebra.cpp ../src/random.cpp

OPENSSL := -lcrypto -DHAVE_OPENSSL
GMP := -lgmp -DHAVE_GMP

default: with-openssl

with-openssl:
	g++ -std=c++20 -O2 $(FILES) $(OPENSSL) -o FastTreeID

with-gmp:
	g++ -std=c++20 -O2 $(FILES) $(GMP) -o FastTreeID

no-prime-check:
	g++ -std=c++20 -O2 $(FILES) -o FastTreeID


clean:
	rm FastTreeID
