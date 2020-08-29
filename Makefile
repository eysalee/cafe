CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -march=native -fpie -fopenmp -g -O3 -DDEBUG -maes -DWITH_SYNC
LDFLAGS  = -L./libs -lpthread -lncomm -lsodium
INCFLAGS = -I./include -I./ncomm/include -I./ztk
SRCDIR = ./source

SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(SRCS:.cpp=.o)

ifeq ($(TIME_ME), 1)
	CXXFLAGS += -DWITH_TIMINGS
endif

ifeq ($(SYNC_ME), 1)
	CXXFLAGS += -DWITH_SYNC
endif

ifeq ($(AES_PRG), 1)
	CXXFLAGS += -DAES_PRG
endif

default: ncomm $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) main.cpp -o prog.exe $(OBJS) $(LDFLAGS)

ncomm:
	cd ncomm && make DEBUG=1
	rm -rf libs
	mkdir -p libs
	cp ncomm/libncomm.a libs/libncomm.a

ncomm-no-debug:
	cd ncomm && make DEBUG=0
	rm -rf libs
	mkdir -p libs
	cp ncomm/libncomm.a libs/libncomm.a

ncomm-clean:
	cd ncomm && make clean

shamir-ds-bm: shamir-simd-ds-bm shamir-flex-ds-bm shamir-regular-ds-bm

f2d-test: ncomm $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) f2d_test.cpp -o f2d_test.exe $(OBJS) $(LDFLAGS)

circuit-test: ncomm $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) circuit_test.cpp -o circuit_test.exe $(OBJS) $(LDFLAGS) -DTESTING

shamir_test: ncomm $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) shamir_test.cpp -o shamir_test.exe $(OBJS) $(LDFLAGS) -DTESTING

shamir-simd-ds-bm: ncomm-no-debug $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) shamir-simd-ds-bm.cpp -o shamir-simd-ds-bm.exe $(OBJS) $(LDFLAGS) -DWITH_TIMINGS

shamir-flex-ds-bm: ncomm-no-debug $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) shamir-flex-ds-bm.cpp -o shamir-flex-ds-bm.exe $(OBJS) $(LDFLAGS) -DWITH_TIMINGS

shamir-regular-ds-bm: ncomm-no-debug $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) shamir-regular-ds-bm.cpp -o shamir-regular-ds-bm.exe $(OBJS) $(LDFLAGS) -DWITH_TIMINGS

shamir-inner-to-flex-bm: ncomm-no-debug $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) shamir-inner-to-flex-bm.cpp -o shamir-inner-to-flex-bm.exe $(OBJS) $(LDFLAGS) -DWITH_TIMINGS -DTESTING

innerprod-local-comp: $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) innerprod_local_comp.cpp -o innerprod_local_comp.exe $(OBJS) $(LDFLAGS)

svm: ncomm $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) svm.cpp -o svm.exe $(OBJS) $(LDFLAGS)

clean:
	rm -rf libs
	rm -f prog.exe
	rm -f circuit_test.exe
	rm -f shamir_test.exe
	find $(SRCDIR) -iname "*.o" -delete

clean-all: clean ncomm-clean

.SUFFIXES: .cpp .o

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $< $(LDFLAGS)

.PHONY: ncomm ncomm-clean \
	clean clean-all \
	circuit-test shamir-test
