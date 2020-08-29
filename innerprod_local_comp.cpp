#include "ztk.hpp"
#include "timing.hpp"
#include <vector>
#include <iostream>
#include <unistd.h>

using namespace ztk;
using namespace std;

template<size_t K>
void innerprod(size_t n) {

    cout << "K=" << K << "\n";

    vector<GR<K, 4>> a0;
    vector<GR<K, 4>> a1;
    vector<GR<K, 4>> b0;
    vector<GR<K, 4>> b1;

    for (size_t i = 0; i < n; i++) {
	a0.emplace_back(GR<K, 4>{Z2k<K>{3}});
	a1.emplace_back(GR<K, 4>{Z2k<K>{3}});

    }

    for (size_t i = 0; i < (n/2)-1; i++) {
	b0.emplace_back(GR<K, 4>{{Z2k<K>{3}, Z2k<K>{3}, Z2k<K>::zero, Z2k<K>::zero}});
	b1.emplace_back(GR<K, 4>{{Z2k<K>::zero, Z2k<K>::zero, Z2k<K>{3}, Z2k<K>{3}}});
    }

    GR<K, 4> v;
    GR<K, 4> w;

    {
	TIMER_BEGIN(naive);
	w = a0[0] * a1[0];

	for (size_t i = 1; i < n; i++) {
	    w += (a0[i] * a1[i]);
	}
	TIMER_END(naive);

    }

    {
	TIMER_BEGIN(innerProd);
	v = b0[0] * b1[0];
	for (size_t i = 0; i < n/2; i++) {
	    v += (b0[i] * b1[i]);
	}
	TIMER_END(innerProd);
    }
}

template<size_t K>
void matmul(size_t n, size_t m) {

    // TODO: encoding stuff with innerprod ftso. matrix multiplications is a bit
    // more annoying than I had initially thought.

    // this makes the innerprod encoding easier
    assert(n%2 == 0);
    assert(m%2 == 0);

    cout << "K=" << K << ", n=" << n << ", m=" << m << "\n";

    vector<vector<GR<K, 4>>> A0 (n);
    vector<vector<GR<K, 4>>> A1 (m);
    vector<vector<GR<K, 4>>> B0 (n);
    vector<vector<GR<K, 4>>> B1 (m);

    cout << "init A0/B0\n";
    int c = 0;
    for (size_t i = 0; i < n; i++) {
	A0[i].resize(m);
	if (i % 2 == 0)
	    B0[i/2].resize(m/2);
	for (size_t j = 0; j < m; j++) {
	    A0[i][j] = GR<K,4>{Z2k<K>{c}};
	    if (j % 2 == 0 && i % 2 == 0) {
	    	B0[i/2][j/2] = GR<K,4>{{Z2k<K>{c}, Z2k<K>{c}, Z2k<K>::zero, Z2k<K>::zero}};
		c++;
	    }
	}
    }

    cout << "init A1/B1\n";
    c = 10000;
    for (size_t j = 0; j < m; j++) {
	A1[j].resize(n);
	if (j % 2 == 0)
	    B1[j/2].resize(n/2);
	for (size_t i = 0; i < n; i++) {
	    A1[j][i] = GR<K,4>{Z2k<K>{c}};
	    if (i % 2 == 0 && j % 2 == 0) {
	    	B1[j/2][i/2] = GR<K,4>{{Z2k<K>::zero, Z2k<K>::zero, Z2k<K>{c}, Z2k<K>{c}}};
		c++;
	    }
	}
    }

    cout << "init AA\n";
    vector<vector<GR<K, 4>>> AA (n);
    for (size_t i = 0; i < n; i++) {
	AA[i].resize(n);
    }

    cout << "init BB\n";
    vector<vector<GR<K, 4>>> BB (n/2);
    for (size_t i = 0; i < n/2; i++)
    	BB[i].resize(n/2);

    cout << "naive\n";
    TIMER_BEGIN(naive);
    for (size_t i = 0; i < n; i++) {
	for (size_t j = 0; j < n; j++) {
	    AA[i][j] = A0[i][0]*A1[0][j];
	    for (size_t k = 1; k < m; k++) {
		AA[i][j] += A0[i][k]*A1[k][j];
	    }

	}
    }
    TIMER_END(naive);

    cout << "innerprod\n";
    TIMER_BEGIN(innerprod);
    for (size_t i = 0; i < n/2; i++) {
	for (size_t j = 0; j < n/2; j++) {
	    BB[i][j] = B0[i][0]*B1[0][j];
	    for (size_t k = 1; k < m/2; k++) {
		BB[i][j] += B0[i][k]*B1[k][j];
	    }
	}
    }
    TIMER_END(innerprod);

    for (size_t i = 0; i < n; i++) {
	for (size_t j = 0; j < n; j++) {
	    cout << AA[i][j][0] << " ";
	}
	cout << "\n";
    }

    cout << "\n\n";

    for (size_t i = 0; i < n/2; i++) {
	for (size_t j = 0; j < n/2; j++) {
	    cout << BB[i][j][3] << " ";
	}
	cout << "\n";
    }
}

int main() {

    size_t n;

#define repeat(__n) for (size_t i = 0; i < (__n); i++)

    // {
    // 	n = 10000;
    // 	cout << "inner product n = " << n << "\n";
    // 	repeat(50) {
    // 	    innerprod<64>(n);
    // 	    sleep(1);
    // 	}

    // 	// matmul<64>(10, 10);
    // }

    cout << "\n\n";
    // return 0;

    {
	n = 100000;
	cout << "inner product n = " << n << "\n";
	repeat(50) {
	    innerprod<64>(n);
	    sleep(1);
	}
    }

    cout << "\n\n";

    // {
    // 	n = 1000000;
    // 	cout << "inner product n = " << n << "\n";
    // 	repeat(100) {
    // 	    innerprod<64>(n);
    // 	    sleep(1);
    // 	}
    // }

    cout << "\n\n";

    // {
    // 	n = 10000000;
    // 	cout << "inner product n = " << n << "\n";
    // 	repeat(50) {
    // 	    innerprod<64>(n);
    // 	    sleep(1);
    // 	}
    // }

    cout << "\n\n";

    // {
    // 	n = 50000000;
    // 	cout << "inner product n = " << n << "\n";
    // 	repeat(50) {
    // 	    innerprod<64>(n);
    // 	    sleep(1);
    // 	}
    // }


}
