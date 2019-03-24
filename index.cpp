#include <iostream>
#include <string>
#include <stdlib.h>
#include <cmath>

using namespace std;

#define N 20
#define N1 21
#define pm(x) (1 - (((x) & 1) << 1))

int main() {
	long double source[N1], target[N1];
	long double EXP2N = exp2l(N);
	long double unit = 6. / N;

	for (int i = 0; i <= N; ++i) {
		source[i] = expl(-2 * unit * i);
		if (i) {
			target[i] = expl(-unit * (2 * i - 1));
		}
	}

	// linear
	long double linear[N1];
	for (int i = 1; i <= N; ++i) {
		linear[i] = (source[i] + source[i - 1]) / 2;
	}

	// lagrange
	long double lagrange[N1], frac[N1], frac2[N1], coef[N1];
	frac[0] = 1; frac2[0] = 1;
	for (int i = 1; i <= N; ++i) {
		frac[i] = frac[i - 1] * i;
		frac2[i] = frac2[i - 1] * (2 * i - 1);
	}
	for (int i = 0; i <= N; ++i) {
		coef[i] = source[i] / frac[i] / frac[N - i] * pm(N - i);
	}
	for (int i = 1; i <= N; ++i) {
		lagrange[i] = 0;
		for (int j = 0; j <= N; ++j) {
			lagrange[i] += coef[j] / (2 * (i - j) - 1);
		}
		lagrange[i] *= frac2[i] * frac2[N - i] * pm(N - i + 1) / EXP2N;
	}

	// spline
	long double spline[N1];
	long double coefb = 1.5 / unit;
	long double A[N1], B[N1], a[N1], b[N1], m[N1];
	a[0] = 0; b[0] = -4;
	a[N] = 1; b[N] = -4 * target[N];
	for (int i = 1; i < N; ++i) {
		a[i] = 0.5;
		b[i] = (source[i + 1] - source[i - 1]) * coefb;
	}
	A[0] = 0; B[0] = -2;
	for (int i = 1; i <= N; ++i) {
		long double deno = 2 + (1 - a[i]) * A[i - 1];
		A[i] = -a[i] / deno;
		B[i] = (b[i] - (1 - a[i]) * B[i - 1]) / deno;
	}
	m[N] = B[N];
	for (int i = N - 1; i >= 0; --i) {
		m[i] = A[i] * m[i + 1] + B[i];
	}
	for (int i = 1; i <= N; ++i) {
		spline[i] = linear[i] - unit / 8 * (m[i] - m[i - 1]);
	}

	cout.precision(8);
	cout.setf(ios::scientific);

	string filename = "out";
	filename += itoa(N);
	filename += ".csv";
	freopen(filename.data(), "w", stdout);
	cout << "lagrange,linear,spline" << endl;

	for (int i = 1; i <= N; ++i) {
		cout << abs(target[i] - lagrange[i]) << ","
				 << abs(target[i] - linear[i]) << ","
				 << abs(target[i] - spline[i]) << endl;
	}

	fclose(stdout);
}
