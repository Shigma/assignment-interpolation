#include <iostream>
#include <string>
#include <stdlib.h>
#include <cmath>

using namespace std;

#define pm(x) (1 - (((x) & 1) << 1))

void interpolation(int N, long double max_error[], bool output = false) {
	int N1 = N + 1;
	long double source[N1], target[N1];
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

	// newton
	// "ad" stands for "average deviation"
	long double newton[N1], ad[N1];

	// initialize ad
	for (int i = 0; i <= N; i++) {
		ad[i] = source[i];
	}

	// calculate ad, O(n^2)
	for (int j = 0; j < N; ++j) {
		for (int i = N - 1; i >= j; --i) {
			ad[i + 1] = (ad[i + 1] - ad[i]) / (j + 1.);
		}
	}

	// get answers
	for (int j = 1; j <= N; ++j) {
		newton[j] = 0;
		for (int i = 0; i <= N; ++i) {
			newton[j] = newton[j] * (j - .5 - N + i) + ad[N - i];
		}
	}

	// lagrange
	long double lagrange[N1], fact[N1], fact2[N1];

	// calculate factorials
	fact[0] = 1; fact2[0] = 1;
	for (int i = 1; i <= N; ++i) {
		fact[i] = fact[i - 1] * i;
		fact2[i] = fact2[i - 1] * (i - .5);
	}

	// get answers
	for (int j = 1; j <= N; ++j) {
		lagrange[j] = 0;
		for (int i = 0; i <= N; ++i) {
			lagrange[j] += source[i]
				* fact2[j]
				* fact2[N - j + 1]
				* pm(i - j + 1)
				/ fact[i]
				/ fact[N - i]
				/ (j - 0.5 - i);
		}
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
	// get answers
	for (int i = 1; i <= N; ++i) {
		spline[i] = linear[i] - unit / 8 * (m[i] - m[i - 1]);
	}

	if (output) {
		string filename = "out_";
		char number[5] = "";
		itoa(N, number, 10);
		filename += number;
		filename += ".csv";
		freopen(filename.data(), "w", stdout);
		cout << "lagrange,newton,linear,spline" << endl;
	}

	max_error[0] = 0;
	max_error[1] = 0;
	max_error[2] = 0;
	max_error[3] = 0;
	for (int i = 1; i <= N; ++i) {
		long double error_lagrange = abs(target[i] - lagrange[i]);
		long double error_newton = abs(target[i] - newton[i]);
		long double error_linear = abs(target[i] - linear[i]);
		long double error_spline = abs(target[i] - spline[i]);

		max_error[0] = max_error[0] > error_lagrange ? max_error[0] : error_lagrange;
		max_error[1] = max_error[1] > error_newton ? max_error[1] : error_newton;
		max_error[2] = max_error[2] > error_linear ? max_error[2] : error_linear;
		max_error[3] = max_error[3] > error_spline ? max_error[3] : error_spline;

		if (output) {
			cout << abs(target[i] - lagrange[i]) << ","
					 << abs(target[i] - newton[i]) << ","
					 << abs(target[i] - linear[i]) << ","
					 << abs(target[i] - spline[i]) << endl;
		}
	}

	if (output) fclose(stdout);
}

int main() {
	cout.precision(8);
	cout.setf(ios::scientific);

	long double errors[10][4];
	interpolation(10, errors[0]);
	interpolation(20, errors[1], true);
	interpolation(30, errors[2], true);
	interpolation(40, errors[3]);
	interpolation(50, errors[4]);
	interpolation(60, errors[5]);
	interpolation(70, errors[6]);
	interpolation(80, errors[7]);
	interpolation(90, errors[8]);

	freopen("errors.csv", "w", stdout);
	cout << "lagrange,newton,linear,spline" << endl;

	for (int i = 0; i < 9; ++i) {
		for (int j = 0; j < 4; ++j) {
			if (j) cout << ",";
			cout << errors[i][j];
		}
		cout << endl;
	}

	fclose(stdout);
}
