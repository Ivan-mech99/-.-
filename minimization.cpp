#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iomanip>
using namespace std;

//Безусловная минимизация функции многих переменных

double eps = 0.000001;
double tau = 0.1 * sqrt(eps);

double func(double* param)
{

	return pow(param[0] - 1, 2) + pow(param[1] - 2, 2);
}

double norm(double* param, int num_param)
{
	double s = 0;
	int i;
	for (i = 0; i < num_param; i++)
	{
		s += (param[i] * param[i]);
	}
	return sqrt(s);
}

double deriv1(double* param, int index)
{
	double result = 0;
	param[index] += tau;
	result += func(param);
	param[index] -= (2 * tau);
	result -= func(param);
	result /= (2 * tau);
	param[index] += tau;
	return result;
}

double deriv2(double* param, int ind1, int ind2)
{
	double result = 0;
	if (ind1 == ind2)
	{
		result -= (2 * func(param));
		param[ind1] += tau;
		result += func(param);
		param[ind1] -= (2 * tau);
		result += func(param);
		param[ind1] += tau;
		result /= (tau * tau);
	}
	else
	{
		result -= (2 * func(param));
		param[ind1] += tau;
		result += func(param);
		param[ind2] -= tau;
		result -= func(param);
		param[ind1] -= tau;
		result += func(param);
		param[ind2] += (2 * tau);
		result += func(param);
		param[ind1] -= tau;
		result -= func(param);
		param[ind2] -= tau;
		result += func(param);
		param[ind1] += tau;
		result /= (tau * tau * 2);
	}
	return result;
}

double deriv12(double* param, int ind1, int ind2)
{
	double result = 0;
	result -= (2 * func(param));
	param[ind1] += tau;
	result += func(param);
	param[ind2] -= tau;
	result -= func(param);
	param[ind1] -= tau;
	result += func(param);
	param[ind2] += (2 * tau);
	result += func(param);
	param[ind1] -= tau;
	result -= func(param);
	param[ind2] -= tau;
	result += func(param);
	param[ind1] += tau;
	result /= (tau * tau * 2);
	return result;
}

void grad(double* a, double* param, int m)
{
	int i;
	for (i = 0; i < m; i++)
	{
		a[i] = deriv1(param, i);
	}
}


void swap_strings(double** a, int ind1, int ind2, int m)
{
	int j; double temp;
	for (j = 0; j < m; j++)
	{
		temp = a[ind2][j];
		a[ind2][j] = a[ind1][j];
		a[ind1][j] = temp;
	}
}

void sum_mult(double** a, int ind1, int ind2, double coeff, int m)
{
	int j;
	for (j = 0; j < m; j++)
	{
		a[ind2][j] = a[ind2][j] + coeff * a[ind1][j];
	}
}

void show_matrix(double** a, int m)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			cout << a[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;
}

void unit_matrix(double** a, int m)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			a[i][j] = 0;
		}
	}
	for (i = 0; i < m; i++)
	{
		a[i][i] = 1;
	}
}

int not_zero(double** a, int ind, int start, int m)
{
	int i;
	cout << "Not zero active" << endl;
	cout << "ind: " << ind << endl;
	cout << "start: " << start << endl;
	for (i = start; i < m; i++)
	{
		if (a[ind][i] != 0)
		{
			return i;
		}
	}
	return -1;
}

void mult_string(double** a, int ind, double coeff, int m)
{
	int j;
	for (j = 0; j < m; j++)
	{
		a[ind][j] = a[ind][j] * coeff;
	}
}

void mult_column(double** a, int ind, int coeff, int m)
{
	int i;
	for (i = 0; i < m; i++)
	{
		a[i][ind] = a[i][ind] * coeff;
	}
}

void eq_2(double** a, double** b, int m)
{
	int i; int j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			a[i][j] = b[i][j];
		}
	}
}

void eq_1(double* a, double* b, int m)
{
	int i;
	for (i = 0; i < m; i++)
	{
		a[i] = b[i];
	}

}

void eq_minus(double* a, double* b, int m)
{
	int i;
	for (i = 0; i < m; i++)
	{
		a[i] = -b[i];
	}

}

void sum_1(double* a, double* b, int m)
{
	int i;
	for (i = 0; i < m; i++)
	{
		a[i] = a[i] + b[i];
	}
}

void dif_1(double* a, double* b, int m)
{
	int i;
	for (i = 0; i < m; i++)
	{
		a[i] = a[i] - b[i];
	}
}

void Hessian(double** a, double* param, int m)
{
	int i; int j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			a[i][j] = deriv2(param, i, j);
		}
	}
}

void mat_mult(double** a, double* b, double* rez, int m)
{
	int i; int j;
	for (i = 0; i < m; i++)
	{
		rez[i] = 0;
		for (j = 0; j < m; j++)
		{
			rez[i] += (a[i][j] * b[j]);
		}
	}
}

void vec_mult(double* a, double coeff, int m)
{
	int i;
	for (i = 0; i < m; i++)
	{
		a[i] = a[i] * coeff;
	}
}

void zero_vec(double* a, int m)
{
	int i;
	for (i = 0; i < m; i++)
	{
		a[i] = 0;
	}
}

void zero_mat(double** a, int m)
{
	int i; int j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			a[i][j] = 0;
		}
	}
}

bool cond1(double* x1, double* x0, double* check, int m, double E)
{
	eq_1(check, x1, m);
	dif_1(check, x0, m);
	double res = norm(check, m);
	if (res <= E)
	{
		return true;
	}
	return false;
}

bool cond2(double* x1, double* x0, int m, double E)
{
	double res1 = 0;
	double res2 = 0;
	res1 = func(x1);
	res2 = func(x0);
	res1 = res1 - res2;
	if (abs(res1) <= E)
	{
		return true;
	}
	return false;
}

bool cond3(double* x0, double* check, int m, double E)
{
	grad(check, x0, m);
	double res = norm(check, m);
	if (res <= E)
	{
		return true;
	}
	return false;
}

void show_vec(double* a, int m)
{
	int i;
	cout << "point: ";
	for (i = 0; i < m; i++)
	{
		cout << fixed << setprecision(6) << a[i] << ' ';
	}
	cout << endl;
}

void inverse(double** a, double** b, int m)
{
	unit_matrix(b, m);
	double hold;
	int k; int not_null; int i;
	for (k = 0; k < m; k++)
	{
		if (a[k][k] == 0)
		{
			not_null = not_zero(a, k, k + 1, m);
			cout << "not_null " << not_null << endl;
			swap_strings(a, k, not_null, m);
			swap_strings(b, k, not_null, m);
		}
		for (i = 0; i < m; i++)
		{
			if (i != k)
			{
				hold = -a[i][k] / a[k][k];
				sum_mult(a, k, i, hold, m);
				sum_mult(b, k, i, hold, m);
			}
		}
		hold = 1 / a[k][k];
		mult_string(a, k, hold, m);
		mult_string(b, k, hold, m);
	}
	eq_2(a, b, m);
}

double tangent(double* xx, double* hh, double* new_args, double a, double b, int m)
{
	zero_vec(new_args, m);
	double fa = 0;
	double fb = 0;
	double ffa = 0;
	double ffb = 0;
	double cn = 0;
	double ffc = 0;
	do {
		eq_1(new_args, hh, m);
		vec_mult(new_args, a, m);
		sum_1(new_args, xx, m);
		fa = func(new_args);
		int i;
		for (i = 0; i < m; i++)
		{
			ffa += hh[i] * deriv1(new_args, i);
		}
		eq_1(new_args, hh, m);
		vec_mult(new_args, b, m);
		sum_1(new_args, xx, m);
		fb = func(new_args);
		for (i = 0; i < m; i++)
		{
			ffb += hh[i] * deriv1(new_args, i);
		}
		cn = (fa - fb + ffb * b - ffa * a) / (ffb - ffa);
		eq_1(new_args, hh, m);
		vec_mult(new_args, cn, m);
		sum_1(new_args, xx, m);
		for (i = 0; i < m; i++)
		{
			ffc += hh[i] * deriv1(new_args, i);
		}
		if (ffc < 0)
		{
			a = cn;
		}
		else if (ffc > 0)
		{
			b = cn;
		}
	} while (abs(ffc) >= eps);
	return cn;
}

double fracture(double* xx, double* hh, double* new_args, int m)
{
	double lambda = 0.5;
	double mu = 2;
	double alpha = 2;
	zero_vec(new_args, m);
	eq_1(new_args, hh, m);
	vec_mult(new_args, alpha, m);
	sum_1(new_args, xx, m);
	if (!(func(new_args) < func(xx)))
	{
		while (!(func(new_args) < func(xx)))
		{
			alpha = lambda * alpha;
			zero_vec(new_args, m);
			eq_1(new_args, hh, m);
			vec_mult(new_args, alpha, m);
			sum_1(new_args, xx, m);
		}
		return alpha;
	}
	else
	{
		while ((func(new_args) < func(xx)))
		{
			alpha = mu * alpha;
			zero_vec(new_args, m);
			eq_1(new_args, hh, m);
			vec_mult(new_args, alpha, m);
			sum_1(new_args, xx, m);
		}
		return alpha;
	}
	return alpha;
}
int main()
{
	double a = 0;
	double b = 1.2;
	double atr = 0;
	int i; int M;
	int num_it = 0;
	double* h0;
	cout << "Enter dimensions: " << endl;
	cin >> M;
	h0 = new double[M];
	double* h00;
	h00 = new double[M];
	double* x0;
	x0 = new double[M];
	double* x1;
	x1 = new double[M];
	double* check;
	check = new double[M];
	double* hh;
	hh = new double[M];
	double* xx;
	xx = new double[M];
	double* new_args;
	new_args = new double[M];
	cout << "Enter x0: " << endl;
	for (i = 0; i < M; i++)
	{
		cin >> x0[i];
	}
	for (i = 0; i < M; i++)
	{
		x1[i] = x0[i];
	}
	do
	{
		eq_1(x0, x1, M);
		zero_vec(x1, M);
		grad(h0, x0, M);
		eq_1(xx, x0, M);
		eq_minus(hh, h0, M);
		atr = fracture(xx, hh, new_args, M);
		vec_mult(h0, atr, M);
		dif_1(x1, h0, M);
		sum_1(x1, x0, M);
		num_it++;
	} while (!(cond1(x1, x0, check, M, sqrt(eps)) & cond2(x1, x0, M, sqrt(eps)) & cond3(x0, check, M, sqrt(eps))));
	double** hesse;
	hesse = new double* [M];
	for (i = 0; i < M; i++)
	{
		hesse[i] = new double[M];
	}
	zero_mat(hesse, M);
	double** help_mat;
	help_mat = new double* [M];
	for (i = 0; i < M; i++)
	{
		help_mat[i] = new double[M];
	}
	zero_mat(help_mat, M);
	eq_1(x0, x1, M);
	cout << "number of iterations: " << num_it << endl;
	show_vec(x1, M);
	cout << "f(x_k): " << func(x1) << endl;
	cout << endl;
	num_it = 0;
	do
	{
		eq_1(x0, x1, M);
		zero_vec(x1, M);
		Hessian(hesse, x0, M);
		inverse(hesse, help_mat, M);
		grad(h0, x0, M);
		mat_mult(hesse, h0, h00, M);
		eq_1(xx, x0, M);
		eq_minus(hh, h00, M);
		atr = tangent(xx, hh, new_args, a, b, M);
		vec_mult(h00, atr, M);
		dif_1(x1, h00, M);
		sum_1(x1, x0, M);
		num_it++;
	} while (!(cond1(x1, x0, check, M, eps) & cond2(x1, x0, M, eps) & cond3(x0, check, M, eps)));
	cout << "number of iterations: " << num_it << endl;
	show_vec(x1, M);
	cout << "f(x_k): " << fixed << setprecision(6) << func(x1);
	delete[] h0;
	delete[] h00;
	delete[] x0;
	delete[] x1;
	delete[] xx;
	delete[] hh;
	delete[] check;
	delete[] new_args;
	for (i = 0; i < M; i++)
	{
		delete[] hesse[i];
	}
	delete[] hesse;
	for (i = 0; i < M; i++)
	{
		delete[] help_mat[i];
	}
	delete[] help_mat;
	return 0;
}

