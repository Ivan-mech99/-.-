#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <iomanip>

using namespace std;

double eps = 0.000001;

double func(double x)
{
	double result = 0;
	result = 1 - x * (exp(x) - exp(-x));                    //5.0 / 6.0 * x;                         //1-x*(exp(x)-exp(-x));
	return result;
}

double K(double x, double t)
{
	double result = 0;
	result = x * x * exp(x * t);                //1.0 / 2.0 * x * t;                             //x*x*exp(x*t);
	return result;
}

double ideal_solution(double x)
{
	double result = 0;
	result = 1;                   //x;                                     //1
	return result;
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

void swap_strings(double* a, int ind1, int ind2)
{
	double temp;
	temp = a[ind2];
	a[ind2] = a[ind1];
	a[ind1] = temp;
}

void mult_string(double** a, int ind, double coeff, int m)
{
	int j;
	for (j = 0; j < m; j++)
	{
		a[ind][j] = a[ind][j] * coeff;
	}
}

void mult_string(double* a, int ind, double coeff)
{
	a[ind] = a[ind] * coeff;
}

void sum_mult(double** a, int ind1, int ind2, double coeff, int m)
{
	int j;
	for (j = 0; j < m; j++)
	{
		a[ind2][j] = a[ind2][j] + coeff * a[ind1][j];
	}
}

void sum_mult(double* a, int ind1, int ind2, double coeff)
{
	a[ind2] = a[ind2] + coeff * a[ind1];
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

void show(double* a, int m)
{
	int i;
	for (i = 0; i < m; i++)
	{
		cout << fixed << setprecision(7) << a[i] << ' ';
	}
	cout << endl;
}

void system_build(double** a, double* fv, double aa, double bb, int m)
{
	int i, j;
	double step = (bb - aa) / (m - 1);
	double* points = new double[m];
	for (i = 0; i < m; i++)
	{
		points[i] = i * step + aa;
	}
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			if (j % 2 == 0)
			{
				a[i][j] = -2.0 / 6.0 * K(points[i], points[j]) * step * 2.0;
			}
			else
			{
				a[i][j] = -4.0 / 6.0 * K(points[i], points[j]) * step * 2.0;
			}
		}
		a[i][0] = -1.0 / 6.0 * K(points[i], points[0]) * step * 2.0;
		a[i][m - 1] = -1.0 / 6.0 * K(points[i], points[m - 1]) * step * 2.0;
		a[i][i] += 1.0;
	}
	for (i = 0; i < m; i++)
	{
		fv[i] = func(points[i]);
	}
	delete[] points;
}

int system_solve(double** a, double* fv, int m)
{
	int i, k; int not_null; double hold;
	for (k = 0; k < m; k++)
	{
		if (a[k][k] == 0)
		{
			not_null = not_zero(a, k, k + 1, m);
			if (not_null == -1)
				return -1;
			cout << "not_null " << not_null << endl;
			swap_strings(a, k, not_null, m);
			swap_strings(fv, k, not_null);
		}
		for (i = 0; i < m; i++)
		{
			if (i != k)
			{
				hold = -a[i][k] / a[k][k];
				sum_mult(a, k, i, hold, m);
				sum_mult(fv, k, i, hold);
			}
		}
		hold = 1 / a[k][k];
		mult_string(a, k, hold, m);
		mult_string(fv, k, hold);
	}
	return 0;
}


void eq_1(double* a, double* b, int m)
{
	int i;
	for (i = 0; i < m; i++)
	{
		a[i] = b[i];
	}

}

double new_func(double* u1, int s1, double* u2, int s2, double x, double aa, double bb)
{
	int i, j;
	double result = 0;
	double* points1 = new double[s1];
	double step1 = (bb - aa) / (s1 - 1);
	double* points2 = new double[s2];
	double step2 = (bb - aa) / (s2 - 1);
	double* coeffs1 = new double[s1];
	double* coeffs2 = new double[s2];
	for (i = 0; i < s1; i++)
	{
		points1[i] = i * step1 + aa;
	}
	for (i = 0; i < s2; i++)
	{
		points2[i] = i * step2 + aa;
	}
	for (i = 0; i < s1; i++)
	{
		if (i % 2 == 0)
		{
			coeffs1[i] = 2.0 / 6.0 * step1 * K(x, points1[i]) * 2.0;
		}
		else
		{
			coeffs1[i] = 4.0 / 6.0 * step1 * K(x, points1[i]) * 2.0;
		}
	}
	coeffs1[0] = 1.0 / 6.0 * step1 * K(x, points1[0]) * 2.0;
	coeffs1[s1 - 1] = 1.0 / 6.0 * step1 * K(x, points1[s1 - 1]) * 2.0;
	for (i = 0; i < s2; i++)
	{
		if (i % 2 == 0)
		{
			coeffs2[i] = 2.0 / 6.0 * step2 * K(x, points2[i]) * 2.0;
		}
		else
		{
			coeffs2[i] = 4.0 / 6.0 * step2 * K(x, points2[i]) * 2.0;
		}
	}
	coeffs2[0] = 1.0 / 6.0 * step2 * K(x, points2[0]) * 2.0;
	coeffs2[s2 - 1] = 1.0 / 6.0 * step2 * K(x, points2[s2 - 1]) * 2.0;
	for (i = 0; i < s1; i++)
	{
		result = result + (u1[i] * coeffs1[i]);
	}
	for (i = 0; i < s2; i++)
	{
		result = result - (u2[i] * coeffs2[i]);
	}
	delete[] coeffs1;
	delete[] coeffs2;
	delete[] points1;
	delete[] points2;
	return result * result;
}

double quadratura3(double* u1, int s1, double* u2, int s2, double aa, double bb, double ra, double rb)
{
	double arg1 = (aa + bb) / 2.0 - sqrt(3.0 / 5.0) * (bb - aa) / 2.0;
	double arg2 = (aa + bb) / 2.0;
	double arg3 = (aa + bb) / 2.0 + sqrt(3.0 / 5.0) * (bb - aa) / 2.0;
	double result = 0;
	result = (bb - aa) / 2.0 * (5.0 / 9.0 * new_func(u1, s1, u2, s2, arg1, ra, rb) + 8.0 / 9.0 * new_func(u1, s1, u2, s2, arg2, ra, rb) + 5.0 / 9.0 * new_func(u1, s1, u2, s2, arg3, ra, rb));
	return result;
}

double adaptive_calc(double* u1, int s1, double* u2, int s2, double aa, double bb)
{
	double h = (bb - aa) / 2;
	int i, j;
	double Ih = 0;
	double Ih2 = 0;
	double p = 1;
	double ro = 0;
	double res = 0;
	double cur_pos = aa;
	do
	{
		h = h * 2;
		do {
			Ih = quadratura3(u1, s1, u2, s2, cur_pos, cur_pos + h, aa, bb);
			Ih2 = quadratura3(u1, s1, u2, s2, cur_pos, cur_pos + h / 2, aa, bb) + quadratura3(u1, s1, u2, s2, cur_pos + h / 2, cur_pos + h, aa, bb);
			ro = (Ih2 - Ih) / (pow(2, p) - 1);
			h = h / 2;
		} while (abs(ro) > eps * h / (bb - aa));
		h = h * 2;
		cur_pos = cur_pos + h;
		res = res + Ih2;
	} while (abs(cur_pos - bb) > eps);
	return sqrt(res);
}

double C_norm(double* u1, double s1, double aa, double bb)
{
	int i; int j;
	double max = 0;
	double* points1 = new double[s1];
	double hold;
	double step = (bb - aa) / (s1 - 1.0);
	for (i = 0; i < s1; i++)
	{
		points1[i] = i * step + aa;
	}
	for (i = 0; i < s1; i++)
	{
		hold = abs(ideal_solution(points1[i]) - u1[i]);
		if (hold > max)
			max = hold;
	}
	return max;
}

double L2_norm_func(double* u2, int s2, double x, double aa, double bb)
{
	int i, j;
	double result = 0;
	double* points2 = new double[s2];
	double step2 = (bb - aa) / (s2 - 1);
	double* coeffs2 = new double[s2];
	for (i = 0; i < s2; i++)
	{
		points2[i] = i * step2 + aa;
	}
	for (i = 0; i < s2; i++)
	{
		if (i % 2 == 0)
		{
			coeffs2[i] = 2.0 / 6.0 * step2 * K(x, points2[i]) * 2.0;
		}
		else
		{
			coeffs2[i] = 4.0 / 6.0 * step2 * K(x, points2[i]) * 2.0;
		}
	}
	coeffs2[0] = 1.0 / 6.0 * step2 * K(x, points2[0]) * 2.0;
	coeffs2[s2 - 1] = 1.0 / 6.0 * step2 * K(x, points2[s2 - 1]) * 2.0;
	result = ideal_solution(x);
	for (i = 0; i < s2; i++)
	{
		result = result - (u2[i] * coeffs2[i]);
	}
	result = result - func(x);
	delete[] coeffs2;
	delete[] points2;
	return result * result;
}

double quadratura3_extra(double* u2, int s2, double aa, double bb, double ra, double rb)
{
	double arg1 = (aa + bb) / 2.0 - sqrt(3.0 / 5.0) * (bb - aa) / 2.0;
	double arg2 = (aa + bb) / 2.0;
	double arg3 = (aa + bb) / 2.0 + sqrt(3.0 / 5.0) * (bb - aa) / 2.0;
	double result = 0;
	result = (bb - aa) / 2.0 * (5.0 / 9.0 * L2_norm_func(u2, s2, arg1, ra, rb) + 8.0 / 9.0 * L2_norm_func(u2, s2, arg2, ra, rb) + 5.0 / 9.0 * L2_norm_func(u2, s2, arg3, ra, rb));
	return result;
}

double adaptive_calc_norm(double* u2, int s2, double aa, double bb)
{
	double h = (bb - aa) / 2;
	int i, j;
	double Ih = 0;
	double Ih2 = 0;
	double p = 1;
	double ro = 0;
	double res = 0;
	double cur_pos = aa;
	do
	{
		h = h * 2;
		do {
			Ih = quadratura3_extra(u2, s2, cur_pos, cur_pos + h, aa, bb);
			Ih2 = quadratura3_extra(u2, s2, cur_pos, cur_pos + h / 2, aa, bb) + quadratura3_extra(u2, s2, cur_pos + h / 2, cur_pos + h, aa, bb);
			ro = (Ih2 - Ih) / (pow(2, p) - 1);
			h = h / 2;
		} while (abs(ro) > eps * h / (bb - aa));
		h = h * 2;
		cur_pos = cur_pos + h;
		res = res + Ih2;
	} while (abs(cur_pos - bb) > eps);
	return sqrt(res);
}

void show(double** a, int m)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			cout << fixed << setprecision(7) << a[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;
}

int main()
{
	int i; int j;
	double adapt = 0;
	int point_num = 2;
	double aa = -1.0;
	double bb = 1.0;
	int vec_size1 = point_num * 2 - 1;
	double* vec1 = new double[vec_size1];
	double* fv = new double[vec_size1];
	double** data = new double* [vec_size1];
	for (i = 0; i < vec_size1; i++)
	{
		data[i] = new double[vec_size1];
	}
	system_build(data, fv, aa, bb, vec_size1);
	system_solve(data, fv, vec_size1);
	eq_1(vec1, fv, vec_size1);
	cout << "Points " << vec_size1 << endl;
	show(vec1, vec_size1);
	delete[] fv;
	for (i = 0; i < vec_size1; i++)
	{
		delete[] data[i];
	}
	delete[] data;
	point_num = point_num + 2;
	int vec_size2 = point_num * 2 - 1;
	double* vec2 = new double[vec_size2];
	fv = new double[vec_size2];
	data = new double* [vec_size2];
	for (i = 0; i < vec_size2; i++)
	{
		data[i] = new double[vec_size2];
	}
	system_build(data, fv, aa, bb, vec_size2);
	system_solve(data, fv, vec_size2);
	eq_1(vec2, fv, vec_size2);
	cout << "Points " << vec_size2 << endl;
	show(vec2, vec_size2);
	adapt = adaptive_calc(vec2, vec_size2, vec1, vec_size1, aa, bb);
	if (adapt > eps)
	{
		while (adapt > eps)
		{
			point_num = point_num + 2;
			vec_size1 = vec_size2;
			delete[] vec1;
			vec1 = new double[vec_size1];
			eq_1(vec1, vec2, vec_size1);
			delete[] vec2;
			delete[] fv;
			for (j = 0; j < vec_size2; j++)
			{
				delete[] data[j];
			}
			delete[] data;
			vec_size2 = point_num * 2 - 1;
			vec2 = new double[vec_size2];
			fv = new double[vec_size2];
			data = new double* [vec_size2];
			for (j = 0; j < vec_size2; j++)
			{
				data[j] = new double[vec_size2];
			}
			system_build(data, fv, aa, bb, vec_size2);
			system_solve(data, fv, vec_size2);
			eq_1(vec2, fv, vec_size2);
			adapt = adaptive_calc(vec2, vec_size2, vec1, vec_size1, aa, bb);
			cout << fixed << setprecision(20) << "||u(n+1)-u(n)|| = " << adapt << endl;
			cout << "Points: " << vec_size2 << endl;
			show(vec2, vec_size2);
		}
	}
	cout << endl;
	cout << "C_norm: " << fixed << setprecision(20) << C_norm(vec2, vec_size2, aa, bb) << endl;
	cout << "L2_norm: " << adaptive_calc_norm(vec2, vec_size2, aa, bb) << endl;
	return 0;
}

