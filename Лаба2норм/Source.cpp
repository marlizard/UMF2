#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>

using namespace std;

typedef double real;

//Порядок функции
#define Function 2

class Method {
public:
	//Число максимальных итераций
	int maxiter = 10000;
	//Значение точности
	real e = 1.0e-14;
	//Значение параметра релаксации
	real w = 1.00;
	//Маркер метода
	int method;

	//Количество конечных элементов
	int n;	//кол-во КЭ

	//Вектор координат
	vector<real> cord;

	//Нижняя диагональ глобальной матрицы
	vector<real> b;

	//Главная диагональ глобальной матрицы
	vector<real> c;

	//Верхняя диагональ глобально матрицы
	vector<real> d;

	//Вектор правой части
	vector<real> f;

	//Вектор весов
	vector<real> q;

	//Вектор весов на предыдущей итерации
	vector<real> q1;

	//Вспомогательный вектор
	vector<real> temp;

	//Вектор значений гамма
	vector<real> Gamma;

	//Вектор значений бета
	vector<real> Beta;

	//Значения искомой функции
	real u(real x, int function);

	//Функция лямбда, зависящая от производной u
	real lambda(real t_aproks, real x);

	//Производная функции лямбда
	real der_lambda(real t_aproks, real x);

	//Функция правой части
	real func(real x, real gamma);

	//Производная искомой функции
	real der_u(real x);

	//Расчет нормаы вектора
	real norm(vector<real> t);

	//Умножение глобальной матрицы на вектор весов
	void multiply();

	//Расчет невязки
	real discrepancy(int k);

	//Генератор сетки
	void grid_generator();

	//Функция ввода
	void input();

	//Обнуление глобальной матрицы и вектора правой части
	void clear_matrix();

	//Добавление коэффициентов для метода Ньютона
	real newton_add(real Xk, real Xk1, real h, real q1, real q2, real arg, int num);

	//Добавление локальной матрицы конечного элемента в глобальную
	void local_matrix(int j);

	//Создание глобальной матрицы
	void global_matrix();

	//Подсчет тетты для вторых краевых
	real tetta_count(int num);

	//Подсчет u бета для третьих краевых
	real u_beta(int num);

	//Учет первых краевых
	void dirichlet_condition(int i);

	//Учет вторых краевых 
	void neyman_condition(int i);

	//Учет третьих краевых
	void newton_condition(int i);

	//Учет краевых условий
	void boundary_conditions();

	//Расчет нового значения вектора весов
	void sweep_method();

	//Преобразование вектора весов с учетом параметра релаксации
	void relaxation();
	
	//Функция вывода
	void output(int iter);

	//Сохранение текущего вектора весов как предыдущего
	void save_previous();
	
	//Реализация метода
	void universal_method();
};

real Method::u(real x, int function)
{
	switch (function)
	{
	case -1: return 1;
	case 0: return 4.675;
	case 1: return x;
	case 2: return x * x;
	case 3: return pow(x, 3);
	case 4: return sin(x);
	}
}

real Method::der_lambda(real t_aproks, real x)
{
	return 1;
}

real Method::lambda(real t_aproks, real x)
{
	switch (Function)
	{
	case 0: return t_aproks + 1;
	case 1: return t_aproks + 1;
	case 2: return t_aproks + 1;
	case 3: return t_aproks + 1;
	case 4: return t_aproks + 1;
	}
}

real Method::func(real x, real gamma)
{
	switch (Function)
	{
	case 0: return gamma * u(x, Function);
	case 1: return gamma * u(x, Function);
	case 2: return -8 * x - 2 + gamma * u(x, Function);
	case 3: return -36 * pow(x, 3) - 6 * x + gamma * u(x, Function);
	case 4: return sin(x) - sin(2 * x) + gamma * u(x, Function);
	}
}

real Method::der_u(real x)
{
	switch (Function)
	{
	case 0: return 0;
	case 1: return 1;
	case 2: return 2 * x;
	case 3: return 3 * x * x;
	case 4: return sin(x);
	}
}


void Method::grid_generator() {
	ifstream input("input.txt");
	input >> n;
	cord.resize(n + 1, 0);
	real a, b, k, h;

	input >> a >> b >> k;

	if (k == 1) h = (b - a) / n;
	else h = (b - a)* (k - 1) / (pow(k, n) - 1);
	
	input.close();

	cord[0] = a;

	for (int i = 1; i < cord.size(); i++) {
		a += h;
		cord[i] = a;
		h *= k;
	}

	cout << "Please, choose the method:" << endl;
	cout << "1. Simple iteration method" << endl;
	cout << "2. Newton method" << endl;
	cin >> method;

	q.resize(n + 1, 0);
	q1.resize(n + 1, 0);

	if (method == 2) {
		for (int i = 0; i < q.size(); i++)
			q[i] = u(cord[i], Function - 1);
	}
	else {
		for (int i = 0; i < q.size(); i++)
			q[i] = 1;
	}
}

void Method::input()
{
	grid_generator();

	b.resize(n + 1, 0);
	c.resize(n + 1, 0);
	d.resize(n + 1, 0);
	f.resize(n + 1, 0);
	temp.resize(n + 1, 0);
	Gamma.resize(n, 0);
	Beta.resize(n, 0);

	for (int i = 0; i < n; i++)
		Gamma[i] = 1;

	for (int i = 0; i < n; i++)
		Beta[i] = 1;
}

real Method::norm(vector<real> t)
{
	real norma = 0;
	for (int i = 0; i <= n; i++)
		norma += pow(t[i], 2);
	return norma;
}

void Method::multiply()
{
	temp[0] = c[0] * q[0] + d[0] * q[1];
	temp[n] = b[n] * q[n - 1] + c[n] * q[n];

	for (int i = 1; i < n; i++) {
		temp[i] = b[i] * q[i - 1];
		temp[i] += c[i] * q[i];
		temp[i] += d[i] * q[i + 1];
	}
}

real Method::discrepancy(int k)
{
	multiply();
	for (int i = 0; i <= n; i++)
		temp[i] = f[i] - temp[i];

	return sqrt(norm(temp) / norm(f));
}

real Method::newton_add(real Xk, real Xk1, real h, real q1, real q2, real der_f, int num)
{
	real dLdQ1_left, dLdQ1_right, dLdQ2_left, dLdQ2_right;
	dLdQ1_left = -der_lambda(der_f, Xk) / h;
	dLdQ1_right = -der_lambda(der_f, Xk1) / h;
	dLdQ2_left = der_lambda(der_f, Xk) / h;
	dLdQ2_right = der_lambda(der_f, Xk1) / h;

	real dA11dQ1, dA11dQ2, dA12dQ1, dA12dQ2, dA21dQ1, dA21dQ2, dA22dQ1, dA22dQ2;

	dA11dQ1 = dA22dQ1 = (dLdQ1_left + dLdQ1_right) / (2 * h);
	dA11dQ2 = dA22dQ2 = (dLdQ2_left + dLdQ2_right) / (2 * h);
	dA12dQ1 = dA21dQ1 = -(dLdQ1_left + dLdQ1_right) / (2 * h);
	dA12dQ2 = dA21dQ2 = -(dLdQ2_left + dLdQ2_right) / (2 * h);

	switch (num) {
	case 1:
		return (dA11dQ1 * q1 + dA12dQ1 * q2);
	case 2:
		return (dA11dQ2 * q1 + dA12dQ2 * q2);
	case 3:
		return (dA21dQ1 * q1 + dA22dQ1 * q2);
	case 4:
		return (dA21dQ2 * q1 + dA22dQ2 * q2);
	case 5:
		return (q1 * (dA11dQ1 * q1 + dA11dQ2 * q2) + q2 * (dA12dQ1 * q1 + dA12dQ2 * q2));
	case 6:
		return (q1 * (dA21dQ1 * q1 + dA21dQ2 * q2) + q2 * (dA22dQ1 * q1 + dA22dQ2 * q2));
	}
}

void Method::local_matrix(int j)
{
	real Xk, Xk1, h, gamma, der_f, q1, q2;

	Xk = cord[j];
	Xk1 = cord[j + 1];
	h = Xk1 - Xk;
	gamma = Gamma[j];
	q1 = q[j];
	q2 = q[j + 1];
	der_f = (q2 - q1) / h;

	if (method == 2) {
		//считаем матрицу жесткости
		c[j] += (lambda(der_f, Xk) + lambda(der_f, Xk1)) / (2 * h) + newton_add(Xk, Xk1, h, q1, q2, der_f, 1);
		d[j] -= (lambda(der_f, Xk) + lambda(der_f, Xk1)) / (2 * h) - newton_add(Xk, Xk1, h, q1, q2, der_f, 2);
		b[j + 1] -= (lambda(der_f, Xk) + lambda(der_f, Xk1)) / (2 * h) - newton_add(Xk, Xk1, h, q1, q2, der_f, 3);
		c[j + 1] += (lambda(der_f, Xk) + lambda(der_f, Xk1)) / (2 * h) + newton_add(Xk, Xk1, h, q1, q2, der_f, 4);

		//считаем вектор правой части
		f[j] += h * (2 * func(Xk, gamma) + func(Xk1, gamma)) / 6 + newton_add(Xk, Xk1, h, q1, q2, der_f, 5);
		f[j + 1] += h * (func(Xk, gamma) + 2 * func(Xk1, gamma)) / 6 + newton_add(Xk, Xk1, h, q1, q2, der_f, 6);
	}
	else {
		c[j] += (lambda(der_f, Xk) + lambda(der_f, Xk1)) / (2 * h);
		d[j] -= (lambda(der_f, Xk) + lambda(der_f, Xk1)) / (2 * h);
		b[j + 1] -= (lambda(der_f, Xk) + lambda(der_f, Xk1)) / (2 * h);
		c[j + 1] += (lambda(der_f, Xk) + lambda(der_f, Xk1)) / (2 * h);

		f[j] += h * (2 * func(Xk, gamma) + func(Xk1, gamma)) / 6;
		f[j + 1] += h * (func(Xk, gamma) + 2 * func(Xk1, gamma)) / 6;
	}

	//считаем матрицу массы
	c[j] += gamma * h / 3;
	d[j] += gamma * h / 6;
	b[j + 1] += gamma * h / 6;
	c[j + 1] += gamma * h / 3;
}

void Method::relaxation()
{
	for (int i = 0; i <= n; i++)
		q[i] = w * q[i] + (1 - w) * q1[i];
}

void Method::clear_matrix()
{
	for (int i = 0; i <= n; i++)
		f[i] = b[i] = c[i] = d[i] = 0;
}

real Method::tetta_count(int num)
{
	real q1, q2, h, x, ret; //возвращаемое значение
	real dx; //значение производной

	if (!num) {
		q1 = q[0];
		q2 = q[1];
		h = cord[1] - cord[0];
		x = cord[0];
		dx = der_u(cord[0]);
	}
	else {
		q1 = q[n - 1];
		q2 = q[n];
		h = cord[n] - cord[n - 1];
		x = cord[n];
		dx = der_u(cord[n]);
	}

	ret = dx * lambda((q2 - q1) / h, x); //ret = lambda*(dU/dx)
	if (!num)//если левая граница
		return -ret;
	else return ret;
}

real Method::u_beta(int num)
{
	real q1, q2, h, x, ret, beta; //возвращаемое значение
	real dx; //значение производной

	//u_beta = (lambda* dU/dn + U|x=x0)/betta
	if (!num) {
		q1 = q[0];
		q2 = q[1];
		beta = Beta[0];
		h = cord[1] - cord[0];
		x = cord[0];
		dx = -der_u(x);
	}
	else {
		q1 = q[n - 1];
		q2 = q[n];
		beta = Beta[n - 1];
		h = cord[n] - cord[n - 1];
		x = cord[n];
		dx = der_u(x);
	}

	ret = 1 / beta * (lambda((q2 - q1) / h, x) * dx + beta * u(x, Function));
	return ret;
}

void Method::dirichlet_condition(int i)
{
	if (!i) {//если левая граница отрезка
		c[0] = 1;
		d[0] = 0;
		f[0] = u(cord[0], Function);
	}
	else {//если правая граница отрезка
		c[n] = 1;
		b[n] = 0;
		f[n] = u(cord[n], Function);
	}
}

void Method::neyman_condition(int i)
{
	switch (i) {
	case 0://если левая граница 
		f[0] += tetta_count(0);
		break;
	case 1://если правая граница
		f[n] += tetta_count(1);
		break;
	}
}

void Method::newton_condition(int i)
{
	if (!i) {//если левая граница
		c[0] += Beta[0];
		f[0] += Beta[0] * u_beta(i);
	}
	else {//если правая граница
		c[n] += Beta[n - 1];
		f[n] += Beta[n - 1] * u_beta(i);
	}
}

void Method::boundary_conditions()
{
	int type;
	ifstream cond("boundary_conditions.txt");
	
	for (int i = 0; i < 2; i++) {
		
		cond >> type;
		switch (type) {
		case 1:
			dirichlet_condition(i);
			break;
		case 2:
			neyman_condition(i);
			break;
		case 3:
			newton_condition(i);
			break;
		}
	}
}

void Method::output(int iter)
{
	ofstream output("output.txt");
	output << "Iteration number = " << iter << endl;

	cout << "Discrepancy value = " << scientific << setprecision(14) << discrepancy(iter) << endl << endl;

	for (int i = 0; i <= n; i++) {
		output << fixed;
		output << setprecision(14) << q[i] << endl;
		cout << fixed;
		cout << setprecision(14) << q[i] << " " << u(cord[i], Function) << " " << scientific << abs(q[i] - u(cord[i], Function)) << endl;
	}
}

void Method::save_previous()
{
	for (int j = 0; j <= n; j++)
		q1[j] = q[j];
}

void Method::global_matrix()
{
	clear_matrix();
	for (int j = 0; j < n; j++) 
		local_matrix(j);	
	boundary_conditions();
}

void Method::sweep_method()
{
	int i;
	vector<real> alpha(n + 1, 0), beta(n + 1, 0); 

	alpha[0] = (-1) * d[0] / c[0];
	beta[0] = f[0] / c[0];
	for (i = 1; i <= n; i++) {
		alpha[i] = (-1) * d[i] / (c[i] + b[i] * alpha[i - 1]);
		beta[i] = (f[i] - b[i] * beta[i - 1]) / (c[i] + b[i] * alpha[i - 1]);
	}

	q[n] = beta[n];
	for (i = n - 1; i > -1; i--) {
		q[i] = alpha[i] * q[i + 1] + beta[i];
	}
}


void Method::universal_method()
{
	int i;
	save_previous();
	global_matrix();
	for (i = 0; i < maxiter && discrepancy(i) >= e; i++) {
		save_previous();
		sweep_method();
		relaxation();
		global_matrix();
	}
	output(i);
}

int main()
{
	Method m;
	m.input();
	m.universal_method();
	return 0;
}
