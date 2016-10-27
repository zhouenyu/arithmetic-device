#define PI 3.14159265358979323846264338327950288419716939937510582097494
#include<stdio.h>
#include<iostream>
#include<string>
#include<algorithm>
#include<math.h>
#include <cstdlib>
using namespace std;
struct complex
{
	double real, imag;
}; //定义复数的结构体,real表示实数部分,imag表示复数部分;
void c_plus(complex a, complex b, complex *c)
{
	c->real = a.real + b.real;
	c->imag = a.imag + b.imag;
}//用于计算复数的加法
void Wn_i(int n, int i, complex *Wn, char flag)
{
	Wn->real = cos(2 * PI*i / n);
	if (flag == 1)
		Wn->imag = -sin(2 * PI*i / n);
	else if (flag == 0)
		Wn->imag = -sin(2 * PI*i / n);
}//用于计算快速傅里叶变换中的wn
void c_sub(complex a, complex b, complex *c)
{
	c->real = a.real - b.real;
	c->imag = a.imag - b.imag;
}//用于计算复数的减法
void c_mul(complex a, complex b, complex *c)
{
	c->real = a.real * b.real - a.imag * b.imag;
	c->imag = a.real * b.imag + a.imag * b.real;
}//用于计算复数的乘法
void conjugate_complex(int n, complex in[], complex out[])
{
	int i = 0;
	for (i = 0; i < n; i++)
	{
		out[i].imag = -in[i].imag;
		out[i].real = in[i].real;
	}
}//用于求复数的共轭
void fft(int N, complex f[])//快速fft算法,用于求解离散傅里叶变换
{
	complex t, wn;//中间变量  
	int i, j, k, m, n, l, r, M;
	int la, lb, lc;
	/*----计算分解的级数M=log2(N)----*/
	for (i = N, M = 1; (i = i / 2) != 1; M++);
	/*----按照倒位序重新排列原信号----*/
	for (i = 1, j = N / 2; i <= N - 2; i++)
	{
		if (i < j)
		{
			t = f[j];
			f[j] = f[i];
			f[i] = t;
		}
		k = N / 2;
		while (k <= j)
		{
			j = j - k;
			k = k / 2;
		}
		j = j + k;
	}
	for (m = 1; m <= M; m++)
	{
		la = pow(2, m); //la=2^m代表第m级每个分组所含节点数       
		lb = la / 2;    //lb代表第m级每个分组所含碟形单元数  
					 //同时它也表示每个碟形单元上下节点之间的距离  
			/*----碟形运算----*/
		for (l = 1; l <= lb; l++)
		{
			r = (l - 1)*pow(2, M - m);
			for (n = l - 1; n < N - 1; n = n + la) //遍历每个分组，分组总数为N/la  
			{
				lc = n + lb;  //n,lc分别代表一个碟形单元的上、下节点编号       
				Wn_i(N, r, &wn, 1);//wn=Wnr  
				c_mul(f[lc], wn, &t);//t = f[lc] * wn复数运算  
				c_sub(f[n], t, &(f[lc]));//f[lc] = f[n] - f[lc] * Wnr  
				c_plus(f[n], t, &(f[n]));//f[n] = f[n] + f[lc] * Wnr  
			}
		}
	}
}
void ifft(int N, complex f[])//用于计算傅里叶逆变换
{
	int i = 0;
	conjugate_complex(N, f, f);
	fft(N, f);
	conjugate_complex(N, f, f);
	for (i = 0; i < N; i++)
	{
		f[i].imag = (f[i].imag) / N;
		f[i].real = (f[i].real) / N;
	}
}
bool ISNUM(string num)//判断输入的数字是否合法
{
	unsigned int i = 0;
	if (num.length() == 0)//判断是否为空
		return false;
	if (num[0] == '+' || num[0] == '-')
		i = 1;
	for (; i < num.length(); i++)
	{
		if (num[i]>'9' || num[i] < '0')
			return false;
	}
	return true;
}
string STANDARDIZE(string l)//将l字符串标准化,去除之前的正负号,和多于的0
{
	int length, i = 0, first = 0;
	length = l.length();
	if (l[i] == '+' || l[i] == '-')//去除正负号
	{
		first++;
		i++;
	}
	for (; i < length; i++)
	{
		if (l[i] == '0')//去除多于的0
			first++;
		else
			break;
	}
	l = l.substr(first, l.length());
	if (l.size() == 0)//如果输入多个0,则最后输出一个0
		l = l + '0';
	return l;
}
bool ISSYMBOL(string symbol)//判断输入的字符是否合格
{
	if (symbol.length() != 1)//判断是否为空
		return false;
	else
		if (symbol[0] == '+' || symbol[0] == '-' || symbol[0] == '*' || symbol[0] == '/')
			return true;
		else
			return false;
}
string ADD(string num_1, string num_2)//加法运算,原理为竖式加法,返回一个string类型的字符串作为结果
{
	string num_return;
	unsigned int i, l1, l2, l_max, d = 0, l_min;
	l1 = num_1.length();
	l2 = num_2.length();
	l_max = max(l1, l2);
	l_min = min(l1, l2);
	num_return.resize(l_max);
	for (i = 1; i <= l_min; i++)//对每一位进行加法计算,d表示进位
	{
		num_return[l_max - i] = (num_1[l1 - i] - '0' + num_2[l2 - i] - '0' + d) % 10 + '0';
		d = (num_1[l1 - i] - '0' + num_2[l2 - i] - '0' + d) / 10;
	}
	for (; i <= l_max; i++)
	{
		if (l1 > l2)
		{
			num_return[l_max - i] = (num_1[l1 - i] - '0' + d) % 10 + '0';
			d = (num_1[l1 - i] - '0' + d) / 10;
		}
		else
		{
			num_return[l_max - i] = (num_2[l2 - i] - '0' + d) % 10 + '0';
			d = (num_2[l2 - i] - '0' + d) / 10;
		}
	}
	if (d == 1)
		num_return = '1' + num_return;
	return num_return;
}
string SUB(string num_1, string num_2)//减法运算,原理为竖式减法,输出为string类型字符串,输出有正负之分
{
	unsigned int  l1, l2, i;
	int d = 0, c = 0;
	string num_return;
	l1 = num_1.length();
	l2 = num_2.length();
	if (l1 < l2)//判断结果的正负来决定是否要输出零
		return '-' + SUB(num_2, num_1);
	if (l1 == l2)
		if (num_1 < num_2)
			return '-' + SUB(num_2, num_1);
		else
			if (num_1 == num_2)
				return "0";
	num_return.resize(l1);//num_return用于存放结果
	for (i = 1; i <= l2; i++)//每一位相减
	{
		c = num_1[l1 - i] - num_2[l2 - i] + d;
		if (c < 0)
		{
			d = -1;
			num_return[l1 - i] = c + 10 + '0';
		}
		else
		{
			num_return[l1 - i] = c + '0';
			d = 0;
		}
	}
	for (; i <= l1; i++)
	{
		c = num_1[l1 - i] - '0' + d;
		if (c < 0)
		{
			d = -1;
			num_return[l1 - i] = c + 10 + '0';
		}
		else
		{
			num_return[l1 - i] = c + '0';
			d = 0;
		}
	}
	return STANDARDIZE(num_return);//将结果标准化后返回
}
string MUL(string num_1, string num_2)
//乘法原理为计算num_1和num_2的傅里叶变化,然后相乘得到num_return的傅里叶变换,然后再进行傅里叶逆变换得到num_return
{
	complex *a, *b, *c;
	int l, i, l_return;
	l = num_1.length() + num_2.length();
	l = (int)pow(2, ceil(log10(l) / log10(2)));//计算傅里叶变换数组的维度,要求为2^n维
	a = new complex[l];
	b = new complex[l];
	c = new complex[l];
	for (i = 1; i <= num_1.length(); i++)//将数组num_1,付给复数数组a的实部,虚部为0
	{
		a[i - 1].real = (double)num_1[num_1.length() - i] - '0';
		a[i - 1].imag = 0;
	}
	for (; i <= l; i++)//剩余项的虚部实部都为0
	{
		a[i - 1].real = 0;
		a[i - 1].imag = 0;
	}
	for (i = 1; i <= num_2.length(); i++)//将数组num_2,付给复数数组b的实部,虚部为0
	{
		b[i - 1].real = (double)num_2[num_2.length() - i] - '0';
		b[i - 1].imag = 0;
	}
	for (; i <= l; i++)//剩余项的虚部实部都为0
	{
		b[i - 1].real = 0;
		b[i - 1].imag = 0;
	}
	fft(l, a);//对a,b数组分别进行傅里叶变换
	fft(l, b);
	for (i = 0; i < l; i++)//将向量a,b相乘得到c
		c_mul(a[i], b[i], &c[i]);
	ifft(l, c);//将向量c进行傅里叶逆变换
	string num_return(l, '0');
	l_return = l;
	int m = 0, n = 0, d = 0, j = 0, k = 0;
	for (i = 0; i < l; i++)//将向量c变换为num_return输出
	{
		j = 0;
		m = (int)(c[i].real + 0.5) % 10;
		d = (int)(c[i].real + 0.5) / 10;
		do {
			if (i + j > l_return - 1)
			{
				num_return = num_return + "0000000000";
				l_return = l_return + 10;
			}
			k = num_return[i + j] - '0' + m;
			d = d + k / 10;
			num_return[i + j] = k % 10 + '0';
			m = d % 10;
			d = d / 10;
			j++;
		} while (m != 0);
	}
	num_return = string(num_return.rbegin(), num_return.rend());
	num_return = STANDARDIZE(num_return);
	return num_return;
}
string DIV(string num_1, string num_2, string &num_remainder)//采用竖式除法的计算方法
{
	int l, l1, l2;
	num_1 = STANDARDIZE(num_1);//标准化
	num_2 = STANDARDIZE(num_2);
	if (num_2[0] == '0')//判断除数是否为零
	{
		cout << "Error input" << endl;
		exit(1);
	}
	l1 = num_1.size();
	l2 = num_2.size();
	l = num_1.size() - num_2.size();
	string s, num_quotient(l + 1, '0'), k;
	if (l < 0)//判断除数和被除数的大小
	{
		num_remainder = num_2;
		return num_quotient = '0';
	}
	s = num_1.substr(0, l2);
	for (int i = 0; i <= l; i++)//进行除法运算
	{
		k = SUB(s, num_2);
		while (k[0] != '-')
		{
			s = k;
			num_quotient[i]++;
			k = SUB(k, num_2);
		}
		if (l2 + i>l1 - 1)
		{
			num_remainder = s;
			break;
		}
		s = s + num_1[l2 + i];
		s = STANDARDIZE(s);
	}
	num_quotient = STANDARDIZE(num_quotient);//输出余数
	num_remainder = STANDARDIZE(num_remainder);
	return num_quotient;//输出商
}
int main()
{
	string num_1, num_2, symbol, num_return, num_remainder_return;
	for (int i = 0;;)
	{
		getline(cin, num_1);
		if (num_1 == "")
			return 0;
		getline(cin, num_2);
		if (num_2 == "")
			return 0;
		getline(cin, symbol);
		if (symbol == "")
			return 0;
		if (ISNUM(num_1) + ISNUM(num_2) != 2)//判断输入数字是否合法
		{
			cout << "Error input" << endl;
			return 1;
		}
		if (!ISSYMBOL(symbol))//判断输出字符是否合法
		{
			cout << "Error input" << endl;
			return 1;
		}
		if (symbol[0] == '+')//加法运算
		{
			if (num_1[0] != '-'&&num_2[0] != '-')//处理输入数字的符号问题
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = ADD(num_1, num_2);
			}
			if (num_1[0] == '-'&&num_2[0] == '-')
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = ADD(num_1, num_2);
				num_return = '-' + num_return;
			}
			if (num_1[0] == '-'&&num_2[0] != '-')
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = SUB(num_2, num_1);
			}
			if (num_1[0] != '-'&&num_2[0] == '-')
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = SUB(num_1, num_2);
			}
		}
		if (symbol[0] == '-')//减法运算
		{
			if (num_1[0] != '-'&&num_2[0] != '-')//处理符号
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = SUB(num_1, num_2);
			}
			if (num_1[0] == '-'&&num_2[0] == '-')
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = SUB(num_2, num_1);
			}
			if (num_1[0] == '-'&&num_2[0] != '-')
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = ADD(num_2, num_1);
				num_return = '-' + num_return;
			}
			if (num_1[0] != '-'&&num_2[0] == '-')
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = ADD(num_1, num_2);
			}
		}
		if (symbol[0] == '*')//乘法运算
		{
			if (num_1[0] != '-'&&num_2[0] != '-')//处理符号
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = MUL(num_1, num_2);
			}
			if (num_1[0] == '-'&&num_2[0] == '-')
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = MUL(num_2, num_1);
			}
			if (num_1[0] == '-'&&num_2[0] != '-')
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = MUL(num_2, num_1);
				num_return = '-' + num_return;
			}
			if (num_1[0] != '-'&&num_2[0] == '-')
			{
				num_1 = STANDARDIZE(num_1);
				num_2 = STANDARDIZE(num_2);
				num_return = MUL(num_1, num_2);
				num_return = '-' + num_return;
			}

		}
		if (symbol[0] == '/')//除法运算
		{
			num_return = DIV(num_1, num_2, num_remainder_return);
			cout << num_return << ' ' << num_remainder_return << endl;//除法输出结果
		}
		else
			cout << num_return << endl;//加减乘输出结果
	}
	return 0;
}





