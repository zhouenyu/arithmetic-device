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
}; //���帴���Ľṹ��,real��ʾʵ������,imag��ʾ��������;
void c_plus(complex a, complex b, complex *c)
{
	c->real = a.real + b.real;
	c->imag = a.imag + b.imag;
}//���ڼ��㸴���ļӷ�
void Wn_i(int n, int i, complex *Wn, char flag)
{
	Wn->real = cos(2 * PI*i / n);
	if (flag == 1)
		Wn->imag = -sin(2 * PI*i / n);
	else if (flag == 0)
		Wn->imag = -sin(2 * PI*i / n);
}//���ڼ�����ٸ���Ҷ�任�е�wn
void c_sub(complex a, complex b, complex *c)
{
	c->real = a.real - b.real;
	c->imag = a.imag - b.imag;
}//���ڼ��㸴���ļ���
void c_mul(complex a, complex b, complex *c)
{
	c->real = a.real * b.real - a.imag * b.imag;
	c->imag = a.real * b.imag + a.imag * b.real;
}//���ڼ��㸴���ĳ˷�
void conjugate_complex(int n, complex in[], complex out[])
{
	int i = 0;
	for (i = 0; i < n; i++)
	{
		out[i].imag = -in[i].imag;
		out[i].real = in[i].real;
	}
}//���������Ĺ���
void fft(int N, complex f[])//����fft�㷨,���������ɢ����Ҷ�任
{
	complex t, wn;//�м����  
	int i, j, k, m, n, l, r, M;
	int la, lb, lc;
	/*----����ֽ�ļ���M=log2(N)----*/
	for (i = N, M = 1; (i = i / 2) != 1; M++);
	/*----���յ�λ����������ԭ�ź�----*/
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
		la = pow(2, m); //la=2^m�����m��ÿ�����������ڵ���       
		lb = la / 2;    //lb�����m��ÿ�������������ε�Ԫ��  
					 //ͬʱ��Ҳ��ʾÿ�����ε�Ԫ���½ڵ�֮��ľ���  
			/*----��������----*/
		for (l = 1; l <= lb; l++)
		{
			r = (l - 1)*pow(2, M - m);
			for (n = l - 1; n < N - 1; n = n + la) //����ÿ�����飬��������ΪN/la  
			{
				lc = n + lb;  //n,lc�ֱ����һ�����ε�Ԫ���ϡ��½ڵ���       
				Wn_i(N, r, &wn, 1);//wn=Wnr  
				c_mul(f[lc], wn, &t);//t = f[lc] * wn��������  
				c_sub(f[n], t, &(f[lc]));//f[lc] = f[n] - f[lc] * Wnr  
				c_plus(f[n], t, &(f[n]));//f[n] = f[n] + f[lc] * Wnr  
			}
		}
	}
}
void ifft(int N, complex f[])//���ڼ��㸵��Ҷ��任
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
bool ISNUM(string num)//�ж�����������Ƿ�Ϸ�
{
	unsigned int i = 0;
	if (num.length() == 0)//�ж��Ƿ�Ϊ��
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
string STANDARDIZE(string l)//��l�ַ�����׼��,ȥ��֮ǰ��������,�Ͷ��ڵ�0
{
	int length, i = 0, first = 0;
	length = l.length();
	if (l[i] == '+' || l[i] == '-')//ȥ��������
	{
		first++;
		i++;
	}
	for (; i < length; i++)
	{
		if (l[i] == '0')//ȥ�����ڵ�0
			first++;
		else
			break;
	}
	l = l.substr(first, l.length());
	if (l.size() == 0)//���������0,��������һ��0
		l = l + '0';
	return l;
}
bool ISSYMBOL(string symbol)//�ж�������ַ��Ƿ�ϸ�
{
	if (symbol.length() != 1)//�ж��Ƿ�Ϊ��
		return false;
	else
		if (symbol[0] == '+' || symbol[0] == '-' || symbol[0] == '*' || symbol[0] == '/')
			return true;
		else
			return false;
}
string ADD(string num_1, string num_2)//�ӷ�����,ԭ��Ϊ��ʽ�ӷ�,����һ��string���͵��ַ�����Ϊ���
{
	string num_return;
	unsigned int i, l1, l2, l_max, d = 0, l_min;
	l1 = num_1.length();
	l2 = num_2.length();
	l_max = max(l1, l2);
	l_min = min(l1, l2);
	num_return.resize(l_max);
	for (i = 1; i <= l_min; i++)//��ÿһλ���мӷ�����,d��ʾ��λ
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
string SUB(string num_1, string num_2)//��������,ԭ��Ϊ��ʽ����,���Ϊstring�����ַ���,���������֮��
{
	unsigned int  l1, l2, i;
	int d = 0, c = 0;
	string num_return;
	l1 = num_1.length();
	l2 = num_2.length();
	if (l1 < l2)//�жϽ���������������Ƿ�Ҫ�����
		return '-' + SUB(num_2, num_1);
	if (l1 == l2)
		if (num_1 < num_2)
			return '-' + SUB(num_2, num_1);
		else
			if (num_1 == num_2)
				return "0";
	num_return.resize(l1);//num_return���ڴ�Ž��
	for (i = 1; i <= l2; i++)//ÿһλ���
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
	return STANDARDIZE(num_return);//�������׼���󷵻�
}
string MUL(string num_1, string num_2)
//�˷�ԭ��Ϊ����num_1��num_2�ĸ���Ҷ�仯,Ȼ����˵õ�num_return�ĸ���Ҷ�任,Ȼ���ٽ��и���Ҷ��任�õ�num_return
{
	complex *a, *b, *c;
	int l, i, l_return;
	l = num_1.length() + num_2.length();
	l = (int)pow(2, ceil(log10(l) / log10(2)));//���㸵��Ҷ�任�����ά��,Ҫ��Ϊ2^nά
	a = new complex[l];
	b = new complex[l];
	c = new complex[l];
	for (i = 1; i <= num_1.length(); i++)//������num_1,������������a��ʵ��,�鲿Ϊ0
	{
		a[i - 1].real = (double)num_1[num_1.length() - i] - '0';
		a[i - 1].imag = 0;
	}
	for (; i <= l; i++)//ʣ������鲿ʵ����Ϊ0
	{
		a[i - 1].real = 0;
		a[i - 1].imag = 0;
	}
	for (i = 1; i <= num_2.length(); i++)//������num_2,������������b��ʵ��,�鲿Ϊ0
	{
		b[i - 1].real = (double)num_2[num_2.length() - i] - '0';
		b[i - 1].imag = 0;
	}
	for (; i <= l; i++)//ʣ������鲿ʵ����Ϊ0
	{
		b[i - 1].real = 0;
		b[i - 1].imag = 0;
	}
	fft(l, a);//��a,b����ֱ���и���Ҷ�任
	fft(l, b);
	for (i = 0; i < l; i++)//������a,b��˵õ�c
		c_mul(a[i], b[i], &c[i]);
	ifft(l, c);//������c���и���Ҷ��任
	string num_return(l, '0');
	l_return = l;
	int m = 0, n = 0, d = 0, j = 0, k = 0;
	for (i = 0; i < l; i++)//������c�任Ϊnum_return���
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
string DIV(string num_1, string num_2, string &num_remainder)//������ʽ�����ļ��㷽��
{
	int l, l1, l2;
	num_1 = STANDARDIZE(num_1);//��׼��
	num_2 = STANDARDIZE(num_2);
	if (num_2[0] == '0')//�жϳ����Ƿ�Ϊ��
	{
		cout << "Error input" << endl;
		exit(1);
	}
	l1 = num_1.size();
	l2 = num_2.size();
	l = num_1.size() - num_2.size();
	string s, num_quotient(l + 1, '0'), k;
	if (l < 0)//�жϳ����ͱ������Ĵ�С
	{
		num_remainder = num_2;
		return num_quotient = '0';
	}
	s = num_1.substr(0, l2);
	for (int i = 0; i <= l; i++)//���г�������
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
	num_quotient = STANDARDIZE(num_quotient);//�������
	num_remainder = STANDARDIZE(num_remainder);
	return num_quotient;//�����
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
		if (ISNUM(num_1) + ISNUM(num_2) != 2)//�ж����������Ƿ�Ϸ�
		{
			cout << "Error input" << endl;
			return 1;
		}
		if (!ISSYMBOL(symbol))//�ж�����ַ��Ƿ�Ϸ�
		{
			cout << "Error input" << endl;
			return 1;
		}
		if (symbol[0] == '+')//�ӷ�����
		{
			if (num_1[0] != '-'&&num_2[0] != '-')//�����������ֵķ�������
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
		if (symbol[0] == '-')//��������
		{
			if (num_1[0] != '-'&&num_2[0] != '-')//�������
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
		if (symbol[0] == '*')//�˷�����
		{
			if (num_1[0] != '-'&&num_2[0] != '-')//�������
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
		if (symbol[0] == '/')//��������
		{
			num_return = DIV(num_1, num_2, num_remainder_return);
			cout << num_return << ' ' << num_remainder_return << endl;//����������
		}
		else
			cout << num_return << endl;//�Ӽ���������
	}
	return 0;
}





