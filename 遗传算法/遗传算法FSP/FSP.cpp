#include "bits/stsc++.h"

using namespace std;

const int gongjian_num = 5;			//��������
const int machine_num =  4;			//��������	
const int maxn = 8;					//��Ⱥ���� 
const double pc = 0.6;				//������� 
const double pm = 0.1;				//������� 
const int iteration = 1000;
int t[gongjian_num][machine_num];
int best_time = 9999;
int best[gongjian_num];
int MAX(int a, int b)
{
	if(a >= b) return a;
	else return b;
}
typedef struct population
{
	int job[gongjian_num];	
	int time;			//�ӹ�ʱ�� 
	double fit;			//��Ӧ�� 
	double p;			//ѡ����� 
	double q;			//�ۼƸ��� 
}pop;
int search(int a, int b[])
{
	for(int k = 0; k < gongjian_num; k++)
	{
		if(a == b[k])
		return k;
	}
	return -1;
}

//��ʼ����Ⱥ 
void Init_population(pop p[])
{
	ifstream fin("data.txt");
	for(int i = 0; i < gongjian_num; i++)
		for(int j = 0; j < machine_num; j++)
			fin >> t[i][j];
	fin.close();
	for(int i = 0; i < maxn; i++)
	{
		vector <int> v;
	for(int j = 0; j < gongjian_num; j++)
		v.push_back(j);
	for(int j = 0; j < gongjian_num; j++)
	{
		int a = rand() % v.size();
		p[i].job[j] = v[a];
		v.erase(v.begin() + a);
	}	
	}
}

//��Ӧ�Ⱥ��� 
void fitness(pop p[])
{
	for(int i = 0; i < maxn; i++)
	{
		int c[gongjian_num][machine_num];
		fill(c[0], c[0] + gongjian_num * machine_num, 0);
		c[p[i].job[0]][0] = t[p[i].job[0]][0];
		for(int j = 1; j < machine_num; j++)
		{
			c[p[i].job[0]][j] = c[p[i].job[0]][j - 1] + t[p[i].job[0]][j];
		}
		for(int j = 1; j < gongjian_num; j++)
		{
			c[p[i].job[j]][0] = c[p[i].job[j - 1]][0] + t[p[i].job[j]][0];
		}
		for(int j = 1; j < gongjian_num; j++)
			for(int k = 1; k < machine_num; k++)
		{
			c[p[i].job[j]][k] = MAX(c[p[i].job[j - 1]][k], c[p[i].job[j]][k - 1]) + t[p[i].job[j]][k];
		}
		p[i].time = c[p[i].job[gongjian_num - 1]][machine_num - 1];
		p[i].fit = (1.0 / p[i].time) ;
	}
}

//ѡ�� 
void Select(pop p[])
{
	double sum_fit = 0;
	double sum_q = 0;
	for(int i = 0; i < maxn; i++)
	{
		sum_fit += p[i].fit;
	}
	for(int i = 0; i < maxn; i++)
	{
		p[i].p = p[i].fit / sum_fit;
		sum_q += p[i].p;
		p[i].q = sum_q;	
	}
	pop temp[maxn];
	for(int i = 0; i < maxn; i++)
	{
		int j = 0;
		double a = (rand() % 1000) * 0.001;
		if(a <= p[0].q)
		{
			temp[i] = p[0];
			continue;
		}
		else
		{
			for(j = 1; j < maxn; j++)
				if(a > p[j - 1].q && a <= p[j].q)
					break;
		}
		temp[i] = p[j];
	}
	for(int i = 0; i < maxn; i++)
		p[i] = temp[i];	
}

//���� 
void Crossover(pop p[])
{
	int n = (int)(pc * 10);
	while(n--)
	{
		int ia = 0, ib = 0;
		while(ia == ib)
		{
			ia = rand() % maxn;
			ib = rand() % maxn;
		}
		int a = 0, b = 0;
		while(a == b)
		{
			a = rand() % gongjian_num;
			b = rand() % gongjian_num; 
		}
		int arr[gongjian_num];
		int brr[gongjian_num];
		fill(arr, arr + gongjian_num, -1);
		fill(brr, brr + gongjian_num, -1);
		for(int j = a; j < b; j++)
		{
			arr[p[ia].job[j]] = p[ib].job[j];
			brr[p[ib].job[j]] = p[ia].job[j];
			int t = p[ia].job[j];
			p[ia].job[j] = p[ib].job[j];
			p[ib].job[j] = t;
		} 
		for(int j = 0; j < a; j++)
		{
			while(search(p[ia].job[j], arr) != -1)
			{
				p[ia].job[j] = search(p[ia].job[j],arr);
			}
		}
		for(int j = b; j < gongjian_num; j++)
		{
			while(search(p[ia].job[j], arr) != -1)
			{
				p[ia].job[j] = search(p[ia].job[j],arr);
			}
		}
		for(int j = 0; j < a; j++)
		{
			while(search(p[ib].job[j], brr) != -1)
			{
				p[ib].job[j] = search(p[ib].job[j],brr);
			}
		}
		for(int j = b; j < gongjian_num; j++)
		{
			while(search(p[ib].job[j], brr) != -1)
			{
				p[ib].job[j] = search(p[ib].job[j],brr);
			}
		}
	}
}



//����(�������) 
void Mutation(pop p[])
{
	int n = (int)(pm * 10);
	while(n--)
	{
		int index = rand() % maxn;
		int a = 0, b = 0;
		while(a == b)
		{
			a = rand() % gongjian_num;
			b = rand() % gongjian_num;
		}
		int temp = p[index].job[a];
		p[index].job[a] = p[index].job[b];
		p[index].job[b] = temp;	
	}
	
}


int main()
{
	pop p[maxn];
	srand((unsigned)(time(0)));
	Init_population(p);
	int cnt = 0;
	while(cnt < iteration)
	{
		fitness(p);
		int jubu_best = 9999;
		for(int i = 0; i < maxn; i++)
			if(p[i].time < jubu_best)
				jubu_best = p[i].time;
		cout<<"��"<<cnt+1<<"��,��̼ӹ�ʱ��Ϊ:"<<jubu_best<<endl;
		Select(p);
		Crossover(p);
		Mutation(p);
		for(int i = 0; i < maxn; i++)
			if(p[i].time <= best_time)
			{
				best_time = p[i].time;
				for(int j = 0; j < gongjian_num; j++)
				best[j] = p[i].job[j];
			}
		cnt++;
	}
	cout<<"��̼ӹ�ʱ��Ϊ�� "<<best_time<<endl;
	cout<<"���żӹ�����Ϊ�� ";
	for(int i = 0; i < gongjian_num; i++)
	{
		if(i == 0) cout<<best[i] + 1;
		else cout<<"-"<<best[i] + 1;
	}
	cout<<endl;	
