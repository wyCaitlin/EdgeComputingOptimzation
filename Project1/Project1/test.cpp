#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <time.h>

using namespace std;

double Time_c[100] = { 0 };  //time task execute on mobile
double Energy_c[100] = { 0 };  //energy task execute on mobile
double Time_tran[100] = { 0 }; //task transmit time
double Time_wait[100] = { 0 };  //task wait time on edge
double Time_e[100] = { 0 };  //time task execute on edge

struct Task {
    int Len; //data length of task
    int Load;  //workload of task
    int fc;  //compute frequency of mobile
    int Limit;  //limit delay time of task
};

/*double gaussrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if (phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;
	X = X + 4;
	if (X < 0) {
		X = 0.0 - X;
	}
	return X;
}*/

void ProduceData(string& file_name, int produce_num) {
	ofstream ostrm(file_name, ios::out);
  for (int i = 0; i < produce_num; ++i) {
		int a = rand() % 20 + 1;
		int b = rand() % 10 + 1;
		int c = rand() % 30 + 1;
		int d = rand() % 10 + 1;
    ostrm << a << " " << b << " " << c << " " << d << endl << endl;
	}
	ostrm.close();
}

double GenExpDistribution(double lambda) {
  std::random_device rd;
  std::mt19937 gen(rd());
  exponential_distribution<double> distribution(lambda);
  return distribution(gen);
}

int main() {
	int m, tran, fe;  //task amount, transmit speed, compute frequency of edge
	cin >> m >> tran >> fe;
	double Time = 0;  //current time
	double proTime = 0;
  int reject_cnt = 0;
  double exp_distribution_lambda = 3;

  string file_name = "E:\\CppWork\\data.txt";
	ProduceData(file_name, m);
	ifstream istrm(file_name, std::ios::in);

	for (int i = 0; i < m; ++i) {
		cout << "task" << i + 1 << ":";
		double delTime = GenExpDistribution(exp_distribution_lambda);  //interval time between serial tasks
		//cout << delTime << endl;
		proTime = delTime + proTime;  //task produce time
		//cout << "task" << i << " proTime is " << proTime << endl;
		Task task_i;
		istrm >> task_i.Len >> task_i.Load >> task_i.fc >> task_i.Limit;

		Time_c[i] = 1.0 * task_i.Len * task_i.Load / task_i.fc;
		//cout << "compute on mobile comsume " << Time_c[i] << endl;
		Time_e[i] = 1.0 * task_i.Len * task_i.Load / fe;
		//cout << "compute on edge comsume " << Time_e[i] << endl;
		//Energy_c[i] = task_i.Len * task_i.Load * task_i.fc * task_i.fc;
		Time_tran[i] = 1.0 * task_i.Len / tran;
		//cout << "trans time is " << Time_tran[i] << endl;

		if (proTime + Time_tran[i] >= Time) {
      double local_cost = Time_c[i];
      double edge_cost = Time_e[i] + Time_tran[i];
      if (task_i.Limit < local_cost && task_i.Limit < edge_cost) {
        cout << "reject" << endl;
        ++reject_cnt;
      } else if (local_cost < edge_cost) {
        cout << "mobile" << endl << "time is " << Time << endl;
      } else {
        Time = Time_tran[i] + Time_e[i] + proTime;
        cout << "edge" << endl << "time is " << Time << endl;
      }
		} else {
      double wait_time = Time - (proTime + Time_tran[i]);
      double local_cost = Time_c[i];
      double edge_cost = Time_e[i] + Time_tran[i] + wait_time;

      if (task_i.Limit < local_cost && task_i.Limit < edge_cost) {
        cout << "reject" << endl;
        ++reject_cnt;
      } else if (local_cost < edge_cost) {
        cout << "mobile" << endl << "time is " << Time << endl;
      } else {
        Time = Time + Time_e[i];
        cout << "edge" << endl << "time is " << Time << endl;
      }
		}
		cout << endl;
	}
	istrm.close();
	cout << "reject: " << reject_cnt << endl;
	return 0;
}