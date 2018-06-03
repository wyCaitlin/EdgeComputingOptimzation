#include<cstdio>
#include<cmath>
#include<iostream>
#include<fstream>
#include<queue>
#include<deque>

using namespace std;

struct Worker {
	int tp;
	int te;
	int tc;
	int tw;
};

int T[1000] = { 0 };

int proData(int N) {
	int s, a = 0, b, c;
	ofstream ostrm("E:\\CppWork\\test.txt", ios::out);
	for (int i = 0; i < N; i++) {
		s = rand() % 5 + 1;
		a = a + s;
		b = rand() % 10 + 1;
		c = rand() % 100 +1;
		while (c <= 5 * b) {
			c = rand() % 100 + 1;
		}
		ostrm << a << " " << b << " " << c << '\n' << endl;
	}
	ostrm.close();
	return 0;
}
  
int proTw(queue<Worker>& Q) {
	int tw = 0;
	int curT = 0;
	while (!Q.empty()) {
		Worker u = Q.front();
		Q.pop();
		if (curT < u.tp) {
			curT = u.tp;
		}
		else {
			u.tw = curT - u.tp;
		}
		curT += u.te;
	}
	return tw;
}

int f(int n, Worker w) {
	int ret = numeric_limits<int>::max();
	for (int i = 0; i < (1 << n); ++i) {
		int cost = 0;
		queue<Worker> Q;
		for (int j = 0; j < n; ++j) {
			if (i & (1 << j)) {
				cost += w.tc;
			}
			else {
				Q.push(w);
			}
		}
		cost += proTw(Q);
	}
	return ret;
}

int main() {
	int m;
	int x[1000] = { 0 };
	cin >> m;
	proData(m);

	ifstream istrm("E:\\CppWork\\test.txt", std::ios::in);
	for (int i = 0; i < m; ++i) {
		Worker w_i;
		istrm >> w_i.tp >> w_i.te >> w_i.tc;
		T[i] = x[i] * (w_i.te + w_i.tw) + (1 - x[i]) * w_i.tc;
	}

	for (int j = 0; j < m; ++j) {
		cout << "task" << j << " " << T[j] << endl;
	}
	
	istrm.close();
	return 0;
}