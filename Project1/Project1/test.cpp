#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <time.h>

const double exp_dist_para = 3.0;

struct TaskDesc {
  int len; // data length of task
  int load;  // workload of task
  int fc;  // compute frequency of mobile
  int limit;  // limit delay time of task

  double time_mobile;
  double time_transmit;
  double time_wait;
  double time_edge;

  double gen_time;
};

void GenJonConf(int task_num, double edge_comp_frequency, double transmit_speed,
    std::vector<TaskDesc>& job_conf) {
  job_conf.resize(task_num);
  double cur_time = 0.0;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::exponential_distribution<double> distribution(exp_dist_para);
  for (int i = 0; i < task_num; ++i) {
    TaskDesc& task_desc = job_conf[i];
    task_desc.len = rand() % 20 + 1;
    task_desc.load = rand() % 10 + 1;
    task_desc.fc = rand() % 30 + 1;
    task_desc.limit = rand() % 10 + 1;

    task_desc.time_mobile = 1.0 * task_desc.len * task_desc.load / task_desc.fc;
    task_desc.time_edge = 1.0 * task_desc.len * task_desc.load / edge_comp_frequency;
    task_desc.time_transmit = 1.0 * task_desc.len / transmit_speed;

    task_desc.gen_time = cur_time;
    cur_time += distribution(gen);
  }
}

void WriteJobConf2File(const std::string& file_name, const std::vector<TaskDesc>& job_conf) {
  std::fstream ostrm(file_name, std::ios::out);
  for (const auto& task_desc : job_conf) {
    ostrm << task_desc.len << " " << task_desc.load << " " << task_desc.fc << " " << task_desc.limit
      << " " << task_desc.time_mobile << " " << task_desc.time_edge << " "
      << task_desc.time_transmit << " " << task_desc.gen_time << std::endl << std::endl;
  }
	ostrm.close();
}

void SolveByOneStepStrategy(const std::vector<TaskDesc>& job_conf) {
  double edge_time = 0;
  int reject_cnt = 0;
  for (int i = 0; i < job_conf.size(); ++i) {
    const TaskDesc& task_desc = job_conf[i];
    if (task_desc.gen_time + task_desc.time_transmit >= edge_time) {
      double local_cost = task_desc.time_mobile;
      double edge_cost = task_desc.time_edge + task_desc.time_transmit;
      if (task_desc.limit < local_cost && task_desc.limit < edge_cost) {
        std::cout << "reject" << std::endl;
        ++reject_cnt;
      } else if (local_cost < edge_cost) {
        std::cout << "mobile" << std::endl << "time is " << edge_time << std::endl;
      } else {
        edge_time = task_desc.time_transmit + task_desc.time_edge + task_desc.gen_time;
        std::cout << "edge" << std::endl << "time is " << edge_time << std::endl;
      }
    } else {
      double wait_time = edge_time - (task_desc.gen_time + task_desc.time_transmit);
      double local_cost = task_desc.time_mobile;
      double edge_cost = task_desc.time_edge  + task_desc.time_transmit + wait_time;
      if (task_desc.limit < local_cost && task_desc.limit < edge_cost) {
       std::cout << "reject" << std::endl;
       ++reject_cnt;
      } else if (local_cost < edge_cost) {
        std::cout << "mobile" << std::endl << "time is " << edge_time << std::endl;
      } else {
        edge_time = edge_time + task_desc.time_edge;
        std::cout << "edge" << std::endl << "time is " << edge_time << std::endl;
      }
    }
    std::cout << std::endl;
  }
  std::cout << "reject: " << reject_cnt << std::endl;
}

int main() {
  int task_num;
  double transmit_speed;
  double edge_comp_frequency;
  std::cin >> task_num >> transmit_speed >> edge_comp_frequency;

  std::string file_name = "/home/zhuyi/job_conf";
  std::vector<TaskDesc> job_conf;
  GenJonConf(task_num, edge_comp_frequency, transmit_speed, job_conf);
  WriteJobConf2File(file_name, job_conf);
  SolveByOneStepStrategy(job_conf);
	return 0;
}
