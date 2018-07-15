#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <sstream>

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

bool GenOneTask(TaskDesc* task_desc, double edge_comp_frequency, double transmit_speed) {
  task_desc->len = std::rand() % 20 + 1;
  task_desc->load = std::rand() % 10 + 1;
  task_desc->fc = std::rand() % 30 + 1;
  // we set limit to inf to avoid reject situation
  task_desc->limit = std::numeric_limits<int>::max();
  // task_desc->limit = rand() % 10 + 1;

  task_desc->time_mobile = 1.0 * task_desc->len * task_desc->load / task_desc->fc;
  task_desc->time_edge = 1.0 * task_desc->len * task_desc->load / edge_comp_frequency;
  task_desc->time_transmit = 1.0 * task_desc->len / transmit_speed;

  if (task_desc->time_mobile > task_desc->limit &&
      (task_desc->time_edge + task_desc->time_transmit) > task_desc->limit) {
    // illegal task desc
    return false;
  } else {
    return true;
  }
}

void GenJonConf(int task_num, double edge_comp_frequency, double transmit_speed,
    std::vector<TaskDesc>& job_conf) {
  job_conf.resize(task_num);
  double cur_time = 0.0;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::exponential_distribution<double> distribution(exp_dist_para);
  std::srand(std::time(nullptr));
  for (int i = 0; i < task_num; ++i) {
    while (!GenOneTask(&job_conf[i], edge_comp_frequency, transmit_speed)) {}
    job_conf[i].gen_time = cur_time;
    cur_time += distribution(gen);
  }
}

void ReadJobConf4File(const std::string& file_name, std::vector<TaskDesc>& job_conf) {
  job_conf.clear();
  std::ifstream istrm(file_name);
  std::string line;
  while (std::getline(istrm, line)) {
    std::istringstream iss(line);
    TaskDesc task_desc;
    iss >> task_desc.len >> task_desc.load >> task_desc.fc >> task_desc.limit >>
      task_desc.time_mobile >> task_desc.time_edge >> task_desc.time_transmit >> task_desc.gen_time;
    job_conf.push_back(task_desc);
  }
  istrm.close();
}

void WriteJobConf2File(const std::string& file_name, const std::vector<TaskDesc>& job_conf) {
  std::fstream ostrm(file_name, std::ios::out);
  for (const auto& task_desc : job_conf) {
    ostrm << task_desc.len << " " << task_desc.load << " " << task_desc.fc << " " << task_desc.limit
      << " " << task_desc.time_mobile << " " << task_desc.time_edge << " "
      << task_desc.time_transmit << " " << task_desc.gen_time << std::endl;
  }
  ostrm.close();
}

void SolveByOneStepStrategy(const std::vector<TaskDesc>& job_conf) {
  double edge_time = 0;
  int reject_cnt = 0;
  double cost = 0;
  for (int i = 0; i < job_conf.size(); ++i) {
    const TaskDesc& task_desc = job_conf[i];
    if (task_desc.gen_time + task_desc.time_transmit >= edge_time) {
      double local_cost = task_desc.time_mobile;
      double edge_cost = task_desc.time_edge + task_desc.time_transmit;
      if (task_desc.limit < local_cost && task_desc.limit < edge_cost) {
        // in fact, this situation will not happen
        std::cout << "reject" << std::endl;
        ++reject_cnt;
      } else if (local_cost < edge_cost) {
        // std::cout << "mobile" << std::endl << "time is " << edge_time << std::endl;
        cost += local_cost;
      } else {
        edge_time = task_desc.time_transmit + task_desc.time_edge + task_desc.gen_time;
        // std::cout << "edge" << std::endl << "time is " << edge_time << std::endl;
        cost += edge_cost;
      }
    } else {
      double wait_time = edge_time - (task_desc.gen_time + task_desc.time_transmit);
      double local_cost = task_desc.time_mobile;
      double edge_cost = task_desc.time_edge  + task_desc.time_transmit + wait_time;
      if (task_desc.limit < local_cost && task_desc.limit < edge_cost) {
        // this situation may happen
        std::cout << "reject" << std::endl;
        ++reject_cnt;
      } else if (local_cost < edge_cost) {
        // std::cout << "mobile" << std::endl << "time is " << edge_time << std::endl;
        cost += local_cost;
      } else {
        edge_time = edge_time + task_desc.time_edge;
        // std::cout << "edge" << std::endl << "time is " << edge_time << std::endl;
        cost += edge_cost;
      }
    }
  }
  // std::cout << "reject: " << reject_cnt << std::endl;
  std::cout << "Cost calculated by One Step Strategy is " << cost << std::endl;
}

double CalcExpectedTimeOnEdge(double prev_gen_time, double time_transmit, double time_edge,
    double prev_end_time, double lambda) {
  if (prev_gen_time + time_transmit >= prev_end_time) {
    return time_transmit + time_edge;
  }
  double threshold = prev_end_time - prev_gen_time - time_transmit;
  // case1 : wait
  double expected_cost1 = (prev_end_time + time_edge - prev_gen_time) * (1 - std::exp(- lambda * threshold))
    + (std::exp(- lambda * threshold) * (lambda * threshold + 1) - 1) / lambda; 
  // case2 : no wait
  double expected_cost2 = (time_transmit + time_edge) * std::exp(- lambda * threshold);
  return expected_cost1 + expected_cost2;
}

void SolveByTwoStepStrategy(const std::vector<TaskDesc>& job_conf) {
  double lambda_1 = 1;
  double lambda_2 = 10;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(lambda_1, lambda_2);
  double lambda = dis(gen);
  double cur_time = 0.0;
  double sum_delta_time = 0.0;
  double prev_end_time = 0.0;
  double total_cost = 0.0;
  for (int i = 0; i < job_conf.size(); ++i) {
    const TaskDesc& task_desc = job_conf[i];
    sum_delta_time += task_desc.gen_time - cur_time;
    cur_time = task_desc.gen_time;
    // update lambda by MAP
    double opt_lambda = 1.0 * (i + 1) / sum_delta_time;
    if (opt_lambda < lambda_1) {
      lambda = lambda_1;
    } else if (opt_lambda > lambda_2) {
      lambda = lambda_2;
    } else {
      lambda = opt_lambda;
    }
    // assume next work share same properties
    // minimize the expected cost of current and next work
    // mobile + mobile
    double t0 = task_desc.time_mobile + task_desc.time_mobile;
    // mobile + edge
    double t1 = task_desc.time_mobile + CalcExpectedTimeOnEdge(task_desc.gen_time,
        task_desc.time_transmit, task_desc.time_edge, prev_end_time, lambda);
    double try_edge_start_time = std::max(prev_end_time, task_desc.gen_time +
        task_desc.time_transmit);
    double try_edge_cost = try_edge_start_time + task_desc.time_edge - task_desc.gen_time;
    // edge + mobile
    double t2 = try_edge_cost + task_desc.time_mobile;
    // edge + edge
    double t3 = try_edge_cost + CalcExpectedTimeOnEdge(task_desc.gen_time,
        task_desc.time_transmit, task_desc.time_edge, try_edge_start_time + task_desc.time_edge,
        lambda);
    if (std::min(t0, t1) < std::min(t2, t3)) {
      total_cost += task_desc.time_mobile;
    } else {
      total_cost += try_edge_cost;
      prev_end_time = try_edge_start_time + task_desc.time_edge;
    }
  }
  std::cout << "Estimated lambda : " << lambda << std::endl;
  std::cout << "Cost calculated by Two Step Strategy is " << total_cost << std::endl;
}

void SolveByBruteForce(const std::vector<TaskDesc>& job_conf) {
  int task_cnt = static_cast<int>(job_conf.size());
  if (task_cnt > 30) {
    std::cout << "Brute force only for problem size <= 30" << std::endl;
    return;
  }
  // reject situation is not considered in current brute force implementation
  double min_cost = std::numeric_limits<double>::max();
  for (int i = 0; i < (1 << task_cnt); ++i) {
    std::queue<TaskDesc> edge_queue;
    double cur_cost = 0;
    int cur_reject_cnt = 0;
    for (int j = 0; j < task_cnt; ++j) {
      if (i & (1 << j)) {
        edge_queue.push(job_conf[j]);
      } else {
        cur_cost += job_conf[j].time_mobile;
      }
    }
    double edge_time = 0;
    while (!edge_queue.empty()) {
      TaskDesc task = edge_queue.front();
      edge_queue.pop();
      double wait_time = std::max(edge_time - (task.gen_time + task.time_transmit), 0.0);
      cur_cost += wait_time + task.time_edge;
      edge_time = task.gen_time + task.time_transmit + wait_time + task.time_edge;
    }
    min_cost = std::min(cur_cost, min_cost);
  }
  std::cout << "Ground truth cost is " << min_cost << std::endl;
}

void SolveByRandom(const std::vector<TaskDesc>& job_conf) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  double total_cost = 0;
  std::queue<TaskDesc> edge_queue;
  for (int i = 0; i < job_conf.size(); ++i) {
    if (dis(gen) < 0.5) {
      edge_queue.push(job_conf[i]);
    } else {
      total_cost += job_conf[i].time_mobile;
    }
  }
  double edge_time = 0;
  while (!edge_queue.empty()) {
    TaskDesc task = edge_queue.front();
    edge_queue.pop();
    double wait_time = std::max(edge_time - (task.gen_time + task.time_transmit), 0.0);
    total_cost+= wait_time + task.time_edge;
    edge_time = task.gen_time + task.time_transmit + wait_time + task.time_edge;
  }
  std::cout << "Cost calculated by Random Strategy is " << total_cost << std::endl;
}

void PrintUsage() {
  std::cout << "Usage of this program: " << std::endl
    << "-g task_num transmit_speed edge_comp_frequency output_filename: "
    << "generate taskconf and store in file" << std::endl
    << "-o input_filename: "
    << "solve taskconf stored in input_filename by one step strategy" << std::endl
    << "-t input_filename: "
    << "solve taskconf stored in input_filename by two step strategy" << std::endl
    << "-b input_filename: "
    << "solve taskconf stored in input_filename by brute force" << std::endl
    << "-r input_filename: "
    << "solve taskconf stored in input_filename by random" << std::endl;
}

int main(int argc, char* argv[]) {
  if (argc == 6 && std::string(argv[1]) == "-g") {
    int task_num = std::stoi(std::string(argv[2]));
    double transmit_speed = std::stod(std::string(argv[3]));
    double edge_comp_frequency = std::stod(std::string(argv[4]));
    std::string output_filename = std::string(argv[5]);
    std::vector<TaskDesc> job_conf;
    GenJonConf(task_num, edge_comp_frequency, transmit_speed, job_conf);
    WriteJobConf2File(output_filename, job_conf);
  } else if (argc == 3) {
    std::string input_filename = std::string(argv[2]);
    std::vector<TaskDesc> job_conf;
    ReadJobConf4File(input_filename, job_conf);
    if (std::string(argv[1]) == "-o") {
      SolveByOneStepStrategy(job_conf);
    } else if (std::string(argv[1]) == "-t") {
      SolveByTwoStepStrategy(job_conf);
    } else if (std::string(argv[1]) == "-b") {
      SolveByBruteForce(job_conf);
    } else if (std::string(argv[1]) == "-r") {
      SolveByRandom(job_conf);
    } else {
      PrintUsage();
    }
  } else {
    PrintUsage();
  }
  return 0;
}
