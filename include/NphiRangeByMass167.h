#ifndef NphiRangeByMass167H  // 防止重复包含头文件的保护符
#define NphiRangeByMass167H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <sstream>
#include <iomanip>

// 常量定义
extern const double NPOT;       // 总模拟次数，表示模拟中每种类型的介子的总数
extern const double Nsim;       // 总模拟次数，表示所有模拟的总次数
extern const double d;          // 探测器距离目标的距离，单位为米
extern const double deltall;    // 探测器的长度，单位为米
extern const double c;          // 光速，单位为米每秒
extern const double ξ_K;        // 几何因子
extern const double P;          // 出现概率

// 函数声明
double p_phi_CM(double m_phi, double m_pi, double m_K);  // 计算动量
double Br_K_to_πφ(const std::string& K_type, double θ, double p_phi_CM_value);  // 计算衰变分支比
std::vector<std::pair<double, double>> read_m_phi_theta_squared_from_csv(const std::string& file_path);  // 读取CSV文件
void calculate_N_phi_in_ranges(const std::string& csv_file);  // 计算每个区间的 N_phi 上下边缘
double calculate_N_phi(double N_K_plus, double N_K_minus, double N_K_L, double θ, double p_phi_CM_value);

#endif  // NphiRangeByMass167H
