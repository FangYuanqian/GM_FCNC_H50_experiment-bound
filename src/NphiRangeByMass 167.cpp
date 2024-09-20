#include "/home/fyq/H5/include/NphiRangeByMass167.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <sstream>
#include <iomanip>
#include <utility>

// 常量定义
const double NPOT = 2e6;       // 总模拟次数，表示模拟中每种类型的介子的总数
const double Nsim = 2e6;       // 总模拟次数，表示所有模拟的总次数
const double d = 128.0;        // 探测器距离目标的距离，单位为米，这里设置为 128 米
const double deltall = 12.0;        // 探测器的长度，单位为米，这里设置为 12 米
const double c = 299792458.0;

// 假设P几乎与模型无关，是一个定量（0.5）
const double ξ_K = 1e-2;       // 几何因子
const double P = 0.5;          // 假设的出现概率，表示介子出现的概率，这里设置为 0.5

 // 计算标量粒子在K介子静止系中的动量
double p_phi_CM(double m_phi, double m_pi, double m_K) {
    double term1 = (m_phi + m_pi) / m_K;
    double term2 = (m_phi - m_pi) / m_K;
    double arg1 = 1 - pow(term1, 2);
    double arg2 = 1 - pow(term2, 2);

    if (arg1 < 0 || arg2 < 0) {
        return 0;  
    }

    return (m_K / 2) * sqrt(arg1 * arg2);
}


// 计算每种介子的衰变分支比 Br(K → πφ)
double Br_K_to_πφ(const std::string& K_type, double θ, double p_phi_CM_value) {
    double m_K;
    if (K_type == "K+") {
        m_K = 0.493677;  // K+ 的质量，单位：GeV/c^2
    } else if (K_type == "K-") {
        m_K = 0.493677;  // K- 的质量，单位：GeV/c^2
    } else if (K_type == "K_L") {
        m_K = 0.497611;  // K_L 的质量，单位：GeV/c^2
    } else {
        return 0;
    }

    // 根据介子类型的质量计算分支比
    if (K_type == "K_L") {
        return 5.7e-3 * (2 * p_phi_CM_value / m_K) * pow(θ, 2);
    } else {
        return 1.6e-3 * (2 * p_phi_CM_value / m_K) * pow(θ, 2);
    }
}

// 计算 N_phi 函数
double calculate_N_phi(double N_K_plus, double N_K_minus, double N_K_L, double θ, double p_phi_CM_value) {
    // 计算每种介子的衰变分支比
    double Br_K_plus = Br_K_to_πφ("K+", θ, p_phi_CM_value);
    double Br_K_minus = Br_K_to_πφ("K-", θ, p_phi_CM_value);
    double Br_K_L = Br_K_to_πφ("K_L", θ, p_phi_CM_value);

    // 计算 N_phi
    double N_phi = (NPOT / Nsim) * (N_K_plus * Br_K_plus + N_K_minus * Br_K_minus + N_K_L * Br_K_L) * ξ_K * P;

    // 返回 N_phi
    return N_phi;
}

// 从 CSV 文件中读取数据
std::vector<std::pair<double, double>> read_m_phi_theta_squared_from_csv(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line;
    std::vector<std::pair<double, double>> data;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> tokens;

        while (std::getline(ss, item, ',')) {
            tokens.push_back(item);
        }

        if (tokens.size() == 2) {
            double m_phi = std::stod(tokens[0]);
            double theta_squared = std::stod(tokens[1]);
            data.emplace_back(m_phi, theta_squared);
        }
    }
    return data;
}

// 计算每个区间内的 N_phi 的上下边缘
void calculate_N_phi_in_ranges(const std::string& csv_file, const std::string& output_file) {
    auto data = read_m_phi_theta_squared_from_csv(csv_file);

    // 获取 m_phi 的最小值和最大值
    double min_mass = std::numeric_limits<double>::max();
    double max_mass = std::numeric_limits<double>::lowest();
    for (const auto& [m_phi, _] : data) {
        min_mass = std::min(min_mass, m_phi);
        max_mass = std::max(max_mass, m_phi);
    }

    int num_points = 20;
    double step = (max_mass - min_mass) / (num_points - 1);

    // 打开输出文件
    std::ofstream out_file(output_file);
    if (!out_file) {
        std::cerr << "无法打开文件 " << output_file << " 进行写入。" << std::endl;
        return;
    }

    out_file << "质量范围下限,质量范围上限,N_phi最小值,N_phi最大值" << std::endl;
    out_file << std::fixed << std::setprecision(10);

    for (int i = 0; i < num_points - 1; ++i) {
        double m_phi_low = min_mass + i * step;
        double m_phi_high = min_mass + (i + 1) * step;

        double min_N_phi = std::numeric_limits<double>::max();
        double max_N_phi = std::numeric_limits<double>::lowest();

        // 定义不同类型 K 介子的数量
        double N_K_plus = NPOT / 3;  // 假设 K+ 的数量为总模拟次数的三分之一
        double N_K_minus = NPOT / 3; // 假设 K- 的数量为总模拟次数的三分之一
        double N_K_L = NPOT / 3;     // 剩余的数量给 K_L

        // 计算每个区间的 N_phi
        for (const auto& [m_phi, θ_squared] : data) {
            if (m_phi >= m_phi_low && m_phi <= m_phi_high) {
                double θ = std::sqrt(θ_squared);  // 从 θ^2 计算 θ
                double p_phi_CM_value = p_phi_CM(m_phi, 0.139, 0.493);  // 计算动量

                // 调用 calculate_N_phi 函数来计算 N_phi
                double N_phi = calculate_N_phi(N_K_plus, N_K_minus, N_K_L, θ, p_phi_CM_value);

                // 检查 N_phi 的有效性
                if (std::isnan(N_phi) || std::isinf(N_phi)) {
                    std::cerr << "检测到无效的 N_phi 值: " << N_phi << std::endl;
                }

                // 更新最小值和最大值
                min_N_phi = std::min(min_N_phi, N_phi);
                max_N_phi = std::max(max_N_phi, N_phi);
            }
        }

        // 写入 CSV 文件
        out_file << m_phi_low << "," << m_phi_high << "," 
                 << min_N_phi << "," << max_N_phi << std::endl;
    }

    out_file.close();  // 关闭文件
}

// 程序入口点
int main() {
    std::string csv_file = "/home/fyq/H5/data/167 REDLINE.csv";  
    std::string output_file = "/home/fyq/H5/data/N_phi_output 167.csv";  // 输出 CSV 文件路径
    calculate_N_phi_in_ranges(csv_file, output_file);
    return 0;
}


