#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <sstream>

// 常量定义
const double NPOT = 2e6;       // 总模拟次数，表示模拟中每种类型的介子的总数
const double Nsim = 2e6;       // 总模拟次数，表示所有模拟的总次数
const double d = 128.0;        // 探测器距离目标的距离，单位为米，这里设置为 128 米
const double deltall = 12.0;        // 探测器的长度，单位为米，这里设置为 12 米
const double c = 299792458.0;

// 假设Nk&ξk比例系数（与实验无关），取为定值
// 假设P几乎与模型无关，是一个定量（0.5）
const double ξ_K = 1e-2;       // 几何因子，用于计算某些物理过程中的几何效应，这里设置为 0.01
const double P = 0.5;          // 假设的出现概率，表示介子出现的概率，这里设置为 0.5

 // 计算标量粒子在K介子静止系中的动量
double p_phi_CM(double m_phi, double m_pi, double m_K) {
    return (m_K / 2) * sqrt((1 - pow((m_phi + m_pi) / m_K, 2)) * (1 - pow((m_phi - m_pi) / m_K, 2)));
}

// 计算每种介子的衰变分支比 Br(K → πφ)
double Br_K_to_πφ(const std::string& K_type, double θ, double p_phi_CM_value) {
    // 定义不同类型 K 介子的质量
    double m_K;
    if (K_type == "K+") {
        m_K = 0.493677;  // K+ 的质量，单位：GeV/c^2
    } else if (K_type == "K-") {
        m_K = 0.493677;  // K- 的质量，单位：GeV/c^2
    } else if (K_type == "K_L") {
        m_K = 0.497614;  // K_L 的质量，单位：GeV/c^2
    } else {
        // 如果介子类型不匹配，返回 0
        return 0;
    }

    // 根据介子类型的质量计算分支比
    if (K_type == "K_L") {
        return 5.7e-3 * (2 * p_phi_CM_value / m_K) * pow(θ, 2);
    } else {
        return 1.6e-3 * (2 * p_phi_CM_value / m_K) * pow(θ, 2);
    }
}

// 从单个 CSV 文件中读取参数
std::vector<std::pair<double, double>> read_csv(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << file_path << std::endl;
        return {};
    }

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
            try {
                double m_phi = std::stod(tokens[0]);
                double θ_squared = std::stod(tokens[1]);
                data.emplace_back(m_phi, θ_squared);
            } catch (const std::exception& e) {
                std::cerr << "解析 CSV 文件中的数据时出错: " << e.what() << std::endl;
            }
        }
    }

    return data;
}

// 计算单个 CSV 文件中 N_phi 的范围
void calculate_N_phi_from_single_csv(const std::string& csv_file) {
    auto data = read_csv(csv_file);  // 从单个 CSV 文件中读取数据

    if (data.empty()) {
        std::cerr << "数据为空或无法读取数据。" << std::endl;
        return;
    }

    double min_N_phi = std::numeric_limits<double>::max();  // 初始化最小值为极大值
    double max_N_phi = std::numeric_limits<double>::lowest();  // 初始化最大值为极小值

    // 定义不同类型 K 介子的数量
    double N_K_plus = NPOT / 3;  // 假设 K+ 的数量为总模拟次数的三分之一
    double N_K_minus = NPOT / 3; // 假设 K- 的数量为总模拟次数的三分之一
    double N_K_L = NPOT / 3; // 剩余的数量给 K_L

    for (const auto& [m_phi, θ_squared] : data) {
        double θ = std::sqrt(θ_squared);  // 从 θ^2 计算 θ
        double p_phi_CM_value = p_phi_CM(m_phi, 0.139, 0.493);  // 计算动量

        // 计算每种介子的衰变分支比
        double Br_K_plus = Br_K_to_πφ("K+", θ, p_phi_CM_value);
        double Br_K_minus = Br_K_to_πφ("K-", θ, p_phi_CM_value);
        double Br_K_L = Br_K_to_πφ("K_L", θ, p_phi_CM_value);

        // 计算 N_phi
        double N_phi = (NPOT / Nsim) * (N_K_plus * Br_K_plus + N_K_minus * Br_K_minus + N_K_L * Br_K_L) * ξ_K * P;

        // 更新最小值和最大值
        if (N_phi < min_N_phi) {
            min_N_phi = N_phi;
        }
        if (N_phi > max_N_phi) {
            max_N_phi = N_phi;
        }
    }

    std::cout << "文件 " << csv_file << " 中 N_phi 的最小值: " << min_N_phi << std::endl;
    std::cout << "文件 " << csv_file << " 中 N_phi 的最大值: " << max_N_phi << std::endl;
}

// 程序入口点
int main() {
    std::vector<std::string> csv_files = {
        "/home/fyq/H5/data/167 REDLINE.csv",
        "/home/fyq/H5/data/167 DOWN LINE.csv",
        "/home/fyq/H5/data/167 UPPER LINE.csv"
    };

    for (const auto& csv_file : csv_files) {
        calculate_N_phi_from_single_csv(csv_file);
    }

    return 0;
}