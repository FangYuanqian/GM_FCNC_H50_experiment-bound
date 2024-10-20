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
const double NPOT = 2e6;       // 总模拟次数
const double Nsim = 2e6;       // 总模拟次数
const double d = 128.0;        // 探测器距离目标的距离，单位为米
const double deltall = 12.0;   // 探测器的长度，单位为米
const double c = 299792458.0;

const double ξ_K = 1e-2;       // 几何因子
const double Pphi = 0.5;       // 假设的出现概率

// 计算标量粒子在 K 介子静止系中的动量
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
double Br_K_to_πφ(const std::string& K_type, double theta_squared, double p_phi_CM_value) {
    double m_K;
    if (K_type == "K+") {
        m_K = 0.493677;
    } else if (K_type == "K-") {
        m_K = 0.493677;
    } else if (K_type == "K_L") {
        m_K = 0.497611;
    } else {
        return 0;
    }

    if (K_type == "K_L") {
        return 5.7e-3 * (2 * p_phi_CM_value / m_K) * theta_squared;
    } else {
        return 1.6e-3 * (2 * p_phi_CM_value / m_K) * theta_squared;
    }
}

// 计算 N_phi 函数
double calculate_N_phi(double N_K_plus, double N_K_minus, double N_K_L, double theta_squared, double p_phi_CM_value) {
    double Br_K_plus = Br_K_to_πφ("K+", theta_squared, p_phi_CM_value);
    double Br_K_minus = Br_K_to_πφ("K-", theta_squared, p_phi_CM_value);
    double Br_K_L = Br_K_to_πφ("K_L", theta_squared, p_phi_CM_value);

    double N_phi = (NPOT / Nsim) * (N_K_plus * Br_K_plus + N_K_minus * Br_K_minus + N_K_L * Br_K_L) * ξ_K * Pphi;

    return N_phi;
}

// 从 CSV 文件中读取数据
std::vector<std::pair<double, double>> read_m_phi_theta_squared_from_csv(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line;
    std::vector<std::pair<double, double>> data;
    std::getline(file, line); // 跳过标题行

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

// 计算每个 CSV 文件的 N_phi，并输出结果（基于 theta^2 范围）
void calculate_N_phi_for_single_file(const std::vector<std::pair<double, double>>& file_data, 
                                     const std::string& file_name, 
                                     int num_points, double N_K_plus, double N_K_minus, double N_K_L) {
    // 获取 theta^2 的最小值和最大值
    double min_theta_squared = std::numeric_limits<double>::max();
    double max_theta_squared = std::numeric_limits<double>::lowest();
    for (const auto& [_, theta_squared] : file_data) {
        min_theta_squared = std::min(min_theta_squared, theta_squared);
        max_theta_squared = std::max(max_theta_squared, theta_squared);
    }

    double step = (max_theta_squared - min_theta_squared) / (num_points - 1);

    std::cout << "文件: " << file_name << " 的 N_phi 结果 (按 theta^2 计算):" << std::endl;

    // 计算每个区间的 N_phi
    for (int i = 0; i < num_points - 1; ++i) {
        double theta_squared_low = min_theta_squared + i * step;
        double theta_squared_high = min_theta_squared + (i + 1) * step;

        double min_N_phi = std::numeric_limits<double>::max();
        double max_N_phi = std::numeric_limits<double>::lowest();

        for (const auto& [m_phi, theta_squared] : file_data) {
            if (theta_squared >= theta_squared_low && theta_squared <= theta_squared_high) {
                double p_phi_CM_value = p_phi_CM(m_phi, 0.139, 0.493);  // 计算动量

                double N_phi = calculate_N_phi(N_K_plus, N_K_minus, N_K_L, theta_squared, p_phi_CM_value);

                // 检查 N_phi 的有效性
                if (!std::isnan(N_phi) && !std::isinf(N_phi)) {
                    min_N_phi = std::min(min_N_phi, N_phi);
                    max_N_phi = std::max(max_N_phi, N_phi);
                }
            }
        }

        // 输出到终端
        std::cout << "theta^2范围: [" << theta_squared_low << ", " << theta_squared_high 
                  << "] - N_phi最小值: " << min_N_phi << ", N_phi最大值: " << max_N_phi << std::endl;
    }
}

// 计算合并后的 N_phi，并输出到 CSV 文件
void calculate_combined_N_phi(const std::vector<std::pair<double, double>>& combined_data, 
                              const std::string& output_file, int num_points, double N_K_plus, 
                              double N_K_minus, double N_K_L) {
    // 获取 theta^2 的最小值和最大值
    double min_theta_squared = std::numeric_limits<double>::max();
    double max_theta_squared = std::numeric_limits<double>::lowest();
    for (const auto& [_, theta_squared] : combined_data) {
        min_theta_squared = std::min(min_theta_squared, theta_squared);
        max_theta_squared = std::max(max_theta_squared, theta_squared);
    }

    double step = (max_theta_squared - min_theta_squared) / (num_points - 1);

    // 打开输出文件
    std::ofstream out_file(output_file);
    out_file << std::fixed << std::setprecision(8);
    out_file << std::scientific;
    out_file << "theta_squared_low,theta_squared_high,min_N_phi,max_N_phi\n";

    // 计算合并数据的 N_phi
    for (int i = 0; i < num_points - 1; ++i) {
        double theta_squared_low = min_theta_squared + i * step;
        double theta_squared_high = min_theta_squared + (i + 1) * step;

        double min_N_phi = std::numeric_limits<double>::max();
        double max_N_phi = std::numeric_limits<double>::lowest();

        for (const auto& [m_phi, theta_squared] : combined_data) {
            if (theta_squared >= theta_squared_low && theta_squared <= theta_squared_high) {
                double p_phi_CM_value = p_phi_CM(m_phi, 0.139, 0.493);  // 计算动量
                
                double N_phi = calculate_N_phi(N_K_plus, N_K_minus, N_K_L, theta_squared, p_phi_CM_value);

                if (!std::isnan(N_phi) && !std::isinf(N_phi)) {
                    min_N_phi = std::min(min_N_phi, N_phi);
                    max_N_phi = std::max(max_N_phi, N_phi);
                }
            }
        }

        // 输出到文件
        out_file << theta_squared_low << "," << theta_squared_high << "," << min_N_phi << "," << max_N_phi << std::endl;
    }

    out_file.close();
}

int main() {
    // 定义 CSV 文件路径
    std::string csv_file1 = "/home/fyq/DIYbyfyq/Beam dump/CHARM/167/data/input/167PS191constraintsupperline.csv";
    std::string csv_file2 = "/home/fyq/DIYbyfyq/Beam dump/CHARM/167/data/input/167PS191constraintsdownline.csv";

    // 定义输出文件路径
    std::string output_file = "/home/fyq/DIYbyfyq/Beam dump/CHARM/167/data/167PS191output/167PS191NphiRangeByTheta_output.csv";

    // 读取 CSV 数据
    auto data_file1 = read_m_phi_theta_squared_from_csv(csv_file1);
    auto data_file2 = read_m_phi_theta_squared_from_csv(csv_file2);

    // 定义 K 介子的数量
    double N_K_plus = NPOT / 3;
    double N_K_minus = NPOT / 3;
    double N_K_L = NPOT / 3;

    // 计算每个文件的 N_phi
    calculate_N_phi_for_single_file(data_file1, "167PS191constraintsupperline.csv", 10, N_K_plus, N_K_minus, N_K_L);
    calculate_N_phi_for_single_file(data_file2, "167PS191constraintsdownline.csv", 10, N_K_plus, N_K_minus, N_K_L);

    // 合并数据并计算最终 N_phi
    std::vector<std::pair<double, double>> combined_data = data_file1;
    combined_data.insert(combined_data.end(), data_file2.begin(), data_file2.end());

    calculate_combined_N_phi(combined_data, output_file, 20, N_K_plus, N_K_minus, N_K_L);

    return 0;
}
