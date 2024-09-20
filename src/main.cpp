#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <iomanip>
#include <cmath>
#include "/home/fyq/H5/include/NphiRangeByMassGM.h"
#include "/home/fyq/H5/include/NphiRangeByMass167.h"

const double NPOT = 2e6;

// 质量范围和实验 N_phi 区间结构体
struct MassRange {
    double m_phi_low;
    double m_phi_high;
    double N_phi_min;
    double N_phi_max;
};

// 从 CSV 文件中读取质量范围和实验 N_phi 区间数据
std::vector<MassRange> read_mass_ranges_from_csv(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line;
    std::vector<MassRange> ranges;

    // 跳过标题行
    std::getline(file, line);

    // 读取数据
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> tokens;

        while (std::getline(ss, item, ',')) {
            tokens.push_back(item);
        }

        if (tokens.size() == 4) {
            double m_phi_low = std::stod(tokens[0]);
            double m_phi_high = std::stod(tokens[1]);
            double N_phi_min = std::stod(tokens[2]);
            double N_phi_max = std::stod(tokens[3]);
            ranges.push_back({m_phi_low, m_phi_high, N_phi_min, N_phi_max});
        }
    }
    return ranges;
}

// 扫描质量和角度范围，并计算 N_phi
void scan_mass_and_theta(const std::vector<MassRange>& ranges, const std::string& output_file) {
    std::ofstream out_file(output_file);
    if (!out_file) {
        std::cerr << "无法打开文件 " << output_file << " 进行写入。" << std::endl;
        return;
    }

    out_file << "质量范围下限,质量范围上限,m_H5,theta_H,N_phi" << std::endl;
    out_file << std::fixed << std::setprecision(10);

     // 定义不同类型 K 介子的数量
        double N_K_plus = NPOT / 3;  // 假设 K+ 的数量为总模拟次数的三分之一
        double N_K_minus = NPOT / 3; // 假设 K- 的数量为总模拟次数的三分之一
        double N_K_L = NPOT / 3;     // 剩余的数量给 K_L

    // 在每个质量范围内扫描
    for (const auto& range : ranges) {
        double m_phi_low = range.m_phi_low;
        double m_phi_high = range.m_phi_high;
        double N_phi_min = range.N_phi_min;
        double N_phi_max = range.N_phi_max;

        // 扫描质量范围内的每个点
        for (double m_H5 = m_phi_low; m_H5 <= m_phi_high; m_H5 += (m_phi_high - m_phi_low) / 10) { // 这里将范围分成了10个点
            // 扫描 θH 的范围
            for (double theta_H = -M_PI; theta_H <= M_PI; theta_H += 0.1) { // 这里将 θH 的范围从 -π 到 π 划分成了0.1的步长
                // 使用 NphiRangeByMass_GM.cpp 中的函数计算 N_phi
                double N_phi = calculate_N_phi(N_K_plus, N_K_minus, N_K_L, theta_H);

                // 检查计算的 N_phi 是否在实验 N_phi 区间内
                if (N_phi >= N_phi_min && N_phi <= N_phi_max) {
                    // 写入到输出文件
                    out_file << m_phi_low << "," << m_phi_high << "," << m_H5 << "," << theta_H << "," << N_phi << std::endl;
                }
            }
        }
    }
}

int main() {
    // 输入文件和输出文件路径
    std::string input_file = "/home/fyq/H5/data/N_phi_output_167.csv";
    std::string output_file = "/home/fyq/H5/data/N_phi_scanned_results.csv";

    // 读取质量范围数据
    auto ranges = read_mass_ranges_from_csv(input_file);

    // 扫描质量和角度范围，并计算 N_phi
    scan_mass_and_theta(ranges, output_file);

    return 0;
}
