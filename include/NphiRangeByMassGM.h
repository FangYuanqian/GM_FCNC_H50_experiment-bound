#ifndef NPHIRANGEBYMASSGM_H
#define NPHIRANGEBYMASSGM_H

// 常量定义
extern const double GF; 
extern const double lambda;
extern const double A_K_pm; 
extern const double f0_K_pm_pi; 
extern const double m_s;  // 单位为 GeV/c^2
extern const double m_d;  // 单位为 GeV/c^2
extern const double v;    // v 的值, 单位：GeV
extern const double Gamma_K_pm; // K^\pm 的衰变宽度，单位为 GeV
extern const double PI;

// 函数声明
// 计算 \xi_{\phi}^{ds}
double calculate_xi_phi_ds(double m_H5, double theta_H);

// 计算 \xi_W^5
double calculate_xi_W(double m_H5, double theta_H);

// 计算分支比 Br(K^\pm \rightarrow \pi^\pm H_5^0)
double calculate_Br_K_pm(double m_K_pm, double m_pi_pm, double m_H5, double theta_H);

// 计算分支比 Br(K_L \to \pi^0 H_5^0)
double calculate_Br_K_L(double m_K_L, double m_pi_0, double m_H5, double theta_H);

// 计算 N_phi 函数
double calculate_N_phi(double N_K_plus, double N_K_minus, double N_K_L, double theta_H);

#endif // NPHIRANGEBYMASSGM_H
