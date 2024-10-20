import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 读取 C++ 生成的 CSV 文件
data_file = "/home/fyq/DIYbyfyq/Beam dump/CHARM/167/data/GMoutput/GM_N_phi.csv"
df = pd.read_csv(data_file)

# 打印 DataFrame 的列名
print("DataFrame 列名:", df.columns.tolist())

# 计算 N_phi 的最小值和最大值
min_N_phi_GM = df['N_phi_GM'].min()
max_N_phi_GM = df['N_phi_GM'].max()

# 创建一个散点图
plt.figure(figsize=(10, 6))

# 绘制散点图，确保 N_phi 用于颜色映射
scatter = plt.scatter(df['m_H5'], df['sin_theta_H'], c=df['N_phi_GM'], cmap='plasma', s=0.1, alpha=0.6, vmin=min_N_phi_GM, vmax=max_N_phi_GM)

# 添加颜色条
cbar = plt.colorbar(scatter, label="N_phi_GM")
cbar.ax.tick_params(labelsize=8)  # 调整颜色条字体大小

# 设置坐标轴为对数坐标
plt.xscale('log')
plt.yscale('log')

# 设置坐标轴范围
plt.xlim(left=1e-2, right=1)  # 设置 x 轴范围
plt.ylim(bottom=1e-4, top=1)   # 设置 y 轴范围

# 设置刻度
plt.xticks([1e-2, 0.1, 1], ['0.01', '0.1', '1']) 
plt.yticks([1e-4, 1e-2, 1], ['1e-4', '1e-2', '1'])

# 设置图表标题和标签
plt.title('N_phi Distribution Plot')
plt.xlabel('m_H5 (GeV)')
plt.ylabel('sin(theta_H)')

# 为每个数量级选择一个点标注 N_phi 数值
log_N_phi_levels = np.arange(np.floor(np.log10(min_N_phi_GM)), np.ceil(np.log10(max_N_phi_GM)) + 1)

for log_level in log_N_phi_levels:
    target_value = 10 ** log_level
    closest_point = df.iloc[(df['N_phi_GM'] - target_value).abs().argmin()]
    
    m_phi = closest_point['m_H5']
    theta_sq = closest_point['sin_theta_H']
    N_phi = closest_point['N_phi_GM']  # 修正为 N_phi_GM
    
    plt.text(m_phi, theta_sq, f'{N_phi:.1e}', fontsize=8, color='red', ha='center', va='bottom')

# 保存图像或显示图表
plt.savefig("/home/fyq/DIYbyfyq/Beamdump/CHARM/167/data/GMoutput/GM_Nphi_scatter.png")  # 保存为图片
plt.show()  # 显示图表
