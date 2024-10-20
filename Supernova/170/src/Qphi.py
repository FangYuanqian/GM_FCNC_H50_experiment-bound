import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.ticker as ticker

# 读取 CSV 文件
input_file_path = '/home/fyq/DIYbyfyq/Experimentbound/Supernova/170/data/data4.csv'
data = pd.read_csv(input_file_path)

# 提取数据
m_phi = data['m_phi(GeV)']
sin_theta_squared = data['sin_theta^2']
energy_loss_rate = data['energy loss rate']
 
# 创建绘图
plt.figure(figsize=(10, 6))

# 使用散点图表示 m_phi 和 sin_theta_squared，点的颜色代表 energy_loss_rate
scatter = plt.scatter(m_phi, sin_theta_squared, c=energy_loss_rate, cmap='viridis', s=50, alpha=0.8, edgecolors='w')

# 设置坐标轴的对数尺度
plt.xscale('log')
plt.yscale('log')

# 设置坐标轴范围
plt.xlim(1e-3, 10)
plt.ylim(1e-14, 1e-3)

# 设置纵坐标的刻度，显示每个对数级别
y_ticks = np.logspace(-14, -3, num=12)  # 生成从 1e-14 到 1e-3 的对数刻度
plt.yticks(y_ticks)  # 设置 y 轴的刻度

# 设置刻度格式，显示每个指数
plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))

# 设置刻度位置
plt.gca().xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
plt.gca().yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=12))

# 添加颜色条来显示 energy_loss_rate 的变化
cbar = plt.colorbar(scatter)
cbar.set_label('Energy Loss Rate (erg/s)', fontsize=14)  # 这里添加单位

# 设置颜色条的刻度为对数刻度
cbar.set_ticks(np.logspace(np.log10(energy_loss_rate.min()), np.log10(energy_loss_rate.max()), num=6))
cbar.ax.set_yticklabels([f'{tick:.1e}' for tick in cbar.get_ticks()])  # 使用科学计数法格式化刻度标签

# 设置颜色条刻度格式为对数刻度
cbar.ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
cbar.ax.yaxis.set_minor_formatter(ticker.NullFormatter())

# 添加标签和标题
plt.xlabel('m_phi (GeV)', fontsize=14)
plt.ylabel('sin_theta^2', fontsize=14)
plt.title('Scatter Plot of m_phi vs sin_theta^2 with Energy Loss Rate', fontsize=16)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# 获取文件名前缀
file_name = os.path.basename(input_file_path)  # 获取文件名
file_prefix = os.path.splitext(file_name)[0]  # 去掉扩展名

# 生成保存的文件名
output_file_name = f'{file_prefix}_scatter_plot_with_energy_loss_rate.png'

# 保存图像到输入文件所在文件夹
output_folder = os.path.dirname(input_file_path)  # 获取输入文件的文件夹路径
plt.savefig(os.path.join(output_folder, output_file_name))

# 显示图形
plt.show()
