import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 读取 C++ 生成的 CSV 文件
data_file = "/home/fyq/DIYbyfyq/H5/data/GMoutput/GM_N_phi.csv"
df = pd.read_csv(data_file)

# 创建一个散点图
plt.figure(figsize=(10, 6))

# 绘制散点图，确保 N_phi 用于颜色映射
scatter = plt.scatter(df['m_H5'], df['sin(theta_H)'], c=df['N_phi'], cmap='plasma', s=1, alpha=0.6)

# 添加颜色条
cbar = plt.colorbar(scatter, label="N_phi")
cbar.ax.tick_params(labelsize=8)  # 调整颜色条字体大小

# 设置坐标轴为对数坐标
plt.xscale('log')
plt.yscale('log')

# 设置坐标轴范围
plt.xlim(left=1e-2, right=1)  # 设置 x 轴范围
plt.ylim(bottom=1e-8, top=1)   # 设置 y 轴范围

# 设置刻度
plt.xticks([1e-2, 0.1, 1], ['0.01', '0.1', '1']) 
plt.yticks([1e-8, 1e-6, 1e-4, 1e-2, 1], 
           ['1e-8', '1e-6', '1e-4', '1e-2', '1'])

# 设置图表标题和标签
plt.title('N_phi Distribution Plot')
plt.xlabel('m_H5 (GeV)')
plt.ylabel('sin(theta_H)')

# 保存图像或显示图表
plt.savefig("/home/fyq/DIYbyfyq/H5/data/GMoutput/GM_Nphi_scatter.png")  # 保存为图片
plt.show()  # 显示图表
