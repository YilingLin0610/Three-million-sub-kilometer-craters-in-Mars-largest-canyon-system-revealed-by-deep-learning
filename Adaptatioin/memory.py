# -*- coding: utf-8 -*-
import re
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

# ====== 修改这里，填你的日志文件路径与起始日期 ======
log_path = r"G:\Mars\Mars_SFS_batch\mem_usage.log"   # 用原始字符串 r"..."
out_png = "memory_usage.png"
START_FROM = datetime(2025, 9, 12, 0, 0, 0)          # 从 2025-09-11 起
# ===================================================

def parse_mem_from_log(path, start_from=None):
    """
    读取日志，解析 [YYYY-MM-DD HH:MM:SS] 当前内存占用：XX(或XX.XX)%
    如提供 start_from（datetime），则仅返回 >= start_from 的记录。
    """
    pat = re.compile(
        r'^\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]\s*当前内存占用[:：]\s*([0-9]+(?:\.[0-9]+)?)%',
        re.UNICODE
    )
    times, vals = [], []
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            m = pat.search(line)
            if not m:
                continue
            t = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S")
            v = float(m.group(2))
            if (start_from is None) or (t >= start_from):
                times.append(t)
                vals.append(v)
    # 按时间排序
    if times:
        times, vals = zip(*sorted(zip(times, vals)))
    return list(times), list(vals)

def plot_mem(times, vals, out_png):
    if not times:
        print("日志中（过滤后）没有找到任何“当前内存占用”记录")
        return
    plt.figure(figsize=(9, 4.5))
    plt.plot(times, vals, marker='o', linewidth=1.5)  # 不指定颜色
    plt.title("内存占比随时间变化（起始：%s）" % times[0].strftime("%Y-%m-%d %H:%M:%S"))
    plt.xlabel("时间")
    plt.ylabel("内存占用 (%)")
    plt.grid(True, alpha=0.3)
    ax = plt.gca()
    ax.xaxis.set_major_formatter(DateFormatter("%H:%M:%S\n%Y-%m-%d"))
    plt.tight_layout()
    #plt.savefig(out_png, dpi=150)
    # 如需弹窗预览再取消注释：
    plt.show()
    print(f"图像已保存到: {out_png}")

if __name__ == "__main__":
    times, vals = parse_mem_from_log(log_path, start_from=START_FROM)
    plot_mem(times, vals, out_png)

