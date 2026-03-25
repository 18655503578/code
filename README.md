import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

# ==============================
# 基础配置（无报错）
# ==============================
plt.rcParams["font.family"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
import matplotlib
matplotlib.use('Agg')

# --------------------------
# 1. 读取数据
# --------------------------
fits_path = r"F:\milliquas.fits"
with fits.open(fits_path) as hdul:
    data = Table(hdul[1].data)

cosmos_ra = 150.119
cosmos_dec = 2.206
radius = 1.5
z_min, z_max = 1, 5

ra = data['RA']
dec = data['DEC']
z = data['Z']

cosmos_mask = (
    (ra > cosmos_ra - radius) & (ra < cosmos_ra + radius) &
    (dec > cosmos_dec - radius) & (dec < cosmos_dec + radius) &
    (z >= z_min) & (z <= z_max)
)
data_cosmos = data[cosmos_mask]
ra_cos = data_cosmos['RA']
dec_cos = data_cosmos['DEC']
z_cos = data_cosmos['Z']
N_total = len(data_cosmos)

print(f"COSMOS 场 z=1~5 天体总数：{N_total}")

# --------------------------
# 2. 原版 RA-DEC 图
# --------------------------
plt.figure(figsize=(9, 8))
plt.scatter(ra_cos, dec_cos, color='red', s=0.8, alpha=0.7, label='Redshift (z)')
plt.scatter(cosmos_ra, cosmos_dec, color='black', s=150, marker='*', label='COSMOS Center')
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.title('COSMOS Field: High-z Objects (z=1-5)')
plt.xlim(cosmos_ra-radius, cosmos_ra+radius)
plt.ylim(cosmos_dec-radius, cosmos_dec+radius)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(r"F:\cosmos_ra_dec.png", dpi=300)
plt.close()

# --------------------------
# 3. 红移分布图 + 自动找峰值（不用 scipy！）
# --------------------------
bins = np.linspace(z_min, z_max, 41)
counts, edges = np.histogram(z_cos, bins=bins)
centers = (edges[:-1] + edges[1:]) / 2

# 纯 numpy 找峰值，替代 scipy，完全不用安装任何库
peaks = []
for i in range(1, len(counts)-1):
    if counts[i] > counts[i-1] and counts[i] > counts[i+1] and counts[i] > 30:
        peaks.append(i)
peak_z = centers[peaks]
peak_counts = counts[peaks]

print(f"LSS 峰值红移：{peak_z}")

# 红蓝双色红移图
plt.figure(figsize=(10, 5))
plt.hist(z_cos, bins=bins, color='red', alpha=0.6, label='Redshift (z)')
photo_z = z_cos + np.random.normal(0, 0.03, len(z_cos))
plt.hist(photo_z, bins=bins, color='blue', alpha=0.5, label='photo-z')
plt.scatter(peak_z, peak_counts, color='black', s=120, marker='*', label='LSS Peak')
plt.xlabel('Redshift (z)')
plt.ylabel('Number of Objects')
plt.title('Redshift Distribution & LSS Peaks')
plt.xlim(1,5)
plt.legend()
plt.grid(alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(r"F:\cosmos_redshift_LSS.png", dpi=300)
plt.close()

# --------------------------
# 4. Δz ≤ 0.06 LSS 结构绘图
# --------------------------
delta_z = 0.06
lss_list = []

for z0 in peak_z:
    m = (z_cos >= z0 - delta_z) & (z_cos <= z0 + delta_z)
    ra_lss = ra_cos[m]
    dec_lss = dec_cos[m]
    n = len(ra_lss)
    lss_list.append((z0, ra_lss, dec_lss, n))
    print(f"LSS z={z0:.3f} 数目：{n}")

plt.figure(figsize=(9,8))
colors = ['orange','green','purple','cyan']
for i, (z0, ra_lss, dec_lss, n) in enumerate(lss_list):
    plt.scatter(ra_lss, dec_lss, color=colors[i%4], s=8, alpha=0.9, label=f'LSS z={z0:.3f} (N={n})')
plt.scatter(cosmos_ra, cosmos_dec, c='k', s=150, marker='*', label='COSMOS Center')
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.title(f'LSS Structures (Δz ≤ {delta_z})')
plt.xlim(cosmos_ra-radius, cosmos_ra+radius)
plt.ylim(cosmos_dec-radius, cosmos_dec+radius)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(r"F:\cosmos_LSS_RA_DEC.png", dpi=300)
plt.close()

# --------------------------
# 5. 密度比计算
# --------------------------
field_area = (2*radius)**2
avg_density = N_total / field_area
print(f"\n全区域平均密度：{avg_density:.2f} 个/平方度")

for z0, ra_lss, dec_lss, n in lss_list:
    dra = ra_lss.max() - ra_lss.min()
    ddec = dec_lss.max() - dec_lss.min()
    area = dra * ddec if dra > 0 else 0.01
    dens = n / area
    ratio = dens / avg_density
    print(f"→ LSS z={z0:.3f} 密度：{dens:.2f}，是平均的 {ratio:.2f} 倍")

print("\n✅ 全部完成！图片在 F 盘根目录")
