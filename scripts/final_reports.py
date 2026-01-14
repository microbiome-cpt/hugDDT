import pandas as pd
import re
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from adjustText import adjust_text
import numpy as np
import matplotlib.ticker as mticker

# --- 1. CheckM2 ---
q = pd.read_csv(sys.argv[1], sep="\t") #quality_report.tsv

# порядок бинов
order = list(q['Name'])

# нужные колонки
q_sub = q[['Name', 'Completeness_Specific', 'Contamination',
           q.columns[5],  # Contig_N50
           q.columns[7]]  # Genome_Size
          ].copy()
q_sub.columns = ['Name', 'Completeness', 'Contamination', 'Contig_N50', 'Genome_Size']

# --- 2. CoverM ---
c = pd.read_csv(sys.argv[2], sep="\t") #genome_cover.tsv

# колонка с покрытием всегда вторая
cover_col = c.columns[1]

# маппинг к порядку из CheckM2
c_map = c.set_index(c.columns[0])[cover_col]
c_ordered = [c_map.get(x, float('nan')) for x in order]
q_sub['Coverage'] = c_ordered

# --- 3. GTDB-Tk ---
g = pd.read_csv(sys.argv[3], sep="\t") #gtdbtk.bac120.summary.tsv

g = g.set_index('user_genome')

# classification → последний таксон
def last_taxon(x):
    if x in ("Unclassified", "Unclassified Bacteria"):
        return x

    parts = []
    for level in x.split(';'):
        if '__' not in level:
            continue
        _, name = level.split('__', 1)
        name = name.strip()
        if name:
            parts.append(name)

    return parts[-1] if parts else ""



g_last = g['classification'].apply(last_taxon)
g_ani = g['closest_genome_ani']

g_last_ordered = [g_last.get(x, "") for x in order]
g_ani_ordered = [g_ani.get(x, float('nan')) for x in order]

q_sub['Taxon'] = g_last_ordered
q_sub['Closest_ANI'] = g_ani_ordered

# --- 4. Abricate summary ---
vf = pd.read_csv(sys.argv[4], sep="\t") #vf_summary.tsv

# очистка имени бина
def clean_bin(x):
    # SemiBin_100.fa_abr_hits.tsv → SemiBin_100
    x = re.sub(r'\.fa.*$', '', x)
    return x

vf['clean'] = vf[vf.columns[0]].apply(clean_bin)

vf = vf.set_index('clean')

# все колонки кроме первой исходной
cols_vf = vf.columns.tolist()
cols_vf.remove(vf.columns[0])  # remove original filename-column
#cols_vf.remove('clean')

# добавить по порядку
for col in cols_vf:
    q_sub[col] = [vf.at[x, col] if x in vf.index else float('nan') for x in order]
# --- финальный вывод ---
q_sub.to_csv("genomes_report.csv", sep=';', index=False)



vlist = pd.read_csv(sys.argv[5], sep='\t')

compl_map = dict(zip(q_sub['Name'], q_sub['Completeness']))
contam_map = dict(zip(q_sub['Name'], q_sub['Contamination']))
taxon_map = dict(zip(q_sub['Name'], q_sub['Taxon']))

vlist['BIN'] = vlist['#FILE'].apply(clean_bin)
vlist['COMPL'] = vlist['BIN'].map(compl_map)
vlist['CONTAM'] = vlist['BIN'].map(contam_map)
vlist['TAXON'] = vlist['BIN'].map(taxon_map)

blist = pd.read_csv(sys.argv[6], sep='\t')

#blist['clean'] = blist['qseqid'].apply(clean_bin)

# qseqid уже = region_id
blist['region_id'] = blist['qseqid']

blist['blast_hit'] = (
    blist['sscinames'].astype(str)
    + ": "
    + blist['pident'].astype(str)
    + "/"
    + blist['qcovs'].astype(str)
)

blast_grouped = (
    blist
    .groupby('region_id')['blast_hit']
    .apply(lambda x: "\n".join(x))
)

# --- формируем region_id в vlist ---
vlist['region_id'] = (
    vlist['#FILE'].astype(str)
    + "_"
    + vlist['SEQUENCE'].astype(str)
    + "_"
    + vlist['START'].astype(str)
    + "-"
    + vlist['END'].astype(str)
)

vlist['BLAST (id/cov)'] = [
    blast_grouped.get(x, "") for x in vlist['region_id']
]

# переставляем колонки так, чтобы COMPL и CONTAM были после #FILE
cols = vlist.columns.tolist()
cols.remove('COMPL')
cols.remove('CONTAM')
cols.remove('TAXON')
cols.remove('BIN')
cols.remove('BLAST (id/cov)')
cols.remove('COVERAGE_MAP')
cols.remove('DATABASE')
cols.remove('RESISTANCE')
cols.remove('region_id')
insert_pos = cols.index('#FILE') + 1
cols = cols[:insert_pos] + ['COMPL', 'CONTAM', 'TAXON', 'BLAST (id/cov)'] + cols[insert_pos:]
vlist = vlist[cols]

# --- 6. Сохраняем итог ---
vlist.to_csv("virulence_list.tsv", sep='\t', index=False)


# --- 7. PLOTS ----

def parse_mags(set):
    df = set
    df = df[["Name", "Taxon", 'Completeness', 'Contamination', 'Genome_Size', 'NUM_FOUND']].copy()
    df['Taxon']=df['Taxon'].str.replace('_[A,B,C,D,E]', '', regex=True)
    return df


df = parse_mags(q_sub)
# масштабирование площади точек
size_scale = 1e-3
sizes = df["Genome_Size"] * size_scale

# качество: высокая полнота + низкая контаминация
quality = df["Completeness"]+10 - df["Contamination"]*1.5

cmap = LinearSegmentedColormap.from_list(
    "qualmap",
    ["#ffcccc", "#fff2b2", "#ccffcc"]  # светло-красный → светло-желтый → светло-зеленый
)

norm = Normalize(vmin=quality.min(), vmax=quality.max())

fig, axes = plt.subplots(1, 2, figsize=(18, 8), gridspec_kw={"width_ratios": [2, 1]})

ax = axes[0]

# базовый слой
scatter = ax.scatter(
    df["Completeness"],
    df["Contamination"],
    s=sizes,
    c=quality,
    cmap=cmap,
    norm=norm,
    edgecolors="black",
    linewidth=0.3,
    alpha=0.8
)

# обводка красным при наличии факторов
has_vf = df["NUM_FOUND"] > 0
text_vf = (df["NUM_FOUND"] > 0) & ((df['Completeness']<80) | (df['Contamination']>5))
ax.scatter(
    df.loc[has_vf, "Completeness"],
    df.loc[has_vf, "Contamination"],
    s=df.loc[has_vf, "Genome_Size"] * size_scale,
    facecolors="none",
    edgecolors="#FF6161FF",
    linewidth=1
)

#подписи (аналог ggrepel)
texts = []
for _, row in df.loc[text_vf].iterrows():
    texts.append(
        ax.text(
            row["Completeness"],
            row["Contamination"],
            row["Taxon"],
            fontsize=16
        )
    )

adjust_text(texts, arrowprops=dict(arrowstyle='- >', lw=1.2), ax=ax, 
            max_move=80, force_text=50, force_points=80, iter_lim=900, only_move={'text':'y'})
            
has_vf = df["NUM_FOUND"] > 1
ax.scatter(
    df.loc[has_vf, "Completeness"],
    df.loc[has_vf, "Contamination"],
    s=df.loc[has_vf, "Genome_Size"] * size_scale,
    facecolors="none",
    edgecolors="#CC0303FF",
    linewidth=1.8
)


ax.grid(True, linestyle='--', linewidth=0.4, alpha=0.5)

ax.set_xlabel("Completeness (%)", fontsize=14)
ax.set_ylabel("Contamination (%)", fontsize=14)
ax.set_xlim(0, 105)
ax.set_ylim(-6, df["Contamination"].max() + 10)
ax.axhline(5,linestyle='--', color='grey')
ax.axvline(80, linestyle='--', color='grey')


#-----------------------------------------------------------------------------

ax2 = axes[1]
df2 = df[(df["Completeness"]>=80) & (df["Contamination"]<=5)]
sizes = df2["Genome_Size"] * size_scale


scatter = ax2.scatter(
    df2["Completeness"],
    df2["Contamination"],
    s=sizes,
    c="#ccffca",
    edgecolors="black",
    linewidth=0.1,
    alpha=0.8
)

# обводка красным при наличии факторов
has_vf = df2["NUM_FOUND"] > 0
ax2.scatter(
    df2.loc[has_vf, "Completeness"],
    df2.loc[has_vf, "Contamination"],
    s=df2.loc[has_vf, "Genome_Size"] * size_scale,
    c="#FF6161FF",
    linewidth=1,
    alpha=0.4
)

#подписи (аналог ggrepel)
texts = []
for _, row in df2.loc[has_vf].iterrows():
    texts.append(
        ax2.text(
            row["Completeness"],
            row["Contamination"],
            row["Taxon"],
            fontsize=18
        )
    )

adjust_text(texts, arrowprops=dict(arrowstyle='- >', lw=1.2, shrinkA=10, shrinkB=5))

# выбираем 3–5 репрезентативных размеров генома
sizes_for_legend = np.linspace(df['Genome_Size'].min(),
                               df['Genome_Size'].max(),
                               5)

def genome_to_area(g):
    # подставь ТУ ЖЕ формулу, что используется в твоём scatter (... s=...)
    return g * size_scale
legend_handles = []
legend_labels = []
for g in sizes_for_legend:
    s = genome_to_area(g)
    h = ax.scatter([], [], s=s,
                   facecolors="lightgray",
                   edgecolors="black",
                   linewidths=1)
    legend_handles.append(h)
    legend_labels.append(f"{int(g)/1000000:.2f}")

ax.legend(
    legend_handles,
    legend_labels,
    title="Genome size (Gbases)",
    title_fontsize=16,
    loc="best",
    frameon=False,
    scatterpoints=1,
    labelspacing=3.6,
    fontsize=14,
    reverse=True,
    handletextpad=3.3,
    borderpad=1.1
)


ax2.grid(True, linestyle='--', linewidth=0.4, alpha=0.5)

ax2.set_xlim(75, 105)
ax2.set_ylim(-0.5, 5)
ax2.xaxis.set_major_locator(mticker.MultipleLocator(10))

plt.tight_layout()
plt.savefig('bins_plot.png', dpi=300)
